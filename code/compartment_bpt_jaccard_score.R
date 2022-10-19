library(mgcv)
library(vroom)
library(Matrix)
library(viridis)
library(data.tree)
library(GenomicRanges)
library(furrr)
library(igraph)
library(seriation)
library(tidyverse)

#--------------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
##Utils. Fn

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(base::load(tmp_file))
  tmp_obj<-names(mget(base::load(tmp_file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,res_num){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat,cluster = 10)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}

full_f_mat<-function(cl_mat,res,var){
  
  range_5kb<-range(as.numeric(unique(c(cl_mat$X1,cl_mat$X2))))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[var]]),symmetric = T)
  return(chr_mat)
}

produce_bin_bpt_map_fn<-function(chr_spec_res,tmp_bins,tmp_res){
  
  chr_bin_GRange<-GRanges(seqnames=chromo,
                          ranges = IRanges(start=tmp_bins,
                                           end=tmp_bins + res_num[tmp_res] -1
                          ))
  plan(multisession,workers=5)
  chr_cl_tbl<-tibble(chr=chromo,cl=names(chr_spec_res$cl_member),bins=chr_spec_res$cl_member) %>% 
    mutate(res=str_split_fixed(cl,"_",2)[,1]) %>% 
    mutate(ends=pmap(list(bins,res),function(bins,res){
      as.numeric(bins) + res_num[res] - 1
    })) %>% 
    
    mutate(GRange=future_pmap(list(chr,bins,ends),function(chr,bins,ends){
      GRanges(seqnames=chr,
              ranges = IRanges(start=as.numeric(bins),
                               end=ends
              ))
    }))
  plan(sequential)
  
  cl_GRange_l<-GRangesList(chr_cl_tbl$GRange)
  bin_cl_inter_tbl<-findOverlaps(chr_bin_GRange,cl_GRange_l) %>% 
    as_tibble
  
  cl_bin_intersect<-pintersect(chr_bin_GRange[bin_cl_inter_tbl$queryHits],cl_GRange_l[bin_cl_inter_tbl$subjectHits])
  
  plan(multisession,workers=5)
  intersect_size<-future_map_int(1:length(cl_bin_intersect),function(i){
    sum(IRanges::width(cl_bin_intersect[[i]]))
  })
  plan(sequential)
  
  bin_cl_inter_tbl<-bin_cl_inter_tbl %>% 
    mutate(inter.size=intersect_size) %>% 
    mutate(bin=tmp_bins[queryHits],
           cl=chr_cl_tbl$cl[subjectHits]) %>% 
    mutate(cl.lvl=node_lvl[cl])
  
  
  max_lvl_tbl<-bin_cl_inter_tbl %>% 
    #  filter(inter.size>=res_num[tmp_res]) %>% 
    group_by(bin) %>% 
    slice_max(inter.size) %>% 
    slice_max(cl.lvl)
  return(max_lvl_tbl %>% 
           dplyr::select(bin,cl,inter.size,cl.lvl))
}
produce_bpt_split_tbl_fn<-function(max_lvl_tbl,node_ancestor,top_parent){
  bin_to_cl_vec<-max_lvl_tbl$cl
  names(bin_to_cl_vec)<-max_lvl_tbl$bin
  
  bin_bpt_tbl<-tibble(bin=names(bin_to_cl_vec)) %>% 
    mutate(cl=bin_to_cl_vec[bin])
  
  target_cl<-unique(c(bin_bpt_tbl$cl))  
  bin_top_parent<-unlist(lapply(target_cl,function(x){
    top_parent[which(top_parent %in% c(x,node_ancestor[[x]]))]
  }))
  names(bin_top_parent)<-target_cl
  bin_bpt_tbl<-bin_bpt_tbl %>% 
    mutate(top.parent=bin_top_parent[cl])
  top_parent<-unique(bin_top_parent)
  bin_bpt_tbl<-bin_bpt_tbl %>% 
    mutate(split=ifelse(top.parent == top_parent[1],1,-1))
  return(bin_bpt_tbl)
}

produce_comp_split_tbl_fn<-function(chr_dat){
  chr_mat<-full_f_mat(chr_dat,res_num[tmp_res],"zscore")
  
  range_bin<-range(unique(c(chr_dat$X1,chr_dat$X2)))
  f_chr_bin<-seq(range_bin[1],range_bin[2],by=res_num[tmp_res])
  dimnames(chr_mat)<-list(f_chr_bin,f_chr_bin)
  
  empty_rows<-which(apply(chr_mat,1,function(x)all(x==0)))
  empty_cols<-which(apply(chr_mat,2,function(x)all(x==0)))
  if(length(empty_rows)>0){
    ok_chr_mat<-chr_mat[-empty_rows,]
    ok_chr_mat<-ok_chr_mat[,-empty_cols]
  } else{
    ok_chr_mat<-chr_mat
  }
  
  cor_mat<-cor(as.matrix(ok_chr_mat))
  
  eig_vec<-svd(cor_mat)$v[,1]
  names(eig_vec)<-colnames(ok_chr_mat)
  
  comp_tbl<-tibble(X1=names(eig_vec)) %>% 
    mutate(bin.eig=eig_vec[X1]) %>% 
    mutate(comp=sign(bin.eig))
  return(comp_tbl)
}
produce_bin_GRange_fn<-function(bin_comp_tbl){
  
  return(GRanges(seqnames=bin_comp_tbl$chr,
                          ranges = IRanges(start=as.numeric(bin_comp_tbl$bin),
                                           end=as.numeric(bin_comp_tbl$bin) + res_num[bin_comp_tbl$res] -1
                          )))
  
}

Jaccard_fn<-function(Grange_A,Grange_B){
  sum(width(GenomicRanges::intersect(Grange_A,Grange_B)))/sum(width(GenomicRanges::union(Grange_A,Grange_B)))
  
}

jaccard_comp_compare_fn<-function(bin_comp_tbl,bin_bpt_tbl,chromo,tmp_res){
  
  comp_GRange_l<-list(
    comp.A=IRanges::reduce(produce_bin_GRange_fn(bin_comp_tbl %>% 
                                                   filter(comp==1)%>%
                                                   dplyr::rename(bin=X1) %>% 
                                                   mutate(chr=chromo,
                                                          res=tmp_res))),
    comp.B=IRanges::reduce(produce_bin_GRange_fn(bin_comp_tbl %>% 
                                                   filter(comp== -1)%>%
                                                   dplyr::rename(bin=X1) %>% 
                                                   mutate(chr=chromo,
                                                          res=tmp_res)))
    
  )
  
  bpt_GRange_l<-list(
    bpt.A=IRanges::reduce(produce_bin_GRange_fn(bin_bpt_tbl %>% 
                                                  filter(split==1)%>%
                                                  mutate(chr=chromo,
                                                         res=tmp_res))),
    bpt.B=IRanges::reduce(produce_bin_GRange_fn(bin_bpt_tbl %>% 
                                                  filter(split== -1)%>%
                                                  mutate(chr=chromo,
                                                         res=tmp_res)))
  )
  
  obs_jaccard<-expand_grid(comp=c("comp.A","comp.B"),bpt=c("bpt.A","bpt.B")) %>% 
    mutate(comp.GRange=comp_GRange_l[comp],
           bpt.GRange=bpt_GRange_l[bpt]) %>% 
    mutate(jaccard=pmap_dbl(list(comp.GRange,bpt.GRange),function(Grange_A,Grange_B){
      Jaccard_fn(Grange_A,Grange_B)
    })) %>% 
    dplyr::select(-c(comp.GRange, bpt.GRange))  %>%
    group_by(comp) %>% slice_max(jaccard) %>%  ungroup()

#  %>%vgroup_by(comp) %>% slice_max(jaccard) %>%  ungroup()
#  %>% summarise(jaccard.score=sum(jaccard)) %>% unlist
  
  return(obs_jaccard)
}
#-----------------------------------------
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
HiC_spec_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
chromo<-"chr7" # chr11 + chr7
chr_set<-str_split_fixed(grep("^chr",list.files(HiC_spec_folder),value=T),"_",2)[,1]

rn_jaccard_res_l<-vector('list',length(chr_set))
obs_jaccard_res_l<-vector('list',length(chr_set))
names(rn_jaccard_res_l)<-chr_set
names(obs_jaccard_res_l)<-chr_set


for(chromo in chr_set){
  message(chromo,": Compute cluster and compartment objects")
  
  chr_spec_res<-get_tbl_in_fn(paste0(HiC_spec_folder,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_lvl<-chr_bpt$Get('level')
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  top_parent<-names(which(node_lvl==2))
  tmp_res<-unique(str_split_fixed(top_parent,"_",2)[,1])
  chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
  tmp_bins<-unique(c(chr_dat$X1,chr_dat$X2))
  
  
  max_lvl_tbl<-produce_bin_bpt_map_fn(chr_spec_res,tmp_bins,tmp_res)
  bin_bpt_tbl<-produce_bpt_split_tbl_fn(max_lvl_tbl,node_ancestor,top_parent)
  
  #-----------------------------------------
  message(chromo,": Obs Jaccard")
  
  bin_comp_tbl<-produce_comp_split_tbl_fn(chr_dat)
  #-----------------------------------------
  ## Compute Jaccard similarity index
  obs_jaccard<-jaccard_comp_compare_fn(bin_comp_tbl,bin_bpt_tbl,chromo,tmp_res)
  #-----------------------------------------
  # Compare to random bpt splits
  message(chromo,": Random Jaccard")
  plan(multisession,workers=5)
  rn_jaccard_score<-future_map(1:20,function(i){
    rn_bin_bpt_tbl<-bin_bpt_tbl %>% 
      mutate(split=sample(split))
    return(jaccard_comp_compare_fn(bin_comp_tbl,rn_bin_bpt_tbl,chromo,tmp_res) %>% mutate(set=paste0("rn.",i)))
  })
  plan(sequential)
  rn_jaccard_res_l[[chromo]]<-do.call(bind_rows,rn_jaccard_score) %>% mutate(chr=chromo)
  obs_jaccard_res_l[[chromo]]<-obs_jaccard %>% mutate(chr=chromo,set="obs")
}
do.call(bind_rows,obs_jaccard_res_l)

gg_comp<-do.call(bind_rows,rn_jaccard_res_l) %>% 
  mutate(chr=fct_relevel(chr,paste0("chr",1:22)),
         set1=str_split_fixed(set,"\\.",2)[,1]) %>% 
  ggplot(.,aes(chr,jaccard,color=comp))+
  geom_violin()+
  geom_point(data=do.call(bind_rows,obs_jaccard_res_l) %>% 
               mutate(set1="obs"))+
  theme_minimal()+
  scale_color_brewer(palette="Set2")+
  facet_grid(comp~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly62/img/comp_jaccard.png",gg_comp)
tibble(score=rn_jaccard_score,set="random") %>%
  ggplot(.,aes(chr,score,fill=set))+
  geom_violin()+
  geom_point(data = tibble(x="chr",score=obs_jaccard,set="observed"),size=5,color="red")


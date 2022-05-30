library(mgcv)
library(vroom)
library(Matrix)
library(viridis)
library(data.tree)
library(GenomicRanges)
library(furrr)
library(igraph)
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

#-----------------------------------------
HiC_dat_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/"
HiC_spec_folder<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
out_file<-"~/data_transfer/GM12878_eig_bpt_dist_100kb.tsv"

chr_set<-str_split_fixed(grep("^chr",list.files(HiC_spec_folder),value=T),"_",2)[,1]

tmp_res<-"100kb"

chr_res_l<-vector("list",length(chr_set))
names(chr_res_l)<-chr_set

for(chromo in chr_set){
  
  message("Build BPT for: ",chromo)
  
  chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
  chr_spec_res<-get_tbl_in_fn(paste0(HiC_spec_folder,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_lvl<-chr_bpt$Get('level')
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  
  message("Map bins to cluster for: ",chromo)
  
  tmp_bins<-unique(c(chr_dat$X1,chr_dat$X2))
  
  chr_bin_GRange<-GRanges(seqnames=chromo,
                          ranges = IRanges(start=tmp_bins,
                                           end=tmp_bins + res_num[tmp_res] -1
                          ))
  plan(multisession,workers=15)
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
  
  plan(multisession,workers=10)
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
    filter(inter.size==res_num[tmp_res]) %>% 
    group_by(bin) %>% 
    slice_max(cl.lvl)
  
  bin_to_cl_vec<-max_lvl_tbl$cl
  names(bin_to_cl_vec)<-max_lvl_tbl$bin
  
  bin_node<-unique(bin_to_cl_vec)
  bin_cl_set<-unique(c(bin_node,unique(unlist(node_ancestor[bin_node])))) 
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  Prune(chr_bpt, function(x) x$name %in% bin_cl_set)
  
  g_bpt<-as.igraph.Node(chr_bpt,directed = T,direction = 'climb')
  
  tmp_d<-distances(g_bpt,bin_node,bin_node)
  
  message("Compute BPT distances for: ",chromo)
  
  bin_inter_tbl<-expand_grid(X1=names(bin_to_cl_vec),X2=names(bin_to_cl_vec)) %>% 
    mutate(cl.A=bin_to_cl_vec[X1],
           cl.B=bin_to_cl_vec[X2])
  
  bin_inter_tbl<-bin_inter_tbl %>% 
    mutate(bpt.d=tmp_d[as.matrix(bin_inter_tbl[,3:4])])
  
  message("Compute eigen distances for: ",chromo)
  
  chr_mat<-full_f_mat(chr_dat,res_num[tmp_res],"zscore")
  
  range_bin<-range(unique(c(chr_dat$X1,chr_dat$X2)))
  f_chr_bin<-seq(range_bin[1],range_bin[2],by=res_num[tmp_res])
  dimnames(chr_mat)<-list(f_chr_bin,f_chr_bin)
  
  empty_rows<-which(apply(chr_mat,1,function(x)all(x==0)))
  empty_cols<-which(apply(chr_mat,2,function(x)all(x==0)))
  
  ok_chr_mat<-chr_mat[-empty_rows,]
  ok_chr_mat<-ok_chr_mat[,-empty_cols]
  
  cor_mat<-cor(as.matrix(ok_chr_mat))
  
  eig_vec<-svd(cor_mat)$v[,1]
  names(eig_vec)<-colnames(ok_chr_mat)
  
  message("Output summary table for: ",chromo)
  
  eig_dist_tbl<- expand_grid(X1=names(eig_vec),X2=names(eig_vec)) %>% 
    mutate(binA.eig=eig_vec[X1],
           binB.eig=eig_vec[X2]) %>% 
    mutate(eig.dist=(abs(binA.eig-binB.eig)),
           same.comp=ifelse(sign(binA.eig)*sign(binB.eig)<0,"diff","same"))
  
  chr_res_l[[chromo]]<-bin_inter_tbl %>% 
    inner_join(.,eig_dist_tbl) %>% 
    filter(!(is.na(eig.dist))) %>%
    mutate(gdist=abs(as.numeric(X1)-as.numeric(X2)),
           chr=chromo)
  
}

out_tbl<-do.call(bind_rows,chr_res_l)
write_tsv(out_tbl,file = out_file)


gg_tmp<-out_tbl %>% 
  filter(!(is.na(eig.dist))) %>%
  ggplot(.,aes(as.factor(bpt.d),eig.dist))+
  #  geom_smooth()+
  #  geom_point(alpha=0.01)#+ facet_grid(.~same.comp,scales="free")
  geom_boxplot(outlier.size=0.1) #+ geom_smooth(mapping = aes(bpt.d,eig.dist),se=F) + facet_grid(.~same.comp,scales="free")

gg_tmp<-out_tbl %>% 
  filter(!(is.na(eig.dist))) %>%
  ggplot(.,aes(bpt.d,same.comp))+
  #  geom_smooth()+
  #  geom_point(alpha=0.01)#+ facet_grid(.~same.comp,scales="free")
  geom_boxplot(outlier.size=0.1) #+ geom_smooth(mapping = aes(bpt.d,eig.dist),se=F) + facet_grid(.~same.comp,scales="free")

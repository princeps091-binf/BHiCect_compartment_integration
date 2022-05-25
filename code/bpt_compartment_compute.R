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

full_bpt_mat<-function(cl_mat,res,var){
  
  range_5kb<-range(as.numeric(unique(c(cl_mat$X1,cl_mat$X2))))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[var]]))
  return(chr_mat)
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
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
HiC_spec_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"



chromo<-"chr19"
tmp_res<-"100kb"

chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
chr_spec_res<-get_tbl_in_fn(paste0(HiC_spec_folder,chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_lvl<-chr_bpt$Get('level')
node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
node_ancestor<-lapply(node_ancestor,'[',-1)

tmp_bins<-unique(c(chr_dat$X1,chr_dat$X2))

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
intersect_size<-unlist(lapply(lapply(pintersect(chr_bin_GRange[bin_cl_inter_tbl$queryHits],cl_GRange_l[bin_cl_inter_tbl$subjectHits]),width),sum))
union_size<-unlist(lapply(lapply(lapply(punion(chr_bin_GRange[bin_cl_inter_tbl$queryHits],cl_GRange_l[bin_cl_inter_tbl$subjectHits]),IRanges::reduce),width),sum))

bin_cl_inter_tbl<-bin_cl_inter_tbl %>% 
  mutate(inter.size=intersect_size,
         union.size=union_size) %>% 
  mutate(jaccard=inter.size/union.size) %>% 
  mutate(bin=tmp_bins[queryHits],
         cl=chr_cl_tbl$cl[subjectHits]) %>% 
  mutate(cl.lvl=node_lvl[cl])

max_jaccard_tbl<-bin_cl_inter_tbl %>% 
  group_by(bin) %>% 
  slice_max(jaccard) %>% 
  arrange(jaccard)

max_lvl_tbl<-bin_cl_inter_tbl %>% 
  filter(inter.size==res_num[tmp_res]) %>% 
  group_by(bin) %>% 
  slice_max(cl.lvl) %>% 
  arrange(jaccard)

bin_to_cl_vec<-max_lvl_tbl$cl
names(bin_to_cl_vec)<-max_lvl_tbl$bin

bin_node<-unique(bin_to_cl_vec)
bin_cl_set<-unique(c(bin_node,unique(unlist(node_ancestor[bin_node])))) 
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
Prune(chr_bpt, function(x) x$name %in% bin_cl_set)

g_bpt<-as.igraph.Node(chr_bpt,directed = T,direction = 'climb')

tmp_d<-distances(g_bpt,bin_node,bin_node)

bin_inter_tbl<-expand_grid(X1=names(bin_to_cl_vec),X2=names(bin_to_cl_vec)) %>% 
  mutate(cl.A=bin_to_cl_vec[X1],
         cl.B=bin_to_cl_vec[X2])

bin_inter_tbl<-bin_inter_tbl %>% 
  mutate(bpt.d=tmp_d[as.matrix(bin_inter_tbl[,3:4])])
chr_bpt_mat<-full_bpt_mat(bin_inter_tbl,res_num[tmp_res],"bpt.d")
image(as.matrix(chr_bpt_mat),col=viridis(100))
d<-as.dist(chr_bpt_mat)
o <- seriate(d,method = "HC")
image(as.matrix(chr_bpt_mat)[get_order(o),get_order(o)],col=viridis(100))

chr_mat<-full_f_mat(chr_dat,res_num[tmp_res],"zscore")
image(as.matrix(chr_mat),col=viridis(100))
image(as.matrix(chr_mat)[get_order(o),get_order(o)],col=viridis(100))

cor_mat<-cor(as.matrix(chr_mat))

image(cor_mat[get_order(o),get_order(o)],col=viridis(100))

png(paste0('./img/',chromo,"_cor_",tmp_res,"_bpt_mat",'.png'), width =40,height = 43,units = 'mm',type='cairo',res=5000)
par(mar = c(0, 0, 0,0))
plot.new()
image(cor_mat[get_order(o),get_order(o)],col=viridis(100))
dev.off()
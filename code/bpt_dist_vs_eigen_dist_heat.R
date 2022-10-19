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

full_f_mat<-function(cl_mat,res,var,symm){
  
  range_5kb<-range(as.numeric(unique(c(cl_mat$X1,cl_mat$X2))))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  
  #chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=sqrt(-log10(cl_mat$pois.pval)),symmetric = T)
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[var]]),symmetric = symm)
  dimnames(chr_mat)<-list(bin_5kb,bin_5kb)
  return(chr_mat)
}

produce_seriate_idx<-function(dist_mat){
  d<-as.dist(dist_mat)
  o <- seriate(d)
  return(get_order(o))
  
}

produce_heat_png_fn<-function(out_file,heat_mat){
  
  png(out_file, width =40,height = 40,units = 'mm',type='cairo',res=5000)
  par(mar = c(0, 0, 0,0))
  plot.new()
  image(heat_mat,col=viridis(100),useRaster=T)
  dev.off()
  
}
#-----------------------------------------
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
HiC_spec_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"



chromo<-"chr22"

chr_spec_res<-get_tbl_in_fn(paste0(HiC_spec_folder,chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_lvl<-chr_bpt$Get('level')
tmp_res<-unique(str_split_fixed(names(which(node_lvl==2)),"_",2)[,1])
chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)

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
# Make parallel
intersect_cpu_l<-as.list(width(pintersect(chr_bin_GRange[bin_cl_inter_tbl$queryHits],cl_GRange_l[bin_cl_inter_tbl$subjectHits])))
intersect_size<-map_int(intersect_cpu_l,sum)

bin_cl_inter_tbl<-bin_cl_inter_tbl %>% 
  mutate(inter.size=intersect_size) %>% 
  mutate(bin=tmp_bins[queryHits],
         cl=chr_cl_tbl$cl[subjectHits]) %>% 
  mutate(cl.lvl=node_lvl[cl])

max_lvl_tbl<-bin_cl_inter_tbl %>% 
  filter(inter.size>=res_num[tmp_res]) %>% 
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

bin_inter_tbl<-expand_grid(X1=names(bin_to_cl_vec),X2=names(bin_to_cl_vec)) %>% 
  mutate(cl.A=bin_to_cl_vec[X1],
         cl.B=bin_to_cl_vec[X2])

bin_inter_tbl<-bin_inter_tbl %>% 
  mutate(bpt.d=tmp_d[as.matrix(bin_inter_tbl[,3:4])])


chr_mat<-full_f_mat(chr_dat,res_num[tmp_res],"zscore",T)

range_bin<-range(unique(c(chr_dat$X1,chr_dat$X2)))
f_chr_bin<-seq(range_bin[1],range_bin[2],by=res_num[tmp_res])
dimnames(chr_mat)<-list(f_chr_bin,f_chr_bin)

empty_rows<-which(apply(chr_mat,1,function(x)all(x==0)))
empty_cols<-which(apply(chr_mat,2,function(x)all(x==0)))
if(length(empty_rows)>1){
  ok_chr_mat<-chr_mat[-empty_rows,]
  ok_chr_mat<-ok_chr_mat[,-empty_cols]
  
} else{
  ok_chr_mat<-chr_mat
}

cor_mat<-cor(as.matrix(ok_chr_mat))

eig_vec<-svd(cor_mat)$v[,1]
names(eig_vec)<-colnames(ok_chr_mat)

eig_dist_tbl<-expand_grid(X1=names(eig_vec),X2=names(eig_vec)) %>% 
  mutate(binA.eig=eig_vec[X1],
         binB.eig=eig_vec[X2]) %>% 
  mutate(eig.dist=(abs(binA.eig-binB.eig)),
         same.comp=ifelse(sign(binA.eig)*sign(binB.eig)<0,"diff","same"))
eig_mat<-full_f_mat(eig_dist_tbl,res_num[tmp_res],"eig.dist",F)
bpt_mat<-full_f_mat(bin_inter_tbl,res_num[tmp_res],"bpt.d",F)

eig_idx<-produce_seriate_idx(eig_mat)
bpt_idx<-produce_seriate_idx(bpt_mat)

image(as.matrix(bpt_mat)[bpt_idx,bpt_idx],col=viridis(100))
image(as.matrix(eig_mat)[eig_idx,eig_idx],col=viridis(100))
# out viz
out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_bpt_mat_reorder.png")
produce_heat_png_fn(out_file,as.matrix(bpt_mat)[bpt_idx,bpt_idx])
out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_eig_mat_reorder.png")
produce_heat_png_fn(out_file,as.matrix(eig_mat)[eig_idx,eig_idx])

cor_bins_combo<-expand_grid(X1=unique(unlist(dimnames(cor_mat))),X2=unique(unlist(dimnames(cor_mat))))
cor_mat_gap<-matrix(0,nrow=nrow(chr_mat),ncol=ncol(chr_mat),dimnames = dimnames(chr_mat))
cor_mat_gap[as.matrix(cor_bins_combo)]<-cor_mat[as.matrix(cor_bins_combo)]

image(cor_mat_gap,col=viridis(100))
image(cor_mat_gap[eig_idx,eig_idx],col=viridis(100))
image(cor_mat_gap[bpt_idx,bpt_idx],col=viridis(100))


chr_mat_raw<-full_f_mat(chr_dat,res_num[tmp_res],"lw",T)
image(as.matrix(chr_mat_raw),col=viridis(100))
image(as.matrix(chr_mat_raw)[eig_idx,eig_idx],col=viridis(100))
image(as.matrix(chr_mat_raw)[bpt_idx,bpt_idx],col=viridis(100))

out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_cor_mat_raw.png")
produce_heat_png_fn(out_file,cor_mat_gap)
out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_raw_mat_raw.png")
produce_heat_png_fn(out_file,as.matrix(chr_mat_raw))

out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_cor_mat_eigen_order.png")
produce_heat_png_fn(out_file,cor_mat_gap[eig_idx,eig_idx])
out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_cor_mat_bpt_order.png")
produce_heat_png_fn(out_file,cor_mat_gap[bpt_idx,bpt_idx])

out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_raw_mat_eigen_order.png")
produce_heat_png_fn(out_file,as.matrix(chr_mat_raw)[eig_idx,eig_idx])
out_file<-paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_raw_mat_bpt_order.png")
produce_heat_png_fn(out_file,as.matrix(chr_mat_raw)[bpt_idx,bpt_idx])



gg_comp<-tibble(
raw.val=as.matrix(chr_mat_raw)[as.matrix(expand_grid(ego=colnames(chr_mat_raw),alter=colnames(chr_mat_raw)))],
cor.val=as.matrix(cor_mat_gap)[as.matrix(expand_grid(ego=colnames(chr_mat_raw),alter=colnames(chr_mat_raw)))]
) %>% 
  ggplot(.,aes(raw.val,cor.val))+
  geom_point(size=0.1,alpha=0.4)
ggsave(paste0("~/Documents/multires_bhicect/weeklies/weekly62/img/",chromo,"_raw_cor_scatter.png"),gg_comp)

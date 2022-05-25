library(mgcv)
library(vroom)
library(Matrix)
library(viridis)
library(tidyverse)

#--------------------------------
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
##Utils. Fn

get_tbl_in_fn<-function(tmp_file){
  out_tbl<-get(load(tmp_file))
  tmp_obj<-names(mget(load(tmp_file)))
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
  
  range_5kb<-range(unique(c(cl_mat$X1,cl_mat$X2)))
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

chromo<-"chr1"
tmp_res<-"100kb"

chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)

chr_mat<-full_f_mat(chr_dat,res_num[tmp_res],"zscore")
image(as.matrix(chr_mat),col=viridis(100))
empty_rows<-which(apply(chr_mat,1,function(x)all(x==0)))
empty_cols<-which(apply(chr_mat,2,function(x)all(x==0)))
ok_chr_mat<-chr_mat[-empty_rows,]
ok_chr_mat<-ok_chr_mat[,-empty_cols]
image(as.matrix(ok_chr_mat),col=viridis(100))

cor_mat<-cor(as.matrix(ok_chr_mat))
image(cor_mat,col=viridis(100))

eigen_idx<-sort(svd(cor_mat)$v[,1],index.return=T)$ix
image(as.matrix(ok_chr_mat)[eigen_idx,eigen_idx],col=viridis(100))
image(cor_mat[eigen_idx,eigen_idx],col=viridis(100))

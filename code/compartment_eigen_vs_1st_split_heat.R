library(mgcv)
library(vroom)
library(Matrix)
library(caret)
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

lp_fn<-function(x){
  
  Dinv=Diagonal(nrow(x),1/Matrix::rowSums(x))
  
  lp_chr1=Diagonal(nrow(x),1)-Dinv %*% x
  dimnames(lp_chr1)<-dimnames(x)
  if(dim(lp_chr1)[1] > 10000){
    return(eigs_sym(lp_chr1,k=2,sigma = 0, which='LM',maxitr=10000))
  }
  else{
    temp<-eigen(lp_chr1)
    return(list(tibble(bin=rownames(lp_chr1),fielder=temp[['vectors']][,length(temp$values)-1],indicator=temp[['vectors']][,length(temp$values)]),values=temp[['values']][c(length(temp$values)-1,length(temp$values))]))
    
    
  }
}

produce_fiedler_fn<-function(chr_dat){
  g_chr1<-graph_from_data_frame(chr_dat %>% dplyr::select(X1,X2,weight),directed = F)
  #eleminate self loop 
  g_chr1<-delete.edges(g_chr1,E(g_chr1)[which(which_loop(g_chr1))])
  chr_mat<-get.adjacency(g_chr1,type='both',attr='weight')
  diag(chr_mat)<-0
  if(any(colSums(chr_mat)==0)){
    out<-which(colSums(chr_mat)==0)
    chr_mat<-chr_mat[-out,]
    chr_mat<-chr_mat[,-out]
  }
  return(lp_fn(chr_mat)[[1]])
  
}
#-----------------------------------------
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
HiC_spec_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

chromo<-"chr11"
chr_spec_res<-get_tbl_in_fn(paste0(HiC_spec_folder,chromo,"_spec_res.Rda"))
chr_bpt<-FromListSimple(chr_spec_res$part_tree)
node_lvl<-chr_bpt$Get('level')
tmp_res<-unique(str_split_fixed(names(which(node_lvl==2)),"_",2)[,1])
chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
preprocessParams <- BoxCoxTrans(chr_dat$X3,na.rm = T)

chr_dat<-chr_dat %>% 
  mutate(weight=predict(preprocessParams, .$X3)) %>% 
  mutate(weight=weight+(1-min(.$weight,na.rm = T)))

fiedler_tbl<-produce_fiedler_fn(chr_dat)


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
eig_bin_tbl<-tibble(bin=names(eig_vec),eigen=eig_vec)
fiedler_tbl %>% 
  inner_join(.,eig_bin_tbl) %>% 
  ggplot(.,aes(fielder,eigen))+
  geom_point()

fielder_idx<-fiedler_tbl %>% 
  arrange(fielder) %>% 
  dplyr::select(bin) %>% 
  unlist
eigen_idx<-eig_bin_tbl %>% 
  arrange(eigen) %>% 
  dplyr::select(bin) %>% 
  unlist
image(as.matrix(cor_mat[eigen_idx,eigen_idx]),col=viridis(100),raster=T)
image(as.matrix(cor_mat[fielder_idx,fielder_idx]),col=viridis(100))

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

build_part_GRange_fn<-function(bin_tbl,tmp_res,chromo){
  GenomicRanges::reduce(GRanges(seqnames=chromo,
          ranges = IRanges(start=as.numeric(bin_tbl$bin),
                           end=as.numeric(bin_tbl$bin) + res_num[tmp_res] -1
          )))
}
#-----------------------------------------
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
HiC_spec_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
CAGE_GRange_file<-"~/Documents/multires_bhicect/Epigenome_enrichment_Bioconductor/data/CAGE_tbl/GM12878_CAGE_TSS_tbl.Rda"

CAGE_tbl<-get_tbl_in_fn(CAGE_GRange_file)

CAGE_GRange<-GRanges(seqnames=CAGE_tbl$chr,
                     ranges = IRanges(start=as.numeric(CAGE_tbl$start),
                                      end=CAGE_tbl$end
                     ))
mcols(CAGE_GRange)<-tibble(m=CAGE_tbl$m)
# For each chromosome
chr_set<-str_split_fixed(grep("^chr",list.files(HiC_spec_folder),value=T),"_",2)[,1]
chr_res_l<-vector('list',length(chr_set))
names(chr_res_l)<-chr_set
for (chromo in chr_set){
  message(chromo)
  chr_spec_res<-get_tbl_in_fn(paste0(HiC_spec_folder,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_lvl<-chr_bpt$Get('level')
  top_bpt_cl<-names(which(node_lvl==2))
  
  # Produce Compartment partition
  # Set resolution based on resolution of top BHiCect partition
  tmp_res<-unique(str_split_fixed(top_bpt_cl,"_",2)[,1])
  # Build z-score
  chr_dat<-compute_chr_res_zscore_fn(HiC_dat_folder,tmp_res,chromo,res_num)
  
  chr_mat<-full_f_mat(chr_dat,res_num[tmp_res],"zscore")
  
  # Build correlation matrix
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
  image(as.matrix(cor_mat),col=viridis(100),useRaster=T)
  # Produce eigen vector for compartment assignment
  eig_vec<-svd(cor_mat)$v[,1]
  names(eig_vec)<-colnames(ok_chr_mat)
  eig_bin_tbl<-tibble(bin=names(eig_vec),eigen=eig_vec) %>% 
    mutate(eigen.part=ifelse(eigen>0,1,2))
  
  # Produce table assigning bins to BHiCect top partitions
  bpt_bin_tbl<-tibble(bin=names(eig_vec)) %>% 
    mutate(bpt.part=ifelse(bin %in% chr_spec_res$cl_member[[top_bpt_cl[1]]],1,
                           ifelse(bin %in% chr_spec_res$cl_member[[top_bpt_cl[2]]],2,NA)))
  # Loop through each partition and compute CAGE-peak content
  
  
  chr_res_l[[chromo]]<-  do.call(bind_rows,lapply(unique(bpt_bin_tbl$bpt.part),function(i){
    tmp_bpt_part_tbl<-bpt_bin_tbl %>% filter(bpt.part==i)
    bpt_part_GRange<-build_part_GRange_fn(tmp_bpt_part_tbl,tmp_res, chromo)
    tmp_comp_part_tbl<-eig_bin_tbl %>% filter(eigen.part==i)
    comp_part_GRange<-build_part_GRange_fn(tmp_comp_part_tbl,tmp_res, chromo)
    bpt_part_dens<-sum(countOverlaps(bpt_part_GRange,CAGE_GRange))/sum(width(bpt_part_GRange))
    comp_part_dens<-sum(countOverlaps(comp_part_GRange,CAGE_GRange))/sum(width(comp_part_GRange))
    
    
    tibble(set="comp",CAGE.dens=comp_part_dens,part=i) %>% 
      bind_rows(.,tibble(set="bpt",CAGE.dens=bpt_part_dens,part=i))
  })) %>% 
    mutate(chr=chromo)
  
}

do.call(bind_rows,chr_res_l) %>% 
  ggplot(.,aes(x=set,CAGE.dens,color=as.factor(part)))+geom_point()+facet_grid(.~chr)


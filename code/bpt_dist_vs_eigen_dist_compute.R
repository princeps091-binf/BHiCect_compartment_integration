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
#-----------------------------------------
HiC_dat_folder<-"~/Documents/multires_bhicect/data/GM12878/"
HiC_spec_folder<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"



chromo<-"chr22"
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
# Make parallel
intersect_size<-unlist(lapply(
  lapply(
    pintersect(chr_bin_GRange[bin_cl_inter_tbl$queryHits],cl_GRange_l[bin_cl_inter_tbl$subjectHits]),
    width),
  sum))

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

bin_inter_tbl<-expand_grid(X1=names(bin_to_cl_vec),X2=names(bin_to_cl_vec)) %>% 
  mutate(cl.A=bin_to_cl_vec[X1],
         cl.B=bin_to_cl_vec[X2])

bin_inter_tbl<-bin_inter_tbl %>% 
  mutate(bpt.d=tmp_d[as.matrix(bin_inter_tbl[,3:4])])


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

eig_dist_tbl<-t(combn(names(eig_vec),2)) %>% 
  as_tibble %>% 
  mutate(binA.eig=eig_vec[V1],
         binB.eig=eig_vec[V2]) %>% 
  mutate(eig.dist=(abs(binA.eig-binB.eig)),
         same.comp=ifelse(sign(binA.eig)*sign(binB.eig)<0,"diff","same"))

bin_inter_tbl %>% 
  left_join(.,eig_dist_tbl,by=c("X1"="V1","X2"="V2")) %>% 
  filter(!(is.na(eig.dist))) %>%
  mutate(gdist=abs(as.numeric(X1)-as.numeric(X2))) %>% 
  ggplot(.,aes(bpt.d,eig.dist))+
  geom_smooth()+
  geom_point(alpha=0.01)+
  facet_grid(.~same.comp,scales="free")

bin_inter_tbl %>% 
  left_join(.,eig_dist_tbl,by=c("X1"="V1","X2"="V2")) %>% 
  filter(!(is.na(eig.dist))) %>%
  mutate(gdist=abs(as.numeric(X1)-as.numeric(X2))) %>% 
  ggplot(.,aes(same.comp,bpt.d))+
  geom_boxplot()

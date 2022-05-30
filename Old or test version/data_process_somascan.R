
library(tidyverse)

outputpath <- "/mnt2/zhijian/decode_ACE_ACE2/"

filename <- "10714_7_ACE_ACE.txt.gz"
protein <- sub(filename,pattern='.txt.gz',replace='')

a <- read.table(file = filename,header=T,stringsAsFactors=T) #about 30 million SNPs

A <- a %>% mutate(Z = Beta/SE) %>% select(,-Name)

#delete NA;delete duplicate
A <- A %>% drop_na(rsids) %>% distinct(rsids,.keep_all = TRUE) #about 25 million SNPs

top.idx <- order(abs(A$Z),decreasing=T)[1]
top.pos <- A[top.idx,'Pos']
cis.chr <- A[top.idx,'Chrom']
cis <- A %>% filter(Chrom == cis.chr) %>% filter(Pos < top.pos + 5E5 & Pos > top.pos - 5E5)

write.table(cis,file = paste0(outputpath,'cis_test.txt'),row.names=F,quote=F,sep='\t')

# lowmaf <- arrange(A,ImpMAF)[1:1E4,]
lowmaf.idx <- order(A$ImpMAF)[1:2E4]
lowmaf <- data.frame(rsids = A[lowmaf.idx,'rsids'],Z = A[lowmaf.idx,'Z'],stringsAsFactors=T)

write.table(lowmaf,file = paste0(outputpath,protein,'_lowmaf.txt'),quote=F,row.names=F,sep = '\t')

A1 <- A[nchar(A$effectAllele) <= 2 & nchar(A$otherAllele) <= 2,] %>% 
		filter(ImpMAF > 0.01)

snplist <- A1[,c('rsids','Chrom','Pos')]
A1 <- A1 %>% dplyr::select(,-Chrom,-Pos,-Pval,-minus_log10_pval,-Z) # 0.986Gb in R

t1 <- Sys.time()
saveRDS(A1,file = paste0(outputpath,protein,'.rds')) # 138M 
t2 <- Sys.time()

print(t2 - t1) # 27.10122 secs

A.TGCA <- A[nchar(A$effectAllele) == 1 & nchar(A$otherAllele) == 1,] %>% 
		filter(ImpMAF > 0.05)

A.TGCA <- A.TGCA %>% dplyr::select(,rsids,N,Z)

saveRDS(A.TGCA,file='test_tgca.rds')


# cis.file.list <- read.delim('/opt/share2Big/somascan_data/cis_list.txt',h=T,stringsAsFactors=F)
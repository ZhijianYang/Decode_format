##

# This script should be ran on server Bigblack

setwd('/Volumes/share2Big/somascan_data/')

args <- commandArgs(trailingOnly = TRUE)
start_num <- args[1]
end_num <- args[2]

library(tidyverse)

sink(file = paste0('process_somascan_',start_num,'_',end_num,'.log'),split = TRUE)

print('MAF threshold 0.01 for rds, 0.05 for TGCA',quote = F)

cis.file.list <- read.delim(file = 'cis_list.txt',header = T,stringsAsFactors = F)
inputpath <- "/Volumes/Seagate\ Basic/decode.is/"
outputpath <- "/Volumes/share2Big/somascan_data/"

SNP_list_all <- SNP_list_intersect <- c()

for (i in start_num:end_num) {

	filename <- cis.file.list[i,'filename']
	symbol <- paste(cis.file.list[i,1],cis.file.list[i,4],sep = '_')
	
	print(Sys.time(),quote = F)
	cat(paste0(i,'  Processing ',symbol,'\n'))
	cat(paste0('reading file ',filename,' ...','\n'))
	t1 <- Sys.time()
	a <- tryCatch(read.table(file = paste0(inputpath,filename),header = T,stringsAsFactors = F),error = function(cond) {message(cond);return(cond)}) # about 30 million SNPs
	t2 <- Sys.time()
	print(t2 - t1)

	if (inherits(a,'error')) {
		cat('\n')
		cat('Failed to open file ',filename,'\n')
		print(a)
		next
	}

	A <- a %>% mutate(Z = Beta/SE) %>% select(,-Name)
	
	#delete NA;delete duplicate
	A <- A %>% drop_na(rsids) %>% distinct(rsids,.keep_all = TRUE) # about 25 million SNPs
	
	top.idx <- order(abs(A$Z),decreasing = T)[1]
	top.pos <- A[top.idx,'Pos']
	cis.chr <- A[top.idx,'Chrom']
	# 1Mb window centred top cis-pQTL
	cis <- A %>% filter(Chrom == cis.chr) %>% filter(Pos < top.pos + 5E5 & Pos > top.pos - 5E5)
	write.table(cis,file = paste0(outputpath,'cis_1Mb/',symbol,'_cis.txt'),row.names = F,quote = F,sep = '\t')

	# lowmaf <- arrange(A,ImpMAF)[1:1E4,]
	lowmaf.idx <- order(A$ImpMAF)[1:2E4]
	lowmaf <- data.frame(rsids = A[lowmaf.idx,'rsids'],Z = A[lowmaf.idx,'Z'],stringsAsFactors = F)
	
	write.table(lowmaf,file = paste0(outputpath,'lowmaf/',symbol,'_lowmaf.txt'),quote = F,row.names = F,sep = '\t')
	
	A1 <- A[nchar(A$effectAllele) <= 2 & nchar(A$otherAllele) <= 2,] %>% 
			filter(ImpMAF > 0.01)
	
	snplist <- A1[,c('rsids','Chrom','Pos','effectAllele','otherAllele','N','ImpMAF')]
	if (is.null(SNP_list_all)) {
		SNP_list_all <- snplist
		SNP_list_intersect <- snplist
	} else {
		SNP_list_all <- full_join(SNP_list_all,snplist,by = c('rsids','Chrom','Pos','effectAllele','otherAllele','N','ImpMAF'))
		SNP_list_intersect <- inner_join(SNP_list_intersect,snplist,by = c('rsids','Chrom','Pos','effectAllele','otherAllele','N','ImpMAF'))
	}
	
	A1 <- A1 %>% dplyr::select(,-Chrom,-Pos,-Pval,-minus_log10_pval,-Z) # 0.986Gb in R
	
	cat(paste0(nrow(A1),' SNPs remained in rds file','\n'))
	saveRDS(A1,file = paste0(outputpath,'RDS/',symbol,'.rds')) # 138M 

	A2 <- A1 %>% dplyr::select(,rsids,Beta,SE)
	saveRDS(A2,file = paste0(outputpath,'Multi_RDS/',symbol,'.rds'))

	A.TGCA <- A[nchar(A$effectAllele) == 1 & nchar(A$otherAllele) == 1,] %>% 
		filter(ImpMAF > 0.05)

	A.TGCA <- A.TGCA %>% dplyr::select(,rsids,N,Z)

	cat(paste0(nrow(A.TGCA),' SNPs remained for TGCA','\n'))
	saveRDS(A.TGCA,file = paste0(outputpath,'TGCA_format/',symbol,'_tgca.rds'))

	# rm(a,A)
	gc()
}

cat('The whole SNP list saved in SNP_list_all.rds and SNPs shared by all protein in SNP_list_intersect.rds \n')
saveRDS(SNP_list_all,file = paste0(outputpath,'SNP_list_all_',start_num,'_',end_num,'.rds'))
saveRDS(SNP_list_intersect,file = paste0(outputpath,'SNP_list_intersect_',start_num,'_',end_num,'.rds'))

sink()



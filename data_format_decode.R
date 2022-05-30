## data process decode
## This script should be ran on LB

# t1 <- Sys.time()
# t2 <- Sys.time()
# print(t2 - t1)

## read input arguments
args <- commandArgs(trailingOnly = TRUE)
start_num <- args[1]
end_num <- args[2]

setwd('/opt/storage/deCODE_formatted/zhijian/somascan_data/')
sink(file = paste0('process_decode_',start_num,'_',end_num,'.log'),split = TRUE)

library(tidyverse)

cat('################################################################','\n')
print('MAF threshold 0.01 for rds file',quote = F)
cat('################################################################','\n')

inputpath <- "/NiuHe/decode.is/"
outputpath <- "/opt/storage/deCODE_formatted/zhijian/somascan_data/"

cis.info <- read.delim('/opt/share2Big/somascan_info/cis_info.txt',header = T, stringsAsFactors = F)
cat('Reading annotated.txt.gz and excluded.txt.gz \n')
cat('\n')
annot <- read.table('/NiuHe/assocvariants.annotated.txt.gz',header = T,stringsAsFactors = F)
exclude <- read.table('/NiuHe/assocvariants.excluded.txt.gz',header = T,stringsAsFactors = F)

annot <- annot %>% select(-Chrom,-Pos) %>% filter(rsids != '.')
# cis.info[cis.info$Chr == 'chrX',] <- 'chr23'
# cis.info <- cis.info %>% mutate(Chr = as.numeric(sub(Chr,pattern = 'chr',replacement = '')))

SNP_list_all <- SNP_list_intersect <- c()

for (i in start_num:end_num) {

	filename <- cis.info[i,'filename']
	symbol <- paste(cis.info[i,'SeqId'],cis.info[i,'Gene'],sep = '_')
	
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

	#delete NA;delete duplicate
	t1 <- Sys.time()

	A <- a %>% drop_na(rsids)  # about 26 million SNPs
	A <- A %>% filter(!Name %in% exclude$Name) %>% inner_join(.,annot,by = c('rsids','effectAllele'))# about 23 million SNPs
	A <- A %>% mutate(Z = Beta/SE) %>% select(,c(-minus_log10_pval,-otherAllele.y,-Name.y))

	t2 <- Sys.time()
	print(t2 - t1)

	colnames(A) <- sub(colnames(A),pattern='.x',replacement='')

	# A[A$Chrom == 'chrX','Chrom'] <- 'chr23'
	# A <- A %>% mutate(Chrom = as.numeric(sub(Chrom,pattern = 'chr',replacement = '')))

	t1 <- Sys.time()
	cat('Saving cis region ... \t')
	cis.chr <- cis.info$Chr[i]
	A.chr <- A[A$Chrom == cis.chr,]
	top.idx <- order(abs(A.chr$Z),decreasing = T)[1]
	top.pos <- A.chr[top.idx,'Pos']
	# 1Mb window centred top cis-pQTL
	cis <- A.chr %>% filter(Pos < top.pos + 5E5 & Pos > top.pos - 5E5)
	write.table(cis,file = paste0(outputpath,'cis_1Mb/',symbol,'_cis.txt'),row.names = F,quote = F,sep = '\t')

	A <- A %>% distinct(rsids,.keep_all = TRUE)

	cat('Saving 2E4 lowmaf SNPs ... \n')
	lowmaf.idx <- order(A$ImpMAF)[1:2E4]
	lowmaf <- data.frame(rsids = A[lowmaf.idx,'rsids'],Z = A[lowmaf.idx,'Z'],stringsAsFactors = F)
	
	write.table(lowmaf,file = paste0(outputpath,'lowmaf/',symbol,'_lowmaf.txt'),quote = F,row.names = F,sep = '\t')

	t2 <- Sys.time()
	print(t2 - t1)
	
	A1 <- A %>% filter(ImpMAF > 0.01)
	
	snplist <- A1[,c('rsids','Chrom','Pos','effectAllele','otherAllele','N','effectAlleleFreq')]
	if (is.null(SNP_list_all)) {
		SNP_list_all <- snplist
		SNP_list_intersect <- snplist
	} else {
		SNP_list_all <- full_join(SNP_list_all,snplist,by = c('rsids','Chrom','Pos','effectAllele','otherAllele','N','effectAlleleFreq'))
		SNP_list_intersect <- inner_join(SNP_list_intersect,snplist,by = c('rsids','Chrom','Pos','effectAllele','otherAllele','N','effectAlleleFreq'))
	}
	
	A1 <- A1 %>% dplyr::select(-Name,-Chrom,-Pos,-Pval,-Z,-ImpMAF) # 0.986Gb in R
	
	cat(paste0(nrow(A1),' SNPs remained in rds file','\n'))
	saveRDS(A1,file = paste0(outputpath,'RDS/',symbol,'.rds')) # ~ 174M

	A1 <- A1 %>% dplyr::select(rsids,Beta,SE)
	saveRDS(A1,file = paste0(outputpath,'B_SE/',symbol,'.rds'))

	cat('\n')

	gc()
}

cat('The whole SNP list saved in SNP_list_all.rds and SNPs shared by all protein in SNP_list_intersect.rds \n')
saveRDS(SNP_list_all,file = paste0(outputpath,'SNP_list_all_',start_num,'_',end_num,'.rds'))
saveRDS(SNP_list_intersect,file = paste0(outputpath,'SNP_list_intersect_',start_num,'_',end_num,'.rds'))

sink()



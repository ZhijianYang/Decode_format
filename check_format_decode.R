##

library(data.table)

# check size of the origin files
filelist <- list.files('/NiuHe/decode.is/',pattern = 'txt.gz')
length(filelist)

table(fls > 940000000) # TRUE 4907

################################################################

setwd('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/')
filelist = list.files('SNP_list/',pattern = '.rds')

for (i in 1:length(filelist)) {
	b <- readRDS(paste0('SNP_list/',filelist[i]))

	cat(sub(filelist[i],pattern = '.rds',replacement = ''),'\n')

	if (nrow(b) == 10293879) {
		cat('rownum: ',nrow(b),', checked. \n')
	} else {
		cat('rownum: ',nrow(b),', check failed. \n')
	}

	if (i == 1) {
		SNP_list <- b$rsids
		a <- b
	} else {
		cat('All rsids matched',all(SNP_list == b$rsids),'\n')
		cat('All effectAlleles matched.',all(a$effectAllele == b$effectAllele),'\n')
		cat('All otherAlleles matched.',all(a$otherAllele == b$otherAllele),'\n')
	}
	
	cat('\n')

}

saveRDS(b,'SNP_list.rds')

################################################################

library(data.table)

setwd('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/')

write.line <- function(x, file, append = TRUE) {
  write.table(x, file, row.names = FALSE, col.names = FALSE, quote = FALSE, append = append)
}

write.line(x = 'protein that rownum does not equal to 10293879',file = 'logs/rownum_fails.txt',append = FALSE)
write.line(x = 'protein that rsids, effectAllele, otherAllele are not matched',file = 'logs/rsids_fails.txt',append = FALSE)

sink('logs/check_format_decode.log',split = TRUE)

## The 10293879 SNPs list
b <- readRDS('SNP_list.rds')

filelist <- list.files('sumstat/',pattern = 'txt.gz')

for (i in 1:length(filelist)) {

	protein <- sub(filelist[i],pattern='.txt.gz',replacement = '')
	cat(i,'\t',protein,'\n')

	a <- fread(paste0('sumstat/',filelist[i]))
	if (nrow(a) == 10293879) {
		cat('rownum: ',nrow(a),', checked. \n')
	} else {
		cat('rownum: ',nrow(a),', check failed. \n')
		write.line(x = protein,file = 'rownum_fails.txt')
	}
	
	if (all(a$rsids == b$rsids) & all(a$effectAllele == b$effectAllele) & all(a$otherAllele == b$otherAllele)) {
		cat('rsids,effectAllele,otherAllele are matched. \n')
	} else {
		cat('rsids,effectAllele,otherAllele are not matched. \n')
		write.line(x = protein,file = 'rsids_fails.txt')
	}
	cat('\n')
}

sink()

################################################################

# shell 

# zhijian@LargeBlack:/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/logs$ grep 'rsids,effectAllele,otherAllele are matched' check_format_decode.log| wc -l
# 1881

# zhijian@LargeBlack:/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/logs$ grep 'checked' check_format_decode.log| wc -l
# 1881

################################################################


library(data.table)

setwd('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/')

filelist <- list.files('sumstat/',pattern = 'txt.gz')

if (!dir.exists("/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/B_SE_N/")){dir.create("/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/B_SE_N/")}

for (i in 1:length(filelist)) {

	print(Sys.time())
	protein <- sub(filelist[i],pattern='.txt.gz',replacement = '')
	cat(i,'\t',protein,'\n')

	if (file.exists(paste0('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/B_SE_N/',protein,'.txt.gz'))) {
		cat('Result file exist, next \n')
		next
	}

	a <- fread(paste0('sumstat/',filelist[i]),data.table = FALSE)

	b <- a[,c('Beta','SE','N')]

	write.table(b,paste0('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/B_SE_N/',protein,'.txt'),row.names = F,quote = F)
	system(paste0('gzip ','/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/B_SE_N/',protein,'.txt'))

}

################################################################


library(data.table)

setwd('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/')

filelist <- list.files('B_SE_N/',pattern = 'txt.gz')

k <- 1

sink(paste0('logs/b_se_check_n_',k,'.log'),split = TRUE)

if (!dir.exists("/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/b_se/")){dir.create("/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/b_se/")}

a <- fread(paste0('B_SE_N/',filelist[1]),data.table = FALSE)
n <- a$N

num <- length(filelist) - 1
start_num <- num*(k-1)/5 + 1
end_num <- num*k/5

if (end_num == length(filelist) - 1) {
	end_num <- end_num + 1
}

cat('processing ',start_num,' to ',end_num,' ... \n')

for (i in start_num:end_num) {

	cat('\n')
	print(Sys.time())
	protein <- sub(filelist[i],pattern='.txt.gz',replacement = '')
	cat(i,'\t',protein,'\n')

	if (file.exists(paste0('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/b_se/',protein,'.txt.gz'))) {
		cat('Result file exist, next \n')
		next
	}

	a <- fread(paste0('B_SE_N/',filelist[i]),data.table = FALSE)

	if (all(a$N == n)) {

		cat('n check:',all(a$N == n),'\n')
		b <- a[,c('Beta','SE')]

		write.table(b,paste0('/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/b_se/',protein,'.txt'),row.names = F,quote = F)
		system(paste0('gzip ','/opt/storage/deCODE_formatted/zhijian/decode_somascan_data/b_se/',protein,'.txt'))
		cat('write out done. \n')

	} else {

		cat('n check:',all(a$N == n),', next \n')
		next


	}
	

}

sink()





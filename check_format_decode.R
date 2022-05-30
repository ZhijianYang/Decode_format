##


# check size of the origin files
filelist <- list.files('/NiuHe/decode.is/',pattern = 'txt.gz')
length(filelist)

table(fls > 940000000) # TRUE 4907

################################################################

setwd('/opt/storage/deCODE_formatted/zhijian/somascan_data')

write.line <- function(x, file, append = TRUE) {
  write.table(x, file, row.names = FALSE, col.names = FALSE, quote = FALSE, append = append)
}

write.line(x = 'protein that rownum does not equal to 10293879',file = 'rownum_fails.txt',append = FALSE)
write.line(x = 'protein that rsids,effectAllele,otherAllele are not matched',file = 'rsids_fails.txt',append = FALSE)

## small typo here
sink('check_foemat_decode.log',split = TRUE)

## The 10293879 SNPs list
b <- readRDS('SNP_list_all_943_1413.rds')

filelist <- list.files('RDS/',pattern = '.rds')

for (i in 1:length(filelist)) {
	protein <- sub(filelist[i],pattern='.rds',replacement = '')
	cat(i,'\t',protein,'\n')

	a <- readRDS(paste0('RDS/',filelist[i]))
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



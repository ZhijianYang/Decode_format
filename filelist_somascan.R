##

# This script shoud be ran on bigblack

setwd('/Volumes/Seagate\ Basic')

filelist <- list.files('decode.is/',pattern = '.gz') # 4907 files

sr <- function(a) {
	a <- gsub(a,pattern = '.txt.gz',replace	= '')
	b <- strsplit(a,split = '_')
	c <- data.frame(seqId = sapply(b,function(x) paste(x[1],x[2],sep='_')),
			gene_protein = sapply(b,function(x) paste(x[-c(1,2)],collapse='_')),
			stringsAsFactors = F)
	return(c)
}

files <- sr(filelist)

files$filename <- filelist

## the file produced from Supplementary table of the paper with EXCEL
info <- read.csv('/Volumes/share2Big/somascan_info/ss_pQTL_asso_light.csv',h = T,stringsAsFactors = F,skip = 2)
info <- info[,c('pQTL_ID...global.','region_ID...prot..','gene...prot..','shortname..prot..','UniProt','SeqId','variant','cis..trans')]
colnames(info) <- sub(colnames(info),pattern='..prot..',replace='')

info <- info[info$cis..trans == 'cis',] #  cis 7572 trans 20619

info1 <- info[!duplicated(info$SeqId),]# 1881 seqs
cis.info <- data.frame(seqId = info1$SeqId, UniProt = info1$UniProt, Gene = info1$gene., Protein = info1$shortname, stringsAsFactors = F)

require(dplyr)
res <- inner_join(files,cis.info,by='seqId') # 1881 seqs

res <- select(res,-gene_protein)

write.table(res,'/Volumes/share2Big/somascan_data/cis_list.txt',sep = "\t",row.names=F,quote=F)

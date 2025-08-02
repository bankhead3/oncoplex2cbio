# 20250801 arb
# get chrom per gene

options(stringsAsFactors=F)

library(dplyr)

inFile1 = '../oncoplex-genes.txt'
inFile2 = '../input/MANE.GRCh38.v1.4.summary.txt'

data1 = read.delim(inFile1)
data2 = read.delim(inFile2)

data2a = select(data2,gene=symbol,chrom1=GRCh38_chr)
data2a$chrom = sub('NC_[0]+','',data2a$chrom1)
data2a$chrom = sub('[.].*$','',data2a$chrom)
data2b = unique(select(data2a,gene,chrom))
data2b$chrom = paste0('chr',data2b$chrom)
data2b$chrom[data2b$chrom == 'chr23'] = 'chrX'

data1a = merge(data1,data2b,by='gene',all.x=T)

data1a$chrom[data1a$gene == 'C19MC'] = 'chr19'

data1a = data1a[order(data1a$gene),]
write.table(data1a,inFile1,quote=F,row.names=F,sep="\t")

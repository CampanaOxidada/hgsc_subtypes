library(data.table)
library(Biobase)


aaces.data <- '1.DataInclusion/Data/Aaces'

# read aaces expression data
pcl.file <- file.path(aaces.data, 'aaces.expr.pcl')

expr.rma.gene.df <- read.table(file = pcl.file, header = FALSE, sep = ' ', skip = 1, row.names = 1)

sample.name <- read.table(file = pcl.file, nrows  = 1L, sep = ' ')

names(expr.rma.gene.df) <- sample.name

# get metadata
dat <- fread(file.path(aaces.data, 'AACESmetaData.csv'))

dat[,V1 := NULL]

mdat <- dat[dat[,use]]

# reorder expr.rma.gene columns to match mdat
expr.rma.gene.df <- expr.rma.gene.df[,as.character(mdat[,sampleid])]

expr.rma.gene <- as.matrix(expr.rma.gene.df)


pdata <- with(mdat, data.frame(row.names = sampleid, sample_type = 'tumor', histological_type = hist, grade = 3))


var.mdat <- data.frame(
		    col.name = c('sample_type','histological_type','grade'),
		    var.class = c('character','character','integer'),
		    allowedvalues = c(
				    'tumor|metastatic|borderline|benign|adjacentnormal|healthy|cellline',
				    'ser|endo|clearcell|mucinous|undifferentiated|other|mix',
				    '0|1|2|3|4'),
		    description = c('healthy should be only from individuals without cancer, adjacentnormal from individuals with cancer, metastatic for non-primary tumors, borderline includes both borderline and LMP tumors, benign for benign tumors.',
				    'ser=serous;endo=endometrioid;clearcell=mixture of ser+endo.  Other includes sarcomatoid, endometroid, papillary serous, unspecified adenocarcinoma, dysgerminoma',
				    'G (0-4): If multiple given, ie 12, 23, use highest given')
		    )

aaces.eset <- ExpressionSet(expr.rma.gene, phenoData = AnnotatedDataFrame(pdata, varMetadata = var.mdat))

save(aaces.eset, file = file.path(aaces.data, 'aaces.eset.RData'))

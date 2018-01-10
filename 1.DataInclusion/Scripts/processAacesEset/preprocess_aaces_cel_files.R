
###########################################################################
#LOAD REQUIRED LIBRARIES
###########################################################################

#Use the following to install required packages
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("oligo", "pd.hta.2.0", "hta20stprobeset.db", "hta20sttranscriptcluster.db", "affxparser", "plyr", "sva"))
#biocLite("hta20transcriptcluster.db")

library(oligo)                        # whence RMA for normalization
library(pd.hta.2.0)                   # HTA 2.0 probe annotations
library(hta20transcriptcluster.db)    #Provides the mapping from probe/transcript IDs to gene SYMBOLS
library(data.table)                   #Convenience package used to handle data frames/matrices
library(sva)                          # whence ComBat
library(WGCNA)                        # whence collapseRows


###########################################################################
#LOAD AND PROCESS CEL FILES
###########################################################################

aaces.data.dir <- '1.DataInclusion/Data/Aaces'

cel.file.dir <- file.path(aaces.data.dir,'expression')

cel.pathname <- c(list.celfiles(cel.file.dir, full.names = TRUE))

#Read the CEL files
eset.cel <- read.celfiles(cel.pathname)

#load metadata from CSV

meta.data.file <- file.path(aaces.data.dir, 'AACESmetaData.csv')

dat <- fread(meta.data.file, key = 'sampleid', colClasses = c(sampleid = 'character'))

dat[,V1 := NULL]

# limit to qualifying (use = TRUE) samples 
mdat <- dat[dat[,use]]

#Limit to only samples that are useable 
filename.use <- intersect(sampleNames(eset.cel), mdat[,exp_file_name])

#oligo::rma normalize on core transcript probes
eset.rma.trans <- oligo::rma(eset.cel[,filename.use], target = 'core')

# sort by file name 
setkey(mdat, exp_file_name) 
sampleNames(eset.rma.trans) <- mdat[sampleNames(eset.rma.trans),sampleid]

#Free up some RAM now that we have normalized expression values
rm(eset.cel)
gc()

# expression matrix
expr.rma.trans <- exprs(eset.rma.trans)

###########################################################################
#REMOVE THE SITE BATCH AFFECT USING COMBAT
###########################################################################
# define batch solely on source of microarray data
# Perform the batch correction therefrom
fac <- factor(mdat[,exp_facility])

modl <- model.matrix(~1, data = as.data.frame(t(expr.rma.trans)))

expr.rma.bat.gene <- ComBat(dat = expr.rma.trans, batch = fac, mod = modl)


###########################################################################
#SUMMARIZE EXPRESSION BY GENE SYMBOL
###########################################################################

# rename probe features after gene symbols
#sym.map <- hta20sttranscriptclusterSYMBOL  #This has been depcreciated and is no longer available for the latest version of R
sym.map <- hta20transcriptclusterSYMBOL
mapped_trans <- unlist(as.list(sym.map[mappedkeys(sym.map)]), use.names = TRUE)
trans.data <- data.table(gene.symbol = mapped_trans, probe.id = names(mapped_trans))

setorder(trans.data, gene.symbol) #Sort alphabetically by gene SYMBOL
expr.rma.genes <- expr.rma.trans[trans.data[,probe.id],]

# collapse rows to single symboles using WGCA's default max average
collapse.list <- with(trans.data,
                      
                      collapseRows(expr.rma.genes, rowGroup = gene.symbol, rowID = probe.id)
)
expr.rma.gene <- collapse.list[['datETcollapsed']]

output.file.pcl <- file.path(aaces.data.dir, "aaces.expr.pcl")
output.file.rda <- file.path(aaces.data.dir, "aaces.eset.RData")

write.table(expr.rma.gene, file = output.file.pcl, row.names = TRUE, col.names = TRUE)

save(expr.rma.gene, file = output.file.rda)

# CEL files from Affymetrix arrays

library(limma)
library(Biobase)
library(Biostrings)
library(genefilter)
library(oligo)
library(ggplot2)

celPath <- "D:/DataScience/Profile/MultipleMyeloma_GPL570/celFiles"

# import CEL files containing raw probe_level into R AffyBatch object
celList <- list.files(celPath, full.names = T)
celData <- read.celfiles(celList)

# retrieve intensities

Expr <- exprs(celData)
# alternative: intensity(celData)
# alternative: pm(celData), reorders the rows by probe numbers
Expr[1:5,]


######### retrieving annotation

# sample annotation
ph <- celData@phenoData
ph@data
pData(celData)

# probe annotation
feat <- celData@featureData
feat@data # featureData is not defined

# experiment annotation
expAnn <- celData@experimentData
expAnn # not defined

# IDs of the probes (numbers of probe sets)
length(featureNames(celData))

# probe names (number of probe names)
length(probeNames(celData))
probeNames(celData)



#############  Quality control

# give the samples informative names
 ph@data[,1] <- c("KMS-11/Cfz.1", "KMS-11/Cfz.2","KMS-11/Cfz.3","KMS-34/Cfz.1","KMS-34/Cfz.2","KMS-34/Cfz.3",
                "KMS-11_1","KMS-11_2","KMS-11_3","KMS-34_1","KMS-34_2","KMS-34_3")
ph@data


### microarray pictures
op = par(mfrow = c(4,3))
op
for (i in 1:12){
image(celData[,i], main = ph@data$index[i])
}


### Chip pseudo image
Pset <- fitProbeLevelModel(celData)
par(mfrow = c(4,3))
for (i in 1:12){
  image(Pset, which = i, type = "residuals", main = ph@data$index[i]) # based on residuals
}# Positive residuals are plotted in red, negative residuals in blue.

# positive: the intensity of the probe is larger than the ideal value according to the model
# negative: the intensity of the probe is smaller than the ideal value according to the model


par(mfrow = c(4,3))
for (i in 1:12){
  image(Pset, which = i, main = ph@data$index[i]) # based on weights
}# Small weights (outliers) are indicated in green on the figure.


### histogram
# prepare dataframe with long intensities in one column and sample names in another column
pmexp <- pm(celData)

sampleName <- vector()
logs <- vector()
for(i in 1:12){
  sampleName <- c(sampleName, rep(ph@data[i,1], dim(pmexp)[1]))
  logs <- c(logs, log2(pmexp)[,i])
}

# single data frame
logData <- data.frame(logInt = logs, sampleName = sampleName)

logHist <- ggplot(logData, aes(logInt, colour = sampleName))
logHist <- logHist + geom_density()
logHist # differences in shape or center means the need of normalization

### Box plot
logBox <- ggplot(logData, aes(sampleName, logInt))
logBox = logBox + geom_boxplot() + theme(axis.text.x = element_text(angle = 45))
logBox # different median and scale of the boxes: needs normalization


### MA plot
# The MA plot shows to what extent the variability in expression depends on the expression level (more variation on high expression values?)

# M is the difference between the intensity of a probe on the array and the median intensity of that probe over all arrays
# M = logPMInt_array - logPMInt_medianarray
# A is the average of the intensity of a probe on that array and the median intesity of that probe over all arrays
# A = (logPMInt_array + logPMInt_medianarray)/2

for(i in 1:12){
  MAplot(celData, which = i)
}# Ideally, the cloud of data points should be centered around M=0 (blue line). Needs normalization.


######## Normalization

## using RMA


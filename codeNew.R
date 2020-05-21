#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation
library(hgu133plus2.db)


#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

#Plotting and color options packages
#library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)



############

celPath <- "D:/DataScience/Profile/MultipleMyeloma_GPL570/celFiles"

# import CEL files containing raw probe_level into R AffyBatch object
celList <- list.files(celPath, full.names = T)
celData <- read.celfiles(celList)

# phenodata
ph <- celData@phenoData
ph@data[,1] <- c("KMS11.Cfz.1", "KMS11.Cfz.2","KMS11.Cfz.3","KMS34.Cfz.1","KMS34.Cfz.2","KMS34.Cfz.3",
                 "KMS11.1","KMS11.2","KMS11.3","KMS34.1","KMS34.2","KMS34.3")
colnames(ph@data)[1] <- "Samples"

ph@data[,2] <- c("KMS.11.Cfz","KMS.11.Cfz","KMS.11.Cfz","KMS.34.Cfz","KMS.34.Cfz","KMS.34.Cfz", "KMS.11","KMS.11","KMS.11","KMS.34","KMS.34","KMS.34")
colnames(ph@data)[2] <- "Groups"

ph@data[,3] <- c("KMS11","KMS11","KMS11","KMS34","KMS34","KMS34", "KMS11","KMS11","KMS11","KMS34","KMS34","KMS34")
colnames(ph@data)[3] <- "CellLine"

ph@data[,4] <- c("Cfz", "Cfz", "Cfz", "Cfz", "Cfz", "Cfz", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl")
colnames(ph@data)[4] <- "Treatment"
###
colnames(celData)<- ph@data$Samples
pData(celData)[,1] <- c("KMS.11.Cfz","KMS.11.Cfz","KMS.11.Cfz","KMS.34.Cfz","KMS.34.Cfz","KMS.34.Cfz", "KMS.11","KMS.11","KMS.11","KMS.34","KMS.34","KMS.34")
pData(celData)[,2] <- c("KMS11","KMS11","KMS11","KMS34","KMS34","KMS34", "KMS11","KMS11","KMS11","KMS34","KMS34","KMS34")
pData(celData)[,3] <- c("Cfz", "Cfz", "Cfz", "Cfz", "Cfz", "Cfz", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl")
names(celData@phenoData@data) <- c("Groups", "Cells", "Treatments")

####### Quality control of the raw data

exprs(celData)[1:5,1:3]

### take log2 of expression and perform a PCA 
exp_raw <- log2(exprs(celData))
PCA_raw <- prcomp(t(exp_raw), sacle. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2]/percentVar[1])

ggData <- data.frame(PC1 = PCA_raw$x[,1],
                     PC2 = PCA_raw$x[,2],
                     Samples = ph@data$Samples,
                     Cells = ph@data$CellLine,
                     Treatment = ph@data$Treatment)

P1 <- ggplot(ggData, aes(PC1, PC2)) + geom_point(aes(shape = Cells, colour = Treatment))+
  ggtitle("PCA of log_transformed raw expression data") + xlab(paste0("PC1, VarExp:", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp:", percentVar[2], "%")) + theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) + scale_shape_manual(values = c(4,15)) + scale_color_manual(values = c("darkorange2", "dodgerblue4"))
P1 # cell lines are differentiated by second component!


### look at log intensity on boxplot from oligo package
set.seed(542)
oligo::boxplot(celData, target = "core", main = "log2-intensities of raw data")  # different median and scale of the boxes: needs normalization

# advanced plots  
arrayQualityMetrics(expressionset = celData, outdir = tempdir(),
                    force = TRUE, do.logtransform = TRUE, intgroup = c("Cells", "Treatments"))


### RLE (Relative Log Expression) data quality analysis
rleData <- rma(celData, normalize = FALSE) # without normalization to get only log2 values

# The RLE is performed by calculating the median log2 intensity of every transcript across all arrays.
row_medians <- Biobase::rowMedians(as.matrix(exprs(rleData)))

#We then substract this transcript median intensity from every transcript intensity via the sweep function.
RLE_Data <- sweep(exprs(rleData), 1, row_medians)
RLE_Data <- as.data.frame(RLE_Data)

Rle_Data_gathered <- gather(RLE_Data, Treatment, log2_expression_deviation)

ggplot(Rle_Data_gathered, aes(Treatment, log2_expression_deviation)) +
  geom_boxplot(outlier.shape = NA) + ylim(c(-2,2)) + 
  theme(axis.text.x = element_text(colour = "darkgreen", angle = 60, size = 6.5, hjust = 1, face = "bold"))


###### RMA calibration
eset_norm <- rma(celData) 


### quality assessment boxplot
set.seed(754)
oligo::boxplot(eset_norm, target = "core", main = "log2-intensities of normalized data")

### quality assessment PCA plot
norm_exp <- exprs(eset_norm)
PCA <- prcomp(t(norm_exp), scale. = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2]/percentVar[1])
ggData <- data.frame(PC1 = PCA$x[,1],
                     PC2 = PCA$x[,2],
                     Cells = pData(eset_norm)$Cells,
                     Treatment = pData(eset_norm)$Treatment)

ggplot(ggData, aes(PC1, PC2)) +
  geom_point(aes(shape = Cells, colour = Treatment)) + 
  ggtitle("PCA plot of normalized data") + 
  xlab(paste0("PC1, VarExp:", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp:", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) + 
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
# now, the first pricipal component separates between cell lines, meaning that the differntial expression between cell lines is dominant.

### Heatmap clustering analysis
treat_name <- ifelse(str_detect(pData(eset_norm)$Treatment, "Cfz"), "Cfz", "Ctrl")
cells_name <- ifelse(str_detect(pData(eset_norm)$Cells, "KMS11"), "KMS11", "KMS34")

ann_heatmap <- data.frame(Cells = cells_name, Treatment = treat_name)
row.names(ann_heatmap) <- row.names(pData(eset_norm))

# calculate distance
dists <- as.matrix(dist(t(norm_exp), method = "manhattan"))

rownames(dists) <- row.names(pData(eset_norm)) 
hmcol <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA

ann_colors <- list(
  Treatment = c(Ctrl = "chartreuse4", Cfz = "burlywood3"),
  Cells = c(KMS11 = "blue4", KMS34 = "cadetblue2")
)

pheatmap(dists, color = (hmcol),
         annotation_row = ann_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), max(dists, na.rm = TRUE)),
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap of normalized samples")

# confirming PCA analysis, that samples are clustering by cell lines.

#### filtering out the low expressed genes

eset_medians <- rowMedians(norm_exp)
hist(eset_medians, 100, col = "cornsilk1", freq = FALSE,
     main = "Histogram of the median intensities",
     border = "antiquewhite4", xlab = "Median intensities")
# set cutoff threshold left to the histogram peak
cut_thereshold <- 4
abline(v=cut_thereshold, col = "coral4", lwd = 2)

# exclude the transcripts left to the cutoff border
no_of_samples <- table(paste0(pData(eset_norm)$Cells, "_", pData(eset_norm)$Treatment))
no_of_samples # equaliy 3, could be different in other setups then we cinsider the min number

sample_cutoff <- min(no_of_samples) # in this case 3
idx_cut_threshold <- apply(norm_exp, 1, function(x){
  sum(x>cut_thereshold) >= sample_cutoff}) # number of arrays where the median intesity passes the threshold is greater than the sample cutoff

table(idx_cut_threshold) # 3437 arrays can be removed

eset_norm_filtered <- subset(eset_norm, idx_cut_threshold)

##### Annotation of the transcript clusters
anno_eset <- AnnotationDbi::select(hgu133plus2.db, keys = (featureNames(eset_norm_filtered)),
                                   columns = c("SYMBOL", "GENENAME"),
                                   keytype = "PROBEID")
anno_eset <- subset(anno_eset, !is.na(SYMBOL)) # filterout probes that do not map to a gene

### removing multiple mappings

# group by probeID
anno_grouped <- group_by(anno_eset, PROBEID)
anno_summarized <- summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1) # filter for PROBEIDs with multiple matches
head(anno_filtered)
nrow(anno_filtered) # 2157 PROBEID with multiple matches, difficult to decide which mapping is correct. Therefoe we exclude these clusters


ids_to_exclude <- (featureNames(eset_norm_filtered) %in% anno_filtered$PROBEID)
table(ids_to_exclude)

# subset to final eset
eset_final <- subset(eset_norm_filtered, !ids_to_exclude)
validObject(eset_final)# TRUE

# we should also remove probe IDs from anno_eset 
head(anno_eset)

# first we store feature data in eset_final, defining PROBEID column
fData(eset_final)$PROBEID <- rownames(fData(eset_final))

# left_join anno_eset
fData(eset_final) <- left_join(fData(eset_final), anno_eset)

# restore row names
rownames(fData(eset_final)) <- fData(eset_final)$PROBEID
validObject(eset_final)


###### Differential Gene Expression
Groups <- as.character(pData(eset_final)$Groups)
Cells <- as.character(pData(eset_final)$Cells)
Treatment <- as.character(pData(eset_final)$Treatment)

# according to the original paper we look for the changes in transcription occur between treatments. we create two designe matrices one for each cell line

i_11 <- Groups[Cells == "KMS11"]
design_eset_11 <- model.matrix(~ 0 + Treatment[Cells == "KMS11"] + i_11)
colnames(design_eset_11)[1:2] <- c("Cfz", "Ctrl")
rownames(design_eset_11) <- i_11

head(design_eset_11)

i_34 <-  Groups[Cells == "KMS34"]
design_eset_34 <- model.matrix(~ 0 + Treatment[Cells == "KMS34"] + i_34)
colnames(design_eset_34)[1:2] <- c("Cfz", "Ctrl")
rownames(design_eset_34) <- i_34

head(design_eset_34)


contrast_matrix_11 <- makeContrasts(Cfz-Ctrl, levels = design_eset_11)
eset_fit_11 <- eBayes(contrasts.fit(lmFit(eset_final[,Cells =="KMS11"],
                                          design = design_eset_11),
                                    contrast_matrix_11))


contrast_matrix_34 <- makeContrasts(Cfz-Ctrl, levels = design_eset_34)
eset_fit_34 <- eBayes(contrasts.fit(lmFit(eset_final[,Cells =="KMS34"],
                                          design = design_eset_34),
                                    contrast_matrix_34))


## extracting results
table_11 <- topTable(eset_fit_11, number = Inf)
head(table_11)

hist(table_11$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Cfz vs Ctrl in KMS11", xlab = "p-values")




table_34 <- topTable(eset_fit_34, number = Inf)
head(table_34)

hist(table_34$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "Cfz vs Ctrl in KMS34", xlab = "p-values")



#### Multiple testing FDR

# we set the significance cutoff at 0.001, adjusted p.value represents FDR
nrow(subset(table_11, P.Value < 0.001)) # 1382 differentially expressed genes
nrow(subset(table_34, P.Value < 0.001)) # 1122 differentially expressed genes


### DE visualization
volcano_names <- ifelse(abs(eset_fit_11$coefficients)>=1,
                        eset_fit_11$genes$SYMBOL, NA)
volcanoplot(eset_fit_11, style = "p-value", highlight = 10, names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch = 16, cex = 0.35)


volcano_names <- ifelse(abs(eset_fit_34$coefficients)>=1,
                        eset_fit_34$genes$SYMBOL, NA)
volcanoplot(eset_fit_34, style = "p-value", highlight = 10, names = volcano_names,
            xlab = "Log2 Fold Change", ylab = NULL, pch = 16, cex = 0.35)

# we can get info of up or downregulated genes depicted on the plot, from genecards.org

##### Gene Ontology (GO) based enrichment analysis
# set FDR cutoff at 10%

DE_genes_11 <- subset(table_11, adj.P.Val < 0.1)$PROBEID
DE_genes_34 <- subset(table_34, adj.P.Val < 0.1)$PROBEID

# matching the background set of genes (genes that are similiar in expression to the DE genes)

back_genes_id11 <- genefinder(eset_final, as.character(DE_genes_11), method = "manhattan", scale = "none")
back_genes_id11 <- sapply(back_genes_id11, function(x)x$indices)

back_genes_id34 <- genefinder(eset_final, as.character(DE_genes_34), method = "manhattan", scale = "none")
back_genes_id34 <- sapply(back_genes_id34, function(x)x$indices)


back_genes11 <- featureNames(eset_final)[back_genes_id11]
back_genes11 <- setdiff(back_genes11, DE_genes_11)
intersect(back_genes11, DE_genes_11)
length(back_genes11)


back_genes34 <- featureNames(eset_final)[back_genes_id34]
back_genes34 <- setdiff(back_genes34, DE_genes_34)
intersect(back_genes34, DE_genes_34)
length(back_genes34)

### multidensity plot
# KMS11
multidensity(list(
  all = table_11[,"AveExpr"],
  fore = table_11[DE_genes_11, "AveExpr"],
  back = table_11[rownames(table_11) %in% back_genes11, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for KMS11-background-matching") # fore and back have similar shape, indicating a sensible background matching


# KMS34
multidensity(list(
  all = table_34[,"AveExpr"],
  fore = table_34[DE_genes_34, "AveExpr"],
  back = table_34[rownames(table_34) %in% back_genes34, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab = "mean expression",
  main = "DE genes for KMS34-background-matching") # fore and back have similar shape, indicating a sensible background matching


### topGO
gene_ID11 <- rownames(table_11)
in_universe <- gene_ID11 %in% c(DE_genes_11, back_genes11)
in_selection <- gene_ID11 %in% DE_genes_11

all_genes11 <- in_selection[in_universe]
all_genes11 <- factor(as.integer(in_selection[in_universe]))
names(all_genes11) <- gene_ID11[in_universe]



gene_ID34 <- rownames(table_34)
in_universe <- gene_ID34 %in% c(DE_genes_34, back_genes34)
in_selection <- gene_ID34 %in% DE_genes_34

all_genes34 <- in_selection[in_universe]
all_genes34 <- factor(as.integer(in_selection[in_universe]))
names(all_genes34) <- gene_ID34[in_universe]

## initialize topGO data set
top_GO_data_11 <- new("topGOdata", ontology = "BP", allGenes = all_genes11,
                      nodeSize = 10, annot = annFUN.db, affyLib = "hgu133plus2.db") # nodeSize=10 specifies a minimum size of GO category, categories with less than 10 genes are not included


top_GO_data_34 <- new("topGOdata", ontology = "BP", allGenes = all_genes34,
                      nodeSize = 10, annot = annFUN.db, affyLib = "hgu133plus2.db")


### testing with Fisher classic and elim algorithm
result11_topGO_elim <- runTest(top_GO_data_11, algorithm = "elim", statistic = "Fisher")
result11_topGO_classic <- runTest(top_GO_data_11, algorithm = "classic", statistic = "Fisher")



result34_topGO_elim <- runTest(top_GO_data_34, algorithm = "elim", statistic = "Fisher")
result34_topGO_classic <- runTest(top_GO_data_34, algorithm = "classic", statistic = "Fisher")


## inspect the results
res_top_GO_11 <- GenTable(top_GO_data_11, Fisher.elim = result11_topGO_elim, Fisher.classic = result11_topGO_classic,
                          orderBy = "Fisher.elim", topNodes = 20)
genes_top_GO_11 <- printGenes(top_GO_data_11, whichTerms = res_top_GO_11$GO.ID, chip = "hgu133plus2.db", geneCutOff = 1000)

res_top_GO_11$sig_genes <- sapply(genes_top_GO_11, function(x){
  str_c(paste0(x[x$"raw p-value" == 2, "Symbol.id"], ";"), collapse = "")
}) # significant genes are denoted with "2" in raw p-value


head(res_top_GO_11[,1:5])


res_top_GO_34 <- GenTable(top_GO_data_34, Fisher.elim = result34_topGO_elim, Fisher.classic = result34_topGO_classic,
                          orderBy = "Fisher.elim", topNodes = 20)
genes_top_GO_34 <- printGenes(top_GO_data_34, whichTerms = res_top_GO_34$GO.ID, chip = "hgu133plus2.db", geneCutOff = 1000)

res_top_GO_34$sig_genes <- sapply(genes_top_GO_34, function(x){
  str_c(paste0(x[x$"raw p-value" == 2, "Symbol.id"], ";"), collapse = "")
})

head(res_top_GO_34[,1:5])


## visualization of GO
showSigOfNodes(top_GO_data_11, score(result11_topGO_elim), firstSigNodes = 3, useInfo = "all")



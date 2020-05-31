#Title: M. capitata spawning gene expession
#last modified 20200530
#Copyright 2020 HM Putnam
#Use of this code must be accompanied by a citation to XXXX
#Data should not be used without permission from HM Putnam
#See Readme

rm(list=ls()) # removes all prior objects

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
#library("GSEABase")
library("goseq")
#library("GO.db")
library("UpSetR")
#library("reshape2")
library("tidyverse")
library("ggpubr")
library("apeglm")
library("multiClust")

##### Build Data Set of normalized expression and sample info #####
counts <- read.csv(file="Data/Mcap_gene_count_matrix.csv", header=T, sep=",", row.names=1) #Load expression matrix 
head(counts)
filt <- filterfun(pOverA(0.75,5)) #set filter values for PoverA, P percent of the samples have counts over A
tfil <- genefilter(counts, filt) #create filter for the counts data
keep <- counts[tfil,] #identify transcripts to keep by count filter
gn.keep <- rownames(keep) #identify transcript list
counts.5x <- as.matrix(counts[which(rownames(counts) %in% gn.keep),]) #data filtered in PoverA, P percent of the samples have counts over A
write.csv(counts.5x, file="Output/filtered_counts.csv")

storage.mode(counts.5x) = "integer" #store counts data as integer
sample.info <- read.csv(file="Data/sample_description.csv", header=T, sep=",", row.names=1) #load sample info

Mcap.annot <- read.csv(file="Data/Mcap-GO-KO-Kegg.tab", header=FALSE, sep="\t") #Load expression matrix f
colnames(Mcap.annot) <- c("Uniprot", "gene", "eval", "Prot.ID", "Rev", "Prot.Name.Long", "Prot.Name.Short", "Taxa", "Num", "GO1", "GO2", "GO3", "GO4", "GO.IDs","KEGG", "KEGG.Path")  
Mcap.annot <- Mcap.annot %>% 
  distinct(gene, .keep_all = TRUE)
Mcap.annot$gene <- gsub("augustus.", "", Mcap.annot$gene)
Mcap.annot$gene <- gsub(".t1", "", Mcap.annot$gene)

#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset
#it is also possible to retrieve the log2 fold changes, p values and adjusted p values of variables other than the last one in the design. While in this case, type is not biologically interesting as it indicates differences across sequencing protocol, for other hypothetical designs, such as ~genotype + condition + genotype:condition, we may actually be interested in the difference in baseline expression across genotype, which is not the last variable in the design.
#In any case, the contrast argument of the function results takes a character vector of length three: the name of the variable, the name of the factor level for the numerator of the log2 ratio, and the name of the factor level for the denominator. The contrast argument can also take other forms, as described in the help page for results and below
data <- DESeqDataSetFromMatrix(countData = counts.5x, colData = sample.info, design = ~ Spawn + Period + Spawn:Period) #create a DESeqDataSet object

##### Expression Visualization #####
rld <- rlog(data, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
head(assay(rld), 3) #view data
sampleDists <- dist(t(assay(rld))) #calculate distance matix
sampleDistMatrix <- as.matrix(sampleDists) #distance matrix
rownames(sampleDistMatrix) <- colnames(rld) #assign row names
colnames(sampleDistMatrix) <- NULL #assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pheatmap(sampleDistMatrix, #plot matrix of expression similarity
         clustering_distance_rows=sampleDists, #cluster rows
         clustering_distance_cols=sampleDists, #cluster columns
         col=colors) #set colors

plotPCA(rld, intgroup = c("Spawn", "Time")) #plot PCA of samples with all data

##### Differential Gene Expression Analysis #####
DEG.int <- DESeq(data) #run differential expression test by group using the wald test 
resultsNames(DEG.int)

SpawnCon.Res <- results(DEG.int, contrast=c("Spawn", "YES", "NO"))
head(SpawnCon.Res)
plotMA(SpawnCon.Res, ylim=c(-7,6))
sum(SpawnCon.Res$padj <0.05, na.rm=T)
SpawnCon.sig <- subset(SpawnCon.Res, padj<0.05,) #identify signficant pvalues with 5%FDR
SpawnCon.sig.list <- data[which(rownames(data) %in% rownames(SpawnCon.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(SpawnCon.sig.list), file="Output/Spawn_compare_DEG.csv")
SpawnCon.rsig <- rlog(SpawnCon.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
plotPCA(SpawnCon.rsig, intgroup = c("Spawn", "Time")) #Plot PCA of all samples for DEG only
mat <- assay(SpawnCon.rsig) #make an expression object
col.order <- c("T4.17","T4.1","T4.6",
               "T5.17","T5.1","T5.6",
               "T7.17","T7.1","T7.6",
               "T4.10","T4.16","T4.8",
               "T5.10","T5.16","T5.8",
               "T7.10","T7.16","T7.8")
mat <- mat[,col.order]
df <- as.data.frame(colData(SpawnCon.rsig)[,c("Spawn","Time")]) #make dataframe
df <- df[order(df$Spawn),]

pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, clustering_method = "average", scale="row",
                         clustering_distance_rows="euclidean", show_rownames =TRUE,fontsize_row = 4, cluster_cols=FALSE, cluster_rows=TRUE, 
                         show_colnames =TRUE) #plot heatmap of all DEG by group

#Time Comparisons
Time_18_20_Con.Res <- results(DEG.int, contrast=c("Period", "PreSpawn", "Set"))
head(Time_18_20_Con.Res)
plotMA(Time_18_20_Con.Res, ylim=c(-7,4))
sum(Time_18_20_Con.Res$padj <0.05, na.rm=T)
Time_18_20_Con.sig <- subset(Time_18_20_Con.Res, padj<0.05,) #identify signficant pvalues with 5%FDR
Time_18_20_Con.sig.list <- data[which(rownames(data) %in% rownames(Time_18_20_Con.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(Time_18_20_Con.sig.list), file="Output/Spawn_compare_DEG_18_20.csv")
Time_18_20_.rsig <- rlog(Time_18_20_Con.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
plotPCA(Time_18_20_.rsig, intgroup = c("Time")) #Plot PCA of all samples for DEG only
mat <- assay(Time_18_20_.rsig) #make an expression object
mat <- mat[,col.order]
df <- as.data.frame(colData(Time_18_20_.rsig)[,c("Spawn","Time")]) #make dataframe
df <- df[order(df$Spawn),]

pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, clustering_method = "average", scale="row",
         clustering_distance_rows="euclidean", show_rownames =TRUE,fontsize_row = 4, cluster_cols=FALSE, cluster_rows=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

Time_20_00_Con.Res <- results(DEG.int, contrast=c("Period", "Set", "PostSpawn"))
head(Time_20_00_Con.Res)
plotMA(Time_20_00_Con.Res, ylim=c(-4,7))
sum(Time_20_00_Con.Res$padj <0.05, na.rm=T)
Time_20_00_Con.sig <- subset(Time_20_00_Con.Res, padj<0.05,) #identify signficant pvalues with 5%FDR
Time_20_00_Con.sig.list <- data[which(rownames(data) %in% rownames(Time_20_00_Con.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(Time_20_00_Con.sig.list), file="Output/Spawn_compare_DEG_20_00.csv")
Time_20_00_.rsig <- rlog(Time_20_00_Con.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
plotPCA(Time_20_00_.rsig, intgroup = c("Time")) #Plot PCA of all samples for DEG only
mat <- assay(Time_20_00_.rsig) #make an expression object
mat <- mat[,col.order]
df <- as.data.frame(colData(Time_20_00_.rsig)[,c("Spawn","Time")]) #make dataframe
df <- df[order(df$Spawn),]

pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, clustering_method = "average", scale="row",
         clustering_distance_rows="euclidean", show_rownames =TRUE,fontsize_row = 4, cluster_cols=FALSE, cluster_rows=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

Time_18_00_Con.Res <- results(DEG.int, contrast=c("Period", "PreSpawn", "PostSpawn"))
head(Time_18_00_Con.Res)
plotMA(Time_18_00_Con.Res, ylim=c(-5,6))
sum(Time_18_00_Con.Res$padj <0.05, na.rm=T)
Time_18_00_Con.sig <- subset(Time_18_00_Con.Res, padj<0.05,) #identify signficant pvalues with 5%FDR
Time_18_00_Con.sig.list <- data[which(rownames(data) %in% rownames(Time_18_00_Con.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(Time_18_00_Con.sig.list), file="Output/Spawn_compare_DEG_18_00.csv")
Time_18_00_.rsig <- rlog(Time_18_00_Con.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
plotPCA(Time_18_00_.rsig, intgroup = c("Time")) #Plot PCA of all samples for DEG only
mat <- assay(Time_18_00_.rsig) #make an expression object
mat <- mat[,col.order]
df <- as.data.frame(colData(Time_18_00_.rsig)[,c("Spawn","Time")]) #make dataframe
df <- df[order(df$Spawn),]

pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, clustering_method = "average", scale="row",
         clustering_distance_rows="euclidean", show_rownames =TRUE,fontsize_row = 4, cluster_cols=FALSE, cluster_rows=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

resultsNames(DEG.int)
# "Intercept"                    "Spawn_YES_vs_NO"             
# "Period_PreSpawn_vs_PostSpawn" "Period_Set_vs_PostSpawn"     
# "SpawnYES.PeriodPreSpawn"      "SpawnYES.PeriodSet" 

Int_Con.Res <- results(DEG.int, contrast=list(c("SpawnYES.PeriodPreSpawn", "SpawnYES.PeriodSet")))
plotMA(Int_Con.Res, ylim=c(-8,12))
sum(Int_Con.Res$padj <0.05, na.rm=T)
Int_Con.sig <- subset(Int_Con.Res, padj<0.05,) #identify signficant pvalues with 5%FDR
Int_Con.sig.list <- data[which(rownames(data) %in% rownames(Int_Con.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(Int_Con.sig.list), file="Output/TxS_Interaction_compare_DEG.csv")
Int.rsig <- rlog(Int_Con.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
plotPCA(Int.rsig, intgroup = c("Time")) #Plot PCA of all samples for DEG only
mat <- assay(Int.rsig) #make an expression object
mat <- mat[,col.order]
df <- as.data.frame(colData(Int.rsig)[,c("Spawn","Time")]) #make dataframe
df <- df[order(df$Spawn),]

pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, clustering_method = "average", scale="row",
         clustering_distance_rows="euclidean", show_rownames =TRUE,fontsize_row = 4, cluster_cols=FALSE, cluster_rows=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

#Unique Genes Heatmap
Time1820.DEGs <- as.data.frame(row.names(Time_18_20_Con.sig.list))
colnames(Time1820.DEGs) <- "DEGs" 
Time1800.DEGs <- as.data.frame(row.names(Time_18_00_Con.sig.list))
colnames(Time1800.DEGs) <- "DEGs" 
Time2000.DEGs <- as.data.frame(row.names(Time_20_00_Con.sig.list))
colnames(Time2000.DEGs) <- "DEGs" 
Spawn.DEGs<- as.data.frame(row.names(SpawnCon.sig.list))
colnames(Spawn.DEGs) <- "DEGs"
#Interact.DEGs <- as.data.frame(row.names(Int_Con.sig.list))
#colnames(Interact.DEGs) <- "DEGs"

DEGS.ALL <- rbind(Spawn.DEGs,Time1820.DEGs,Time1800.DEGs,Time2000.DEGs)
DEGS.ALL <- unique(DEGS.ALL)
length(t(unique(DEGS.ALL)))

Unique.sig.list <- data[which(rownames(data) %in% DEGS.ALL$DEGs),] #subset list of sig transcripts from original count data
Unique.rsig <- rlog(counts(Unique.sig.list), blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
write.csv(counts(Unique.sig.list), file="Output/Unique_DEG.csv")

PCA.plot <- plotPCA(Unique.rsig, intgroup = c("Spawn", "Period")) #Plot PCA of all samples for DEG only
PCA.plot #view plot
PC.info <-PCA.plot$data #extract plotting data
dev.off()
jpeg(file="Output/Unique_PCA.DEG.jpg")
plot(PC.info$PC1, PC.info$PC2, xlab="PC1 50%", ylab="PC2 25%", pch = c(15, 16)[as.numeric(sample.info$Spawn)], col=c("lightpink2","steelblue1","yellow3")[sample.info$Period], cex=1.3)
legend(x="right", 
       bty="n",
       legend = c("PostSpawn", "PreSpawn", "Set", "Not Spawning", "Spawning"),
       text.col = c("lightpink2","steelblue1","yellow3", "black", "black"),
       pch = c(15, 15, 15, 15, 16),
       col = c("white","white","white", "black", "black"),
       cex=1)
dev.off()


df <- as.data.frame(colData(Unique.rsig)[,c("Spawn","Time")]) #make dataframe
col.order <- c("T4.17","T4.1","T4.6",
               "T5.17","T5.1","T5.6",
               "T7.17","T7.1","T7.6",
               "T4.10","T4.16","T4.8",
               "T5.10","T5.16","T5.8",
               "T7.10","T7.16","T7.8")

Unique.DEG.annot <- as.data.frame(counts(Unique.sig.list))
Unique.DEG.annot$gene <- row.names(Unique.DEG.annot)
Unique.DEG.annot <- left_join(Unique.DEG.annot, Mcap.annot)
write.csv(Unique.DEG.annot, file="Output/Unique_DEG_annotated.csv")

Unique.DEG.annot$ann.row <- paste0(Unique.DEG.annot$gene," ", Unique.DEG.annot$Prot.ID)
rownames(Unique.DEG.annot) <- Unique.DEG.annot$ann.row
mat <- as.matrix(Unique.DEG.annot[,1:18]) #make a matrix
mat <- mat[,col.order]

out <- pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, scale="row",
                show_rownames =T, fontsize_row = 4, cluster_cols = FALSE,
                show_colnames =F) #plot heatmap of all DEG by group

dev.off()
pdf(file="Output/Unique_Heatmap.DEG.Annotated.pdf") #save file
pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, scale="row",
         show_rownames =T, fontsize_row = 4, cluster_cols = FALSE,
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()

#K means clustering Elbow Method
set.seed(123) #Specify which set of random numbers to use. This will make sure that the same set of random numbers is generated each timea random set of numbers is called. 
wss <- function(k) {kmeans(Unique.rsig, k, iter.max = 100, nstart = 25 )$tot.withins} # function to compute total within-cluster sum of square 
k.values <- 1:15# Compute and plot wss for k = 1 to k = 15
wss_values <- map_dbl(k.values, wss) # extract wss for 2-15 clusters
pdf(file="Output/Mcap_elbow.pdf")
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters (K)",
     ylab="Total within-clusters sum of squares")
dev.off()

##### Extract the clusters
clusts <- as.data.frame(sort(cutree(out$tree_row, k=6)))
clusts$gene <- rownames(clusts)
colnames(clusts) <- c("Cluster", "gene")











##### Analysis by Group #####
sample.info$group <- paste0(sample.info$Period,"_", sample.info$Spawn)

data2 <- DESeqDataSetFromMatrix(countData = counts.5x, colData = sample.info, design = ~ group) #create a DESeqDataSet object
DEG.int <- DESeq(data2) #run differential expression test by group using the wald test 
results(DEG.int)
DEG.int.res <- results(DEG.int) #save DE results
head(DEG.int.res)
plotMA(DEG.int.res, ylim=c(-7,6))
sum(DEG.int.res$padj <0.05, na.rm=T)
DEG.sig <- subset(DEG.int.res, padj<0.05,) #identify signficant pvalues with 5%FDR
plotMA(DEG.sig, ylim=c(-7,6))
DEG.sig.list <- data2[which(rownames(data2) %in% rownames(DEG.sig)),] #subset list of sig transcripts from original count data
write.csv(counts(DEG.sig.list), file="Output/Group_compare_DEG.csv")
DEG.rsig <- rlog(DEG.sig.list, blind=FALSE) #apply a regularized log transformation to minimize effects of small counts and normalize wrt library size
PCA.plot <- plotPCA(DEG.rsig, intgroup = c("Spawn", "Period")) #Plot PCA of all samples for DEG only
PCA.plot #view plot

PC.info <-PCA.plot$data #extract plotting data
pdf(file="Output/PCA.DEG.pdf")
plot(PC.info$PC1, PC.info$PC2, xlim=c(-12,12), ylim=c(-9, 9), xlab="PC1 44%", ylab="PC2 28%", col = c("lightpink2", "steelblue1","yellow3")[as.numeric(PC.info$Time)], pch=c(16, 17)[as.numeric(PC.info$Spawn)], cex=1.3)
legend(x="topleft", 
       bty="n",
       legend = c("00:00", "18:00", "20:00", "Not Spawning", "Spawning"),
       text.col = c("lightpink2","steelblue1","yellow3", "black", "black"),
       pch = c(15, 15, 15, 16, 17),
       col = c("white","white","white", "black", "black"),
       cex=1)
dev.off()

mat <- assay(DEG.rsig) #make an expression object
col.order <- c("T4.17","T4.1","T4.6",
               "T5.17","T5.1","T5.6",
               "T7.17","T7.1","T7.6",
               "T4.10","T4.16","T4.8",
               "T5.10","T5.16","T5.8",
               "T7.10","T7.16","T7.8")
mat <- mat[,col.order]
df <- as.data.frame(colData(DEG.rsig)[,c("Spawn","Period")]) #make dataframe
df <- df[order(df$Spawn),]

DEG.Heat.res <- pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, clustering_method = "average", scale="row",
         clustering_distance_rows="euclidean", show_rownames =FALSE, cluster_cols=FALSE, cluster_rows=TRUE,
         show_colnames =TRUE) #plot heatmap of all DEG by group

#Pairwise results
resultsNames(DEG.int)
#"Intercept"                           "group_PostSpawn_YES_vs_PostSpawn_NO"
# "group_PreSpawn_NO_vs_PostSpawn_NO"   "group_PreSpawn_YES_vs_PostSpawn_NO" 
# "group_Set_NO_vs_PostSpawn_NO"        "group_Set_YES_vs_PostSpawn_NO"
#Int_Con.Res <- results(DEG.int, contrast=c("group", "Set_NO_vs_PostSpawn_NO","Set_YES_vs_PostSpawn_NO"))






##### Examine Transcript Cluster Patterns #####
#check for optimal number of clusters for reporting expression patterns
#K means clustering Elbow Method
set.seed(123) #Specify which set of random numbers to use. This will make sure that the same set of random numbers is generated each timea random set of numbers is called. 
wss <- function(k) {kmeans(trans.data, k, iter.max = 100, nstart = 25 )$tot.withins} # function to compute total within-cluster sum of square 
k.values <- 1:15# Compute and plot wss for k = 1 to k = 15
wss_values <- map_dbl(k.values, wss) # extract wss for 2-15 clusters
pdf(file="Output/Mcap_elbow.pdf")
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters (K)",
     ylab="Total within-clusters sum of squares")
dev.off()

#Permutational analysis to test cluster against a null model (i.e., chance along)

#set optimal cluster number based on results above
knum <-10
DEG.clust <- cbind(mat, cluster = cutree(DEG.Heat.res$tree_row, k = knum))[DEG.Heat.res$tree_row[["order"]]]
DEG.clust <- as.data.frame(DEG.clust)
cluster.order <- cbind(mat[c(DEG.Heat.res$tree_row[["order"]]),DEG.Heat.res$tree_col[["order"]]],cluster=cutree(DEG.Heat.res$tree_row, k=knum)[DEG.Heat.res$tree_row[["order"]]])

cluster.order <- as.data.frame(cluster.order)
DEG.clust1 <- cluster.order[cluster.order$cluster == '1',]
DEG.clust1 <- DEG.clust1[,-19] 
DEG.clust2 <- cluster.order[cluster.order$cluster == '2',]
DEG.clust2 <- DEG.clust2[,-19] 
DEG.clust3 <- cluster.order[cluster.order$cluster == '3',]
DEG.clust3 <- DEG.clust3[,-19] 
DEG.clust4 <- cluster.order[cluster.order$cluster == '4',]
DEG.clust4 <- DEG.clust4[,-19] 
DEG.clust5 <- cluster.order[cluster.order$cluster == '5',]
DEG.clust5 <- DEG.clust5[,-19] 
DEG.clust6 <- cluster.order[cluster.order$cluster == '6',]
DEG.clust6 <- DEG.clust6[,-19] 
DEG.clust7 <- cluster.order[cluster.order$cluster == '7',]
DEG.clust7 <- DEG.clust7[,-19] 
DEG.clust8 <- cluster.order[cluster.order$cluster == '8',]
DEG.clust8 <- DEG.clust8[,-19] 
DEG.clust9 <- cluster.order[cluster.order$cluster == '9',]
DEG.clust9 <- DEG.clust9[,-19] 
DEG.clust10 <- cluster.order[cluster.order$cluster == '10',]
DEG.clust10 <- DEG.clust10[,-19]
DEG.clust11 <- cluster.order[cluster.order$cluster == '11',]
DEG.clust11 <- DEG.clust11[,-19] 
DEG.clust12 <- cluster.order[cluster.order$cluster == '12',]
DEG.clust12 <- DEG.clust12[,-19] 
DEG.clust13 <- cluster.order[cluster.order$cluster == '13',]
DEG.clust13 <- DEG.clust13[,-19] 
DEG.clust14 <- cluster.order[cluster.order$cluster == '14',]
DEG.clust14 <- DEG.clust14[,-19] 
DEG.clust15 <- cluster.order[cluster.order$cluster == '15',]
DEG.clust15 <- DEG.clust15[,-19] 
DEG.clust16 <- cluster.order[cluster.order$cluster == '16',]
DEG.clust16 <- DEG.clust16[,-19] 
DEG.clust17 <- cluster.order[cluster.order$cluster == '17',]
DEG.clust17 <- DEG.clust17[,-19] 
DEG.clust18 <- cluster.order[cluster.order$cluster == '18',]
DEG.clust18 <- DEG.clust18[,-19] 
DEG.clust19 <- cluster.order[cluster.order$cluster == '19',]
DEG.clust19 <- DEG.clust19[,-19] 
DEG.clust20 <- cluster.order[cluster.order$cluster == '20',]
DEG.clust20 <- DEG.clust20[,-19] 
DEG.clust21 <- cluster.order[cluster.order$cluster == '21',]
DEG.clust21 <- DEG.clust21[,-19] 
DEG.clust22 <- cluster.order[cluster.order$cluster == '22',]
DEG.clust22 <- DEG.clust22[,-19] 
DEG.clust23 <- cluster.order[cluster.order$cluster == '23',]
DEG.clust23 <- DEG.clust23[,-19] 
DEG.clust24 <- cluster.order[cluster.order$cluster == '24',]
DEG.clust24 <- DEG.clust24[,-19] 
DEG.clust25 <- cluster.order[cluster.order$cluster == '25',]
DEG.clust25 <- DEG.clust25[,-19] 

sample.info <- sample.info[order(row.names(sample.info)), ]
CL1 <- t(DEG.clust1)
CL1 <- CL1[ order(row.names(CL1)), ]
CL1 <- cbind(CL1, sample.info)

CL2 <- t(DEG.clust2)
CL2 <- CL2[ order(row.names(CL2)), ]
CL2 <- cbind(CL2, sample.info)

CL3 <- t(DEG.clust3)
CL3 <- CL3[ order(row.names(CL3)), ]
CL3 <- cbind(CL3, sample.info)

CL4 <- t(DEG.clust4)
CL4 <- CL4[ order(row.names(CL4)), ]
CL4 <- cbind(CL4, sample.info)

CL5 <- t(DEG.clust5)
CL5 <- CL5[ order(row.names(CL5)), ]
CL5 <- cbind(CL5, sample.info)

CL6 <- t(DEG.clust6)
CL6 <- CL6[ order(row.names(CL6)), ]
CL6 <- cbind(CL6, sample.info)

CL7 <- t(DEG.clust7)
CL7 <- CL7[ order(row.names(CL7)), ]
CL7 <- cbind(CL7, sample.info)

CL8 <- t(DEG.clust8)
CL8 <- CL8[ order(row.names(CL8)), ]
CL8 <- cbind(CL8, sample.info)

CL9 <- t(DEG.clust9)
CL9 <- CL9[ order(row.names(CL9)), ]
CL9 <- cbind(CL9, sample.info)

CL10 <- t(DEG.clust10)
CL10 <- CL10[ order(row.names(CL10)), ]
CL10 <- cbind(CL10, sample.info)

CL11 <- t(DEG.clust11)
CL11 <- CL11[ order(row.names(CL11)), ]
CL11 <- cbind(CL11, sample.info)

CL12 <- t(DEG.clust12)
CL12 <- CL12[ order(row.names(CL12)), ]
CL12 <- cbind(CL12, sample.info)

CL13 <- t(DEG.clust13)
CL13 <- CL13[ order(row.names(CL13)), ]
CL13 <- cbind(CL13, sample.info)

CL14 <- t(DEG.clust14)
CL14 <- CL14[ order(row.names(CL14)), ]
CL14 <- cbind(CL14, sample.info)

CL15 <- t(DEG.clust15)
CL15 <- CL15[ order(row.names(CL15)), ]
CL15 <- cbind(CL15, sample.info)

CL16 <- t(DEG.clust16)
CL16 <- CL16[ order(row.names(CL16)), ]
CL16 <- cbind(CL16, sample.info)

CL17 <- t(DEG.clust17)
CL17 <- CL17[ order(row.names(CL17)), ]
CL17 <- cbind(CL17, sample.info)

CL18 <- t(DEG.clust18)
CL18 <- CL18[ order(row.names(CL18)), ]
CL18 <- cbind(CL18, sample.info)

CL19 <- t(DEG.clust19)
CL19 <- CL19[ order(row.names(CL19)), ]
CL19 <- cbind(CL19, sample.info)

CL20 <- t(DEG.clust20)
CL20 <- CL20[ order(row.names(CL20)), ]
CL20 <- cbind(CL20, sample.info)

CL21 <- t(DEG.clust21)
CL21 <- CL21[ order(row.names(CL21)), ]
CL21 <- cbind(CL21, sample.info)

CL22 <- t(DEG.clust22)
CL22 <- CL22[ order(row.names(CL22)), ]
CL22 <- cbind(CL22, sample.info)

CL23 <- t(DEG.clust23)
CL23 <- CL23[ order(row.names(CL23)), ]
CL23 <- cbind(CL23, sample.info)

CL24 <- t(DEG.clust24)
CL24 <- CL24[ order(row.names(CL24)), ]
CL24 <- cbind(CL24, sample.info)

CL25 <- t(DEG.clust25)
CL25 <- CL25[ order(row.names(CL25)), ]
CL25 <- cbind(CL25, sample.info)

cluster1 <- melt(CL1, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster1$Cluster <- "Cluster1"
cluster2 <- melt(CL2, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster2$Cluster <- "Cluster2"
cluster3 <- melt(CL3, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp") 
cluster3$Cluster <- "Cluster3"
cluster4 <- melt(CL4, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster4$Cluster <- "Cluster4"
cluster5 <- melt(CL5, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster5$Cluster <- "Cluster5"
cluster6 <- melt(CL6, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster6$Cluster <- "Cluster6"
cluster7 <- melt(CL7, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster7$Cluster <- "Cluster7"
cluster8 <- melt(CL8, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster8$Cluster <- "Cluster8"
cluster9 <- melt(CL9, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster9$Cluster <- "Cluster9"
cluster10 <- melt(CL10, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster10$Cluster <- "Cluster10"

cluster11 <- melt(CL11, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster11$Cluster <- "Cluster11"
cluster12 <- melt(CL12, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster12$Cluster <- "Cluster12"
cluster13 <- melt(CL13, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp") 
cluster13$Cluster <- "Cluster13"
cluster14 <- melt(CL14, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster14$Cluster <- "Cluster14"
cluster15 <- melt(CL15, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster15$Cluster <- "Cluster15"
cluster16 <- melt(CL16, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster16$Cluster <- "Cluster16"
cluster17 <- melt(CL17, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster17$Cluster <- "Cluster17"
cluster18 <- melt(CL18, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster18$Cluster <- "Cluster18"
cluster19 <- melt(CL19, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster19$Cluster <- "Cluster19"
cluster20 <- melt(CL20, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster20$Cluster <- "Cluster20"

cluster21 <- melt(CL21, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster21$Cluster <- "Cluster21"
cluster22 <- melt(CL22, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster22$Cluster <- "Cluster22"
cluster23 <- melt(CL23, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp") 
cluster23$Cluster <- "Cluster23"
cluster24 <- melt(CL24, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster24$Cluster <- "Cluster24"
cluster25 <- melt(CL25, id=c("Spawn", "Time","TimePoint", "Stage", "Rep", "group"), value.name="Rel.Exp")
cluster25$Cluster <- "Cluster25"

C1.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster1, mean, na.rm = TRUE)
C2.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster2, mean, na.rm = TRUE)
C3.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster3, mean, na.rm = TRUE)
C4.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster4, mean, na.rm = TRUE)
C5.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster5, mean, na.rm = TRUE)
C6.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster6, mean, na.rm = TRUE)
C7.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster7, mean, na.rm = TRUE)
C8.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster8, mean, na.rm = TRUE)
C9.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster9, mean, na.rm = TRUE)
C10.means <- aggregate(Rel.Exp ~ Spawn + TimePoint, data = cluster10, mean, na.rm = TRUE)

#Want to plot expression plots for each cluster with all data points and lines of mean of cluster
#YES="gray", NO="black"
dev.off()
CLP1 <- ggplot(cluster1, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster1")+
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP2 <- ggplot(cluster2, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster2")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP3 <- ggplot(cluster3, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster3")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP4 <- ggplot(cluster4, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster4")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP5 <- ggplot(cluster5, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster5")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP6 <- ggplot(cluster6, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster6")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP7 <- ggplot(cluster7, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster7")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP8 <- ggplot(cluster8, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster8")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP9 <- ggplot(cluster9, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster9")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP10 <- ggplot(cluster10, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster10")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP11 <- ggplot(cluster11, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster11")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP12 <- ggplot(cluster12, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster12")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP13 <- ggplot(cluster13, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster13")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP14 <- ggplot(cluster14, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster14")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP15 <- ggplot(cluster15, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster15")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP16 <- ggplot(cluster16, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster16")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP17 <- ggplot(cluster17, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster17")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP18 <- ggplot(cluster18, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster18")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP19 <- ggplot(cluster19, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster19")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP20 <- ggplot(cluster20, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster20")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP21 <- ggplot(cluster21, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster21")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP22 <- ggplot(cluster22, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster22")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP23 <- ggplot(cluster23, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster23")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP24 <- ggplot(cluster24, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster24")+
  theme_classic() +
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

CLP25 <- ggplot(cluster25, aes(x=TimePoint, y=Rel.Exp, fill=Spawn)) +
  geom_boxplot() +
  scale_fill_manual(values=c("black", "gray")) +
  labs(title="Cluster25")+
  theme_classic() +
  theme(legend.position = "none")

figure <- ggarrange(CLP2, CLP4, CLP3, CLP5, CLP9,
                    CLP14, CLP12, CLP20, CLP8, CLP7,
                    CLP6, CLP1, CLP13, CLP15, CLP17,
                    CLP19, CLP16, CLP18, CLP11, CLP10,
                    CLP25, CLP22, CLP23, CLP21, CLP24,
                    ncol = 1, nrow = 25)
figure

ggsave("Output/ggclusters.pdf", device = "pdf", 
       width = 2,
       height = 25,
       units = c("in"))

##### HEATMAP OUTPUT #####

row.annots <- as.data.frame(cluster.order$cluster)
rownames(row.annots) <- rownames(cluster.order)

#set annotation colors on heatmap factors
ann_colors <- list(Spawn = c(YES="gray", NO="black"), Time = c("18:00"="steelblue1", "20:00"="yellow3", "0:00"="lightpink2"))
#"lightpink2", "yellow3", "steelblue1"

dev.off()
pdf(file="Output/DEG_Heatmap.pdf", width=4, height=16)
pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, annotation_colors =ann_colors, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =FALSE,  cutree_rows=25, cluster_cols=F,
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()

pdf(file="Output/Annotated_DEG_Heatmap.pdf", width=4, height=16)
pheatmap(mat, color = colorRampPalette(c("darkblue", "cyan", "white", "yellow", "red"))(100), annotation_col=df, annotation_colors =ann_colors, clustering_method = "average", 
         clustering_distance_rows="euclidean", show_rownames =TRUE, fontsize_row = 4,  cutree_rows=25, cluster_cols=F,
         show_colnames =F) #plot heatmap of all DEG by group
dev.off()


# ##### Generate Annotated Transcript Expression Cluster Table #####
# 
# All.Clusters <- rbind(unique.C1, unique.C2, unique.C3, unique.C4, unique.C5, unique.C6) #bind clusters
# All.Clusters <- All.Clusters[,c(1,21, 23, 26, 30, 31, 32, 46)] #keep selected columns
# write.csv(All.Clusters, file="Output/All.Clusters_unique_annotated_full.csv")
# All.Clusters$gene_ontology_blast2 <- All.Clusters$gene_ontology_blast #add another copy of the GO information
# All.Clusters <- separate(data = All.Clusters, col = sprot_Top_BLASTX_hit, into = c("Protein.Name", "Full.Name"), sep = "\\^", extra="merge") #split columns to get protein ID
# All.Clusters <- separate(data = All.Clusters, col = Full.Name, into = c( "Full.Name", "Name"), sep = "\\=", extra="merge") #split columns to get protein name
# All.Clusters <- separate(data = All.Clusters, col = Name, into = c( "Full.Name"), sep = "\\;", extra="drop") #split columns to get protein name
# All.Clusters <- separate(data = All.Clusters, col = gene_ontology_blast, into = c("CC", "MF"), sep = "GO\\:[0-9]*\\^molecular_function\\^", extra="merge") #split columns to get CC
# All.Clusters <- separate(data = All.Clusters, col = gene_ontology_blast2, into = c("Funct", "BP"), sep = "GO\\:[0-9]*\\^biological_process\\^", extra="merge") #split columns to get BP
# All.Clusters <- separate(data = All.Clusters, col = CC, into = c("CC"), sep = "GO\\:[0-9]*\\^biological_process\\^", extra="drop") #split columns to get CC
# All.Clusters <- separate(data = All.Clusters, col = MF, into = c("MF"), sep = "GO\\:[0-9]*\\^biological_process\\^", extra="drop") #split columns to get MF
# All.Clusters <- subset(All.Clusters, select=-c(Funct)) #remove unwanted text
# 
# All.Clusters$CC <- gsub("GO\\:[0-9]*\\^cellular_component\\^", "", All.Clusters$CC) #remove GO ids
# All.Clusters$CC <- gsub("\\`",", ", All.Clusters$CC) #remove extra symbols
# All.Clusters$CC
# 
# All.Clusters$MF <- gsub("GO\\:[0-9]*\\^molecular_function\\^", "", All.Clusters$MF)#remove GO ids
# All.Clusters$MF <- gsub("\\`",", ", All.Clusters$MF)#remove extra symbols
# All.Clusters$MF
# 
# All.Clusters$BP <- gsub("GO\\:[0-9]*\\^biological_process\\^", "", All.Clusters$BP)#remove GO ids
# All.Clusters$BP <- gsub("\\`",", ", All.Clusters$BP)#remove extra symbols
# All.Clusters$BP
# 
# All.Clusters <- cbind(All.Clusters$Protein.Name, All.Clusters$Full.Name, All.Clusters$CC, All.Clusters$MF, All.Clusters$BP, All.Clusters$Cluster, All.Clusters$Transcript.ID) #choose columns for table
# write.csv(All.Clusters, file="Output/All.Clusters_unique_annotated_Table2.csv") #write table to file
# 
# ##### Pairwise Contrasts of Spawning and Non-Spawnign Corals at each Time Point ##### 
# #pairwise comparisons of DEGs for SP and NSP corals at each collection point
# DEG.18 <- results(DEG.int, contrast=c("group", "YES18:00", "NO18:00")) #contrasts between SP NSP at each time
# DEG.18 #view results
# sig.num.18 <- sum(DEG.18$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.18 <- subset(DEG.18, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.18 <- data[which(rownames(data) %in% rownames(sig.18)),] #subset list of sig transcripts from original count data
# 
# DE18 <- as.data.frame(sig.18)
# 
# DEG.20 <- results(DEG.int, contrast=c("group", "YES20:00", "NO20:00")) #contrasts between SP NSP at each time
# DEG.20 #view results
# sig.num.20 <- sum(DEG.20$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.20 <- subset(DEG.20, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.20 <- data[which(rownames(data) %in% rownames(sig.20)),] #subset list of sig transcripts from original count data
# 
# DE20 <- as.data.frame(sig.20)
# 
# DEG.00 <- results(DEG.int, contrast=c("group", "YES0:00", "NO0:00")) #contrasts between SP NSP at each time
# DEG.00 #view results
# sig.num.00 <- sum(DEG.00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.00 <- subset(DEG.00, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.00 <- data[which(rownames(data) %in% rownames(sig.00)),] #subset list of sig transcripts from original count data
# 
# DE00 <- as.data.frame(sig.00)
# 
# ##### Visualize time points  in Venn Diagram or Intersection Plot #####
# L18 <- data.frame(rownames(sig.list.18)) #identify DE transcript list from 18:00
# L18$L18.Count <- 1 #assign all transcripts count of 1
# colnames(L18) <- c("Transcript.ID", "L18.Count") #name columns
# L20 <- data.frame(rownames(sig.list.20)) #identify DE transcript list from 20:00
# L20$L20.Count <- 1 #assign all transcripts count of 1
# colnames(L20) <- c("Transcript.ID", "L20.Count") #name columns
# L00 <- data.frame(rownames(sig.list.00)) #identify DE transcript list from 00:00
# L00$L00.Count <- 1 #assign all transcripts count of 1
# colnames(L00) <- c("Transcript.ID", "L00.Count") #name columns
# L18.20 <- merge(L18, L20, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
# L18.20.00 <- merge(L18.20, L00, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
# L18.20.00$L18.Count[is.na(L18.20.00$L18.Count)] <- 0 #assign NA=0 count i.e., not DE
# L18.20.00$L20.Count[is.na(L18.20.00$L20.Count)] <- 0 #assign NA=0 count i.e., not DE
# L18.20.00$L00.Count[is.na(L18.20.00$L00.Count)] <- 0 #assign NA=0 count i.e., not DE
# 
# # Plot Venn Diagram of shared DE between time points
# vennData<-L18.20.00[,2:4] #set venn data to DE counts only
# a <- vennCounts(vennData) # compute classification counts
# colnames(a) <- c('18:00', '20:00', '00:00', "Counts") #set catagory names
# #pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.pdf")
# #vennDiagram(a, main='DEG between Spawning and NonSpawning Corals', circle.col=c("yellow", "red", "green")) #draw venn diagram
# #dev.off()
# 
# #Group Intersections
# colnames(vennData) <- c("18:00","20:00", "00:00") #set catagory names
# pdf(file="Output/SP_NSP_DEG_Intersections.pdf", useDingbats = TRUE)
# upset(vennData, sets = c('18:00','20:00', '00:00'), sets.bar.color = c("lightpink2", "yellow3", "steelblue1"), main.bar.color = c("darkgrey", "black", "black", "black", "steelblue1", "yellow3", "lightpink2"), order.by = "degree")
# dev.off()
# 
# #Look for match between all lists
# shared.all <- Reduce(intersect, list(rownames(sig.list.18),rownames(sig.list.20),rownames(sig.list.00))) #idenitfy shared DE transcripts between each time point 
# shared.all #view results
# shared.all <- as.data.frame(shared.all)
# colnames(shared.all)[1] <- "Transcript.ID"
# write.csv(shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_18_20_00.csv")
# 
# ##### DEG lists between spawning and not spawning corals unique to time points and shared between time points ##### 
# intersections <- subset(L18.20.00, L18.Count==1 & L20.Count==1 & L00.Count==1) #list DEG found in all time points
# 
# #DEG between spawning and not spawning corals Unique to each time point
# unique.18 <- subset(L18.20.00, L18.Count==1 & L20.Count==0 & L00.Count==0) #list genes unique to 18:00
# unique.20 <- subset(L18.20.00, L18.Count==0 & L20.Count==1 & L00.Count==0) #list genes unique to 20:00
# unique.00 <- subset(L18.20.00, L18.Count==0 & L20.Count==0 & L00.Count==1) #list genes unique to 00:00
# 
# #DEG between spawning and not spawning corals Shared unique between pairs
# unique.18.20 <- subset(L18.20.00, L18.Count==1 & L20.Count==1 & L00.Count==0) #list genes common to 18:00 and 20:00
# unique.18.00 <- subset(L18.20.00, L18.Count==1 & L20.Count==0 & L00.Count==1) #list genes common to 18:00 and 00:00
# unique.20.00 <- subset(L18.20.00, L18.Count==0 & L20.Count==1 & L00.Count==1) #list genes common to 20:00 and 00:00
# 
# #combine ID and rlog normalized counts of genes found in all time points
# norm.sig.counts <- as.data.frame(mat)
# norm.sig.counts$Transcript.ID <- rownames(norm.sig.counts)
# inters <- merge(intersections, norm.sig.counts, by="Transcript.ID") 
# 
# # Combine with DE list with annotation data
# anot.inters <- merge(inters,  annot, by="Gene.ID") #merge sig intersection list and annotations
# # anot.inters$KG <- sapply(strsplit(as.character(anot.inters$Kegg), split="\\`"), "[", 1)
# # anot.inters$KG <- gsub(".*:","",anot.inters$KG)
# # anot.inters$KO <- sapply(strsplit(as.character(anot.inters$Kegg), split="\\`"), "[", 2)
# # anot.inters$KO <- gsub(".*:","",anot.inters$KO)
# 
# annot.intersect_18_20_00 <- merge(shared.all, annot, by="Transcript.ID")
# write.csv(annot.intersect_18_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_00.csv")
# 
# annot.intersect_18_20 <- merge(unique.18.20, annot, by="Transcript.ID")
# unique(annot.intersect_18_20$Transcript.ID)
# write.csv(annot.intersect_18_20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20.csv")
# 
# annot.intersect_18_00 <- merge(unique.18.00, annot, by="Transcript.ID")
# unique(annot.intersect_18_00$Transcript.ID)
# write.csv(annot.intersect_18_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_00.csv")
# 
# annot.intersect_20_00 <- merge(unique.20.00, annot, by="Transcript.ID")
# unique(annot.intersect_20_00$Transcript.ID)
# write.csv(annot.intersect_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_20_00.csv")
# 
# annot.unique.18 <- merge(unique.18, annot, by="Transcript.ID")
# write.csv(annot.unique.18, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_18.csv")
# 
# annot.unique.20 <- merge(unique.20, annot, by="Transcript.ID")
# write.csv(annot.unique.20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_20.csv")
# 
# annot.unique.00 <- merge(unique.00, annot, by="Transcript.ID")
# write.csv(annot.unique.00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_00.csv")
# 
# 
# ##### GO ENRICHMENT ANALYSIS #####
# #Analysis for All DEGs (both Host and Sym together, only 1 sym gene in DEG list)
# #read in data files
# 
# ##### GO enrichment of shared DE genes #####
# ALL<-row.names(data) #set the all transcripts list
# DEG <- as.character(intersections$Transcript.ID) #set the enrichment test list
# LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
# gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
# colnames(gn.keep) <- c("V1") #name columns
# lens <- as.data.frame(LENGTH) #set lengths
# LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
# GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
# GO <- as.data.frame(GO) #filter GO to subset same count filter as above
# GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
# GOs[is.na(GOs)] <- "unknown" #list NA as unknown
# GOs <- merge(GOs, LENGTH, by="V1") #merge
# GOs <- GOs[,1:2] #remove length 
# splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
# GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
# IDs <- row.names(data) #set ID names
# 
# #change contig lists to vectors
# ALL.vector <-c(t(ALL))
# DEG.vector <-c(t(DEG))
# ID.vector <- LENGTH$V1
# LENGTH.vector <-LENGTH$V2
# 
# gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
# names(gene.vector)=ALL.vector #set names
# DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene
# 
# #Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
# 
# #How many enriched GO terms do we have
# class(GO.wall)
# head(GO.wall)
# tail(GO.wall)
# nrow(GO.wall)
# 
# #Find only enriched GO terms that are statistically significant at cutoff 
# enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
# enriched.GO.05<-data.frame(enriched.GO.05.a)
# colnames(enriched.GO.05) <- c("category")
# enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")
# 
# MF.INT <- subset(enriched.GO.05, ontology=="MF")
# MF.INT <- MF.INT[order(-MF.INT$numDEInCat),]
# CC.INT <- subset(enriched.GO.05, ontology=="CC")
# CC.INT <- CC.INT[order(-CC.INT$numDEInCat),]
# BP.INT <- subset(enriched.GO.05, ontology=="BP")
# BP.INT <- BP.INT[order(-BP.INT$numDEInCat),]
# write.csv(MF.INT , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_INT.csv")
# write.csv(CC.INT , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_INT.csv")
# write.csv(BP.INT , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_INT.csv")
# 
# 
# ##### GO enrichment of DE genes Unique to 18:00 #####
# ALL<-row.names(data)
# DEG <- as.character(unique.18$Transcript.ID) #set the enrichment test list
# LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
# gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
# colnames(gn.keep) <- c("V1") #name columns
# lens <- as.data.frame(LENGTH) #set lengths
# LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
# GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
# GO <- as.data.frame(GO) #filter GO to subset same count filter as above
# GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
# GOs[is.na(GOs)] <- "unknown" #list NA as unknown
# GOs <- merge(GOs, LENGTH, by="V1") #merge
# GOs <- GOs[,1:2] #remove length 
# splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
# GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
# IDs <- row.names(data) #set ID names
# 
# #change contig lists to vectors
# ALL.vector <-c(t(ALL))
# DEG.vector <-c(t(DEG))
# ID.vector <- LENGTH$V1
# LENGTH.vector <-LENGTH$V2
# 
# gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
# names(gene.vector)=ALL.vector #set names
# DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene
# 
# #Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
# 
# #How many enriched GO terms do we have
# class(GO.wall)
# head(GO.wall)
# tail(GO.wall)
# nrow(GO.wall)
# 
# #Find only enriched GO terms that are statistically significant at cutoff 
# enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
# enriched.GO.05<-data.frame(enriched.GO.05.a)
# colnames(enriched.GO.05) <- c("category")
# enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")
# 
# MF.18 <- subset(enriched.GO.05, ontology=="MF")
# MF.18 <- MF.18[order(-MF.18$numDEInCat),]
# CC.18 <- subset(enriched.GO.05, ontology=="CC")
# CC.18 <- CC.18[order(-CC.18$numDEInCat),]
# BP.18 <- subset(enriched.GO.05, ontology=="BP")
# BP.18 <- BP.18[order(-BP.18$numDEInCat),]
# write.csv(MF.18 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_18.csv")
# write.csv(CC.18 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_18.csv")
# write.csv(BP.18 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_18.csv")
# 
# 
# 
# ##### GO enrichment of DE genes Unique to 20:00 #####
# ALL<-row.names(data)
# DEG <- as.character(unique.20$Transcript.ID) #row.names(sig) 
# LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
# gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
# colnames(gn.keep) <- c("V1") #name columns
# lens <- as.data.frame(LENGTH) #set lengths
# LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
# GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
# GO <- as.data.frame(GO) #filter GO to subset same count filter as above
# GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
# GOs[is.na(GOs)] <- "unknown" #list NA as unknown
# GOs <- merge(GOs, LENGTH, by="V1") #merge
# GOs <- GOs[,1:2] #remove length 
# splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
# GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
# IDs <- row.names(data) #set ID names
# 
# #change contig lists to vectors
# ALL.vector <-c(t(ALL))
# DEG.vector <-c(t(DEG))
# ID.vector <- LENGTH$V1
# LENGTH.vector <-LENGTH$V2
# 
# gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
# names(gene.vector)=ALL.vector #set names
# DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene
# 
# #Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
# 
# #How many enriched GO terms do we have
# class(GO.wall)
# head(GO.wall)
# tail(GO.wall)
# nrow(GO.wall)
# 
# #Find only enriched GO terms that are statistically significant at cutoff 
# enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
# enriched.GO.05<-data.frame(enriched.GO.05.a)
# colnames(enriched.GO.05) <- c("category")
# enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")
# 
# MF.20 <- subset(enriched.GO.05, ontology=="MF")
# MF.20 <- MF.20[order(-MF.20$numDEInCat),]
# CC.20 <- subset(enriched.GO.05, ontology=="CC")
# CC.20 <- CC.20[order(-CC.20$numDEInCat),]
# BP.20 <- subset(enriched.GO.05, ontology=="BP")
# BP.20 <- BP.20[order(-BP.20$numDEInCat),]
# write.csv(MF.20 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_20.csv")
# write.csv(CC.20 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_20.csv")
# write.csv(BP.20 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_20.csv")
# 
# ##### GO enrichment of DE genes Unique to 00:00 #####
# ALL<-row.names(data)
# DEG <- as.character(unique.00$Transcript.ID) #row.names(sig) 
# LENGTH <-read.table("Trinity.fasta.seq_lens", sep = "", header=F) #reads in a list of gene lengths
# gn.keep <- as.data.frame(gn.keep) #filter length to subset same count filter as above
# colnames(gn.keep) <- c("V1") #name columns
# lens <- as.data.frame(LENGTH) #set lengths
# LENGTH <- merge(as.data.frame(LENGTH), gn.keep, by="V1") #merge lengths and transcripts
# GO <- read.table("go_annotations.txt", sep = "", header=F, stringsAsFactors=FALSE) #id GO annotations
# GO <- as.data.frame(GO) #filter GO to subset same count filter as above
# GOs <- merge(as.data.frame(GO), gn.keep, by="V1", all=T) #merge GO, ID
# GOs[is.na(GOs)] <- "unknown" #list NA as unknown
# GOs <- merge(GOs, LENGTH, by="V1") #merge
# GOs <- GOs[,1:2] #remove length 
# splitted <- strsplit(as.character(GOs$V2.x), ",") #slit into multiple GO ids
# GO.terms <- data.frame(v1 = rep.int(GOs$V1, sapply(splitted, length)), v2 = unlist(splitted)) #list all GOs with each assigned transcript
# IDs <- row.names(data) #set ID names
# 
# #change contig lists to vectors
# ALL.vector <-c(t(ALL))
# DEG.vector <-c(t(DEG))
# ID.vector <- LENGTH$V1
# LENGTH.vector <-LENGTH$V2
# 
# gene.vector=as.integer(ALL.vector%in%DEG.vector) #Construct new vector with 1 for DEG and 0 for others
# names(gene.vector)=ALL.vector #set names
# DEG.pwf<-nullp(gene.vector, ID.vector, bias.data=LENGTH.vector) #weight vector by length of gene
# 
# #Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# GO.wall<-goseq(DEG.pwf, ID.vector, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
# 
# #How many enriched GO terms do we have
# class(GO.wall)
# head(GO.wall)
# tail(GO.wall)
# nrow(GO.wall)
# 
# #Find only enriched GO terms that are statistically significant at cutoff 
# enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
# enriched.GO.05<-data.frame(enriched.GO.05.a)
# colnames(enriched.GO.05) <- c("category")
# enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")
# 
# MF.00 <- subset(enriched.GO.05, ontology=="MF")
# MF.00 <- MF.00[order(-MF.00$numDEInCat),]
# CC.00 <- subset(enriched.GO.05, ontology=="CC")
# CC.00 <- CC.00[order(-CC.00$numDEInCat),]
# BP.00 <- subset(enriched.GO.05, ontology=="BP")
# BP.00 <- BP.00[order(-BP.00$numDEInCat),]
# write.csv(MF.00 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_Sig_Enriched_GO.05_00.csv")
# write.csv(CC.00 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_Sig_Enriched_GO.05_00.csv")
# write.csv(BP.00 , file = "/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_Sig_Enriched_GO.05_00.csv")
# 
# 
# ##### GO SLIM #####
# ##### Slim for Biological Processes #####
# BP.Ids <- as.character(BP.INT$category)
# myCollection <- GOCollection(BP.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# BP.slim.INT <- goSlim(myCollection, slim, "BP")
# BP.slim.INT <-BP.slim.INT[,c(1,3)]
# 
# BP.Ids <- as.character(BP.18$category)
# myCollection <- GOCollection(BP.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# BP.slim.18 <- goSlim(myCollection, slim, "BP")
# BP.slim.18 <-BP.slim.18[,c(1,3)]
# 
# BP.Ids <- as.character(BP.20$category)
# myCollection <- GOCollection(BP.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# BP.slim.20 <- goSlim(myCollection, slim, "BP")
# BP.slim.20 <-BP.slim.20[,c(1,3)]
# 
# BP.Ids <- as.character(BP.00$category)
# myCollection <- GOCollection(BP.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# BP.slim.00 <- goSlim(myCollection, slim, "BP")
# BP.slim.00 <-BP.slim.00[,c(1,3)]
# 
# BP.slim.all <- merge(BP.slim.INT, BP.slim.18, by="Term", all = T )
# BP.slim.all <- merge(BP.slim.all, BP.slim.20, by="Term", all = T )
# colnames(BP.slim.all) <- c("Term", "Intersection", "18:00", "20:00")
# BP.slim.all <- merge(BP.slim.all, BP.slim.00, by="Term", all = T )
# colnames(BP.slim.all) <- c("Term", "Intersection", "18:00", "20:00", "00:00")
# BP.slim.all 
# BP.slim.all  <- BP.slim.all[order(BP.slim.all$Intersection),]
# BP.slim.all$Term <- gsub("biological_process", "unknown", BP.slim.all$Term)
# BP.slim.all
# BP.slim.all <- BP.slim.all[rowSums(BP.slim.all[, -1])>3, ]
# BP.slim.all <-head(BP.slim.all,-1)
# 
# 
# pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/BP_enrichment_SPvsNO.pdf")
# par(mfrow=c(1,4))
# par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
# barplot(BP.slim.all$Intersection, horiz = T, col="blue", names.arg=BP.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,105), main = "ALL")
# barplot(BP.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,105), main = "18:00")
# barplot(BP.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,105), main = "20:00")
# barplot(BP.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,105), main = "00:00")
# dev.off()
# 
# ##### Slim for Molecular Function #####
# MF.Ids <- as.character(MF.INT$category)
# myCollection <- GOCollection(MF.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# MF.slim.INT <- goSlim(myCollection, slim, "MF")
# MF.slim.INT <-MF.slim.INT[,c(1,3)]
# 
# MF.Ids <- as.character(MF.18$category)
# myCollection <- GOCollection(MF.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# MF.slim.18 <- goSlim(myCollection, slim, "MF")
# MF.slim.18 <-MF.slim.18[,c(1,3)]
# 
# MF.Ids <- as.character(MF.20$category)
# myCollection <- GOCollection(MF.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# MF.slim.20 <- goSlim(myCollection, slim, "MF")
# MF.slim.20 <-MF.slim.20[,c(1,3)]
# 
# MF.Ids <- as.character(MF.00$category)
# myCollection <- GOCollection(MF.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# MF.slim.00 <- goSlim(myCollection, slim, "MF")
# MF.slim.00 <-MF.slim.00[,c(1,3)]
# 
# MF.slim.all <- merge(MF.slim.INT, MF.slim.18, by="Term", all = T )
# MF.slim.all <- merge(MF.slim.all, MF.slim.20, by="Term", all = T )
# colnames(MF.slim.all) <- c("Term", "Intersection", "18:00", "20:00")
# MF.slim.all <- merge(MF.slim.all, MF.slim.00, by="Term", all = T )
# colnames(MF.slim.all) <- c("Term", "Intersection", "18:00", "20:00", "00:00")
# MF.slim.all 
# MF.slim.all  <- MF.slim.all [order(MF.slim.all$Intersection),]
# MF.slim.all$Term <- gsub("molecular_function", "unknown", MF.slim.all$Term)
# MF.slim.all
# MF.slim.all <- MF.slim.all[rowSums(MF.slim.all[, -1])>1, ]
# MF.slim.all <-head(MF.slim.all,-1)
# 
# pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/MF_enrichment_SPvsNO.pdf")
# par(mfrow=c(1,4))
# par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
# barplot(MF.slim.all$Intersection, horiz = T, col="blue", names.arg=MF.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,65), main = "ALL")
# barplot(MF.slim.all$`18:00`, horiz = T, col="yellow", xlim=c(0,65), main = "18:00")
# barplot(MF.slim.all$`20:00`, horiz = T, col="red", xlim=c(0,65), main = "20:00")
# barplot(MF.slim.all$`00:00`, horiz = T, col="green", xlim=c(0,65), main = "00:00")
# dev.off()
# 
# ##### Slim for Cellular Component #####
# CC.Ids <- as.character(CC.INT$category)
# myCollection <- GOCollection(CC.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# CC.slim.INT <- goSlim(myCollection, slim, "CC")
# CC.slim.INT <-CC.slim.INT[,c(1,3)]
# 
# CC.Ids <- as.character(CC.18$category)
# myCollection <- GOCollection(CC.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# CC.slim.18 <- goSlim(myCollection, slim, "CC")
# CC.slim.18 <-CC.slim.18[,c(1,3)]
# 
# CC.Ids <- as.character(CC.20$category)
# myCollection <- GOCollection(CC.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# CC.slim.20 <- goSlim(myCollection, slim, "CC")
# CC.slim.20 <-CC.slim.20[,c(1,3)]
# 
# CC.Ids <- as.character(CC.00$category)
# myCollection <- GOCollection(CC.Ids)
# fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
# slim <- getOBOCollection(fl, evidenceCode="TAS")
# CC.slim.00 <- goSlim(myCollection, slim, "CC")
# CC.slim.00 <-CC.slim.00[,c(1,3)]
# 
# CC.slim.all <- merge(CC.slim.INT, CC.slim.18, by="Term", all = T )
# CC.slim.all <- merge(CC.slim.all, CC.slim.20, by="Term", all = T )
# colnames(CC.slim.all) <- c("Term", "Intersection", "18:00", "20:00")
# CC.slim.all <- merge(CC.slim.all, CC.slim.00, by="Term", all = T )
# colnames(CC.slim.all) <- c("Term", "Intersection", "18:00", "20:00", "00:00")
# CC.slim.all 
# CC.slim.all  <- CC.slim.all[order(CC.slim.all$Intersection),]
# CC.slim.all$Term <- gsub("cellular_component", "unknown", CC.slim.all$Term)
# CC.slim.all
# CC.slim.all <- CC.slim.all[rowSums(CC.slim.all[, -1])>1, ]
# CC.slim.all <-head(CC.slim.all,-1)
# 
# pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/CC_enrichment_SPvsNO.pdf")
# par(mfrow=c(1,4))
# par(mar=c(3,1,1,0), oma = c(0.1, 10, 0.1, 0.5))
# barplot(CC.slim.all$Intersection, horiz = T, col="black", names.arg=CC.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,20), main = "ALL")
# barplot(CC.slim.all$`18:00`, horiz = T, col="steelblue1", xlim=c(0,30), main = "18:00")
# barplot(CC.slim.all$`20:00`, horiz = T, col="yellow3", xlim=c(0,30), main = "20:00")
# barplot(CC.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,30), main = "00:00")
# dev.off()
# 
# ##### Visualizing GO Enrichment #####
# #Plot all Enriched GOs
# pdf(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/GO_enrichment_SPvsNO.pdf")
# par(mfrow=c(2,4))
# par(mar=c(3,1,1,0), oma = c(0.1, 15, 0.1, 0.5))
# barplot(BP.slim.all$Intersection, horiz = T, col="black", names.arg=BP.slim.all$Term, las=1, cex.names = 0.4, xlim=c(0,25), main = "ALL")
# text(labels = 'Biological Process', xpd = NA, srt = 90, x=-26, y=15, cex=1)
# barplot(BP.slim.all$`18:00`, horiz = T, col="steelblue1", xlim=c(0,25), main = "18:00")
# barplot(BP.slim.all$`20:00`, horiz = T, col="yellow3", xlim=c(0,25), main = "20:00")
# barplot(BP.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,25), main = "00:00")
# 
# barplot(MF.slim.all$Intersection, horiz = T, col="black", names.arg=MF.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,25), main = "ALL")
# text(labels = 'Molecular Function', xpd = NA, srt = 90, x=-26, y=8, cex=1)
# barplot(MF.slim.all$`18:00`, horiz = T, col="steelblue1", xlim=c(0,25), main = "18:00")
# barplot(MF.slim.all$`20:00`, horiz = T, col="yellow3", xlim=c(0,25), main = "20:00")
# barplot(MF.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,25), main = "00:00")
# 
# # barplot(CC.slim.all$Intersection, horiz = T, col="black", names.arg=CC.slim.all$Term, las=1, cex.names = 0.6, xlim=c(0,25), main = "ALL")
# # text(labels = 'Cellular Component', xpd = NA, srt = 90, x=-26, y=8, cex=1)
# # barplot(CC.slim.all$`18:00`, horiz = T, col="steelblue1", xlim=c(0,25), main = "18:00")
# # barplot(CC.slim.all$`20:00`, horiz = T, col="yellow3", xlim=c(0,25), main = "20:00")
# # barplot(CC.slim.all$`00:00`, horiz = T, col="lightpink2", xlim=c(0,25), main = "00:00")
# dev.off()
# 
# 
# ##### Pairwise Contrasts of Non_Spawning Corals Between Times #####
# resultsNames(DEG.int)
# #"Intercept"     "groupNO0.00"   "groupNO18.00"  "groupNO20.00"  "groupYES0.00"  "groupYES18.00" "groupYES20.00"
# 
# # Time point Yes No contrasts
# #NO = DEG between different time points for NONspawning corals 
# DEG.NO_18_20 <- results(DEG.int, contrast=c("group", "NO18:00", "NO20:00")) #contrasts between SP NSP at each time
# DEG.NO_18_20 #view results
# 
# sig.num.NO_18_20 <- sum(DEG.NO_18_20$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.NO_18_20 <- subset(DEG.NO_18_20, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.NO_18_20 <- data[which(rownames(data) %in% rownames(sig.NO_18_20)),] #subset list of sig transcripts from original count data
# sig.list.NO_18_20
# 
# DEG.NO_18_00 <- results(DEG.int, contrast=c("group", "NO18:00", "NO0:00")) #contrasts between SP NSP at each time
# DEG.NO_18_00 #view results
# 
# sig.num.NO_18_00 <- sum(DEG.NO_18_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.NO_18_00 <- subset(DEG.NO_18_00, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.NO_18_00 <- data[which(rownames(data) %in% rownames(sig.NO_18_00)),] #subset list of sig transcripts from original count data
# sig.list.NO_18_00
# 
# DEG.NO_20_00 <- results(DEG.int, contrast=c("group", "NO20:00", "NO0:00")) #contrasts between SP NSP at each time
# DEG.NO_20_00 #view results
# 
# sig.num.NO_20_00 <- sum(DEG.NO_20_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.NO_20_00 <- subset(DEG.NO_20_00, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.NO_20_00 <- data[which(rownames(data) %in% rownames(sig.NO_20_00)),] #subset list of sig transcripts from original count data
# sig.list.NO_20_00
# 
# times_18_20 <- merge(intersections, norm.sig.counts, by="Transcript.ID")
# 
# ##### Visualize Intersections of Differental Expression of Non-Spawning Corals at each time #####
# NO_18_20 <- data.frame(rownames(sig.list.NO_18_20)) #identify DE transcript list from NO betwee 18:00 and 20
# NO_18_20$NO_18_20.Count <- 1 #assign all transcripts count of 1
# colnames(NO_18_20) <- c("Transcript.ID", "NO_18_20.Count") #name columns
# 
# NO_18_00 <- data.frame(rownames(sig.list.NO_18_00)) #identify DE transcript list from NO betwee 18:00 and 00
# NO_18_00$NO_18_00.Count <- 1 #assign all transcripts count of 1
# colnames(NO_18_00) <- c("Transcript.ID", "NO_18_00.Count") #name columns
# 
# NO_20_00 <- data.frame(rownames(sig.list.NO_20_00)) #identify DE transcript list from NO betwee 20:00 and 00
# NO_20_00$NO_20_00.Count <- 1 #assign all transcripts count of 1
# colnames(NO_20_00) <- c("Transcript.ID", "NO_20_00.Count") #name columns
# 
# L_18_20 <- merge(NO_18_20, NO_18_00, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
# L_18_20_00 <- merge(L_18_20, NO_20_00, by="Transcript.ID", all=T ) #merge lists sequentially with no removal
# L_18_20_00$NO_18_20.Count[is.na(L_18_20_00$NO_18_20.Count)] <- 0 #assign NA=0 count i.e., not DE
# L_18_20_00$NO_18_00.Count[is.na(L_18_20_00$NO_18_00.Count)] <- 0 #assign NA=0 count i.e., not DE
# L_18_20_00$NO_20_00.Count[is.na(L_18_20_00$NO_20_00.Count)] <- 0 #assign NA=0 count i.e., not DE
# 
# # Plot Venn Diagram of shared DE between time points
# vennData.NO<-L_18_20_00[,2:4] #set venn data to DE counts only
# a <- vennCounts(vennData.NO) # compute classification counts
# colnames(a) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00', "Counts") #set catagory names
# jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.No.Spawn.jpg")
# vennDiagram(a, main='DEG between Time Points for Nonspawning Corals') #draw venn diagram
# dev.off()
# 
# colnames(vennData.NO) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00') #set catagory names
# jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/NOSP_Timepoints_DEG_Intersections.jpg")
# upset(vennData.NO, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
# dev.off()
# 
# #Look for match between all lists
# NO.shared.all <- Reduce(intersect, list(rownames(sig.list.NO_18_20),rownames(sig.list.NO_18_00),rownames(sig.list.NO_20_00))) #idenitfy shared DE transcripts between each time point 
# NO.shared.all #view results
# NO.shared.all <- as.data.frame(NO.shared.all)
# colnames(NO.shared.all)[1] <- "Transcript.ID"
# write.csv(NO.shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_shared_NOSpawn_18_20_00.csv")
# 
# #DEG between spawning and not spawning corals Found in all time points
# intersections.NO <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==1 & NO_20_00.Count==1) #list DEG found in all time points
# 
# #DEG between not spawning corals Unique to each time point
# unique.NO_18_20 <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==0 & NO_20_00.Count==0) #list genes unique to 18:00
# unique.NO_18_00 <- subset(L_18_20_00, NO_18_20.Count==0 & NO_18_00.Count==1 & NO_20_00.Count==0) #list genes unique to 20:00
# unique.NO_20_00 <- subset(L_18_20_00, NO_18_20.Count==0 & NO_18_00.Count==0 & NO_20_00.Count==1) #list genes unique to 00:00
# 
# #DEG between spawning and not spawning corals Shared unique between pairs
# shared.NO_18_20 <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==1 & NO_20_00.Count==0) #list genes common to 18:00 and 20:00
# shared.NO_18_00 <- subset(L_18_20_00, NO_18_20.Count==1 & NO_18_00.Count==0 & NO_20_00.Count==1) #list genes common to 18:00 and 00:00
# shared.NO_20_00 <- subset(L_18_20_00, NO_18_20.Count==0 & NO_18_00.Count==1 & NO_20_00.Count==1) #list genes common to 20:00 and 00:00
# 
# #Combine with DE list with annotation data 
# 
# annot.intersect_18_20_00 <- merge(NO.shared.all, annot, by="Transcript.ID")
# write.csv(annot.intersect_18_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_00_NO.csv")
# 
# annot.intersect_18_20 <- merge(unique.NO_18_20, annot, by="Transcript.ID")
# unique(annot.intersect_18_20$Transcript.ID)
# write.csv(annot.intersect_18_20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_18_20_NO.csv")
# 
# annot.intersect_18_00 <- merge(unique.NO_18_00, annot, by="Transcript.ID")
# unique(annot.intersect_18_00$Transcript.ID)
# write.csv(annot.intersect_18_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_18_00_NO.csv")
# 
# annot.intersect_20_00 <- merge(unique.NO_20_00, annot, by="Transcript.ID")
# unique(annot.intersect_20_00$Transcript.ID)
# write.csv(annot.intersect_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_20_00_NO.csv")
# 
# ##### Pairwise Contrasts of Spawning Corals Between Times #####
# # YES = DEG between different time points for spawning corals
# DEG.YES_18_20 <- results(DEG.int, contrast=c("group", "YES18:00", "YES20:00")) #contrasts between SP NSP at each time
# DEG.YES_18_20 #view results
# 
# sig.num.YES_18_20 <- sum(DEG.YES_18_20$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.YES_18_20 <- subset(DEG.YES_18_20, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.YES_18_20 <- data[which(rownames(data) %in% rownames(sig.YES_18_20)),] #subset list of sig transcripts from original count data
# sig.list.YES_18_20
# 
# DEG.YES_18_00 <- results(DEG.int, contrast=c("group", "YES18:00", "YES0:00")) #contrasts between SP NSP at each time
# DEG.YES_18_00 #view results
# 
# sig.num.YES_18_00 <- sum(DEG.YES_18_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.YES_18_00 <- subset(DEG.YES_18_00, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.YES_18_00 <- data[which(rownames(data) %in% rownames(sig.YES_18_00)),] #subset list of sig transcripts from original count data
# sig.list.YES_18_00
# 
# DEG.YES_20_00 <- results(DEG.int, contrast=c("group", "YES20:00", "YES0:00")) #contrasts between SP NSP at each time
# DEG.YES_20_00 #view results
# 
# sig.num.YES_20_00 <- sum(DEG.YES_20_00$padj <0.05, na.rm=T) #the number of significant p values with 1%FDR
# sig.YES_20_00 <- subset(DEG.YES_20_00, padj<0.05,) #identify signficant pvalues with 1%FDR
# sig.list.YES_20_00 <- data[which(rownames(data) %in% rownames(sig.YES_20_00)),] #subset list of sig transcripts from original count data
# sig.list.YES_20_00
# 
# 
# ##### Visualize Intersections of Differental Expression of Spawning Corals at each time #####
# YES_18_20 <- data.frame(rownames(sig.list.YES_18_20)) #identify DE transcript list from YES betwee 18:00 and 20
# YES_18_20$YES_18_20.Count <- 1 #assign all transcripts count of 1
# colnames(YES_18_20) <- c("Transcript.ID", "YES_18_20.Count") #name columns
# 
# YES_18_00 <- data.frame(rownames(sig.list.YES_18_00)) #identify DE transcript list from YES betwee 18:00 and 00
# YES_18_00$YES_18_00.Count <- 1 #assign all transcripts count of 1
# colnames(YES_18_00) <- c("Transcript.ID", "YES_18_00.Count") #name columns
# 
# YES_20_00 <- data.frame(rownames(sig.list.YES_20_00)) #identify DE transcript list from YES betwee 20:00 and 00
# YES_20_00$YES_20_00.Count <- 1 #assign all transcripts count of 1
# colnames(YES_20_00) <- c("Transcript.ID", "YES_20_00.Count") #name columns
# 
# Y_18_20 <- merge(YES_18_20, YES_18_00, by="Transcript.ID", all=T ) #merge lists sequentially with YES removal
# Y_18_20_00 <- merge(Y_18_20, YES_20_00, by="Transcript.ID", all=T ) #merge lists sequentially with YES removal
# Y_18_20_00$YES_18_20.Count[is.na(Y_18_20_00$YES_18_20.Count)] <- 0 #assign NA=0 count i.e., YESt DE
# Y_18_20_00$YES_18_00.Count[is.na(Y_18_20_00$YES_18_00.Count)] <- 0 #assign NA=0 count i.e., YESt DE
# Y_18_20_00$YES_20_00.Count[is.na(Y_18_20_00$YES_20_00.Count)] <- 0 #assign NA=0 count i.e., YESt DE
# 
# # Plot Venn Diagram of shared DE between time points
# vennData.YES<-Y_18_20_00[,2:4] #set venn data to DE counts only
# a <- vennCounts(vennData.YES) # compute classification counts
# #a<- a[rowSums(a[, -1])>0, ]
# colnames(a) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00', "Counts") #set catagory names
# #a<- a[,-4]
# jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Venn.DEG.YES.Spawn.jpg")
# vennDiagram(a, main='DEG between Time Points for YES spawning Corals') #draw venn diagram
# dev.off()
# 
# colnames(vennData.YES) <- c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00') #set catagory names
# jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/YESSP_Timepoints_DEG_Intersections.jpg")
# upset(vennData.YES, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
# dev.off()
# 
# 
# # jpeg(file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Timepoints_DEG_Intersections.jpg")
# # layout(mat=matrix(c(1,2),nrow=1,ncol=2,byrow=FALSE))
# # upset(vennData.NO, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
# # upset(vennData.YES, sets = c('18:00 vs 20:00', '18:00 vs 00:00', '20:00 vs 00:00'), order.by = "degree")
# # dev.off()
# 
# #Look for match between all lists
# YES.shared.all <- Reduce(intersect, list(rownames(sig.list.YES_18_20),rownames(sig.list.YES_18_00),rownames(sig.list.YES_20_00))) #idenitfy shared DE transcripts between each time point 
# YES.shared.all #view results
# YES.shared.all <- as.data.frame(YES.shared.all)
# colnames(YES.shared.all)[1] <- "Transcript.ID"
# write.csv(YES.shared.all, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/DEG_All_DEG_18_20_00_Yes.csv")
# 
# #DEG between spawning and YESt spawning corals Found in all time points
# intersections.YES <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==1 & YES_20_00.Count==1) #list DEG found in all time points
# 
# #DEG between YESt spawning corals Unique to each time point
# unique.YES_18_20 <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==0 & YES_20_00.Count==0) #list genes unique to 18:00
# unique.YES_18_00 <- subset(Y_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==1 & YES_20_00.Count==0) #list genes unique to 20:00
# unique.YES_20_00 <- subset(Y_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==0 & YES_20_00.Count==1) #list genes unique to 00:00
# 
# #DEG between spawning and YESt spawning corals Shared unique between pairs
# shared.YES_18_20 <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==1 & YES_20_00.Count==0) #list genes common to 18:00 and 20:00
# shared.YES_18_00 <- subset(Y_18_20_00, YES_18_20.Count==1 & YES_18_00.Count==0 & YES_20_00.Count==1) #list genes common to 18:00 and 00:00
# shared.YES_20_00 <- subset(Y_18_20_00, YES_18_20.Count==0 & YES_18_00.Count==1 & YES_20_00.Count==1) #list genes common to 20:00 and 00:00
# 
# #Combine with DE list with annotation data 
# 
# annot.intersect_18_20_00 <- merge(YES.shared.all, annot, by="Transcript.ID")
# write.csv(annot.intersect_18_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_shared_DEG_18_20_00_YES.csv")
# 
# annot.intersect_18_20 <- merge(unique.YES_18_20, annot, by="Transcript.ID")
# unique(annot.intersect_18_20$Transcript.ID)
# write.csv(annot.intersect_18_20, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_18_20_YES.csv")
# 
# annot.intersect_18_00 <- merge(unique.YES_18_00, annot, by="Transcript.ID")
# unique(annot.intersect_18_00$Transcript.ID)
# write.csv(annot.intersect_18_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_18_00_YES.csv")
# 
# annot.intersect_20_00 <- merge(unique.YES_20_00, annot, by="Transcript.ID")
# unique(annot.intersect_20_00$Transcript.ID)
# write.csv(annot.intersect_20_00, file="/Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Output/Annotated_unique_DEG_20_00_YES.csv")
# 
# #####

##################################################################
#                        Assignment of FGT                      #
#           Author:B224911    April 2023   version 5            #
##################################################################

## This script is the pipeline of gene expression differential analysis between 
## microarray datasets.
## All codes of this script (except for the volcanoplot and heatmap part)
## are based  on the turtorials in the course 
## Functional Genomic Technologies of Edinburgh University.

## step1 Load the package and data 
##load the required libraries for the tutorials...
library(affy)
library(limma)
library(mouse4302.db)
library(annotate)
library(ggplot2)
library(ggrepel)
library(pheatmap) 
(.packages()) # check out the current package

# Load the target file into an AnnotatedDataFrame object
Target <- read.AnnotatedDataFrame("ICA_targets.txt",header=TRUE,row.names=1,as.is=TRUE)
target <-read.table ("ICA_targets.txt", header= T,as.is = T)
View(target)

# quickly load all CEL files in the R working directory
rawdata <- ReadAffy()
# View a summary of the example data
rawdata



##Step2 Quality control of raw data
# Quality control plots
png("qc_hist_1.png")
hist(rawdata,
     # labels= TRUE,
     main="QC Histogram of microarray data")
dev.off()
# View(rawdata)

# And a boxplot with different colour per sample group
colours <- c(rep("green",3),rep("purple",3))
png("qc_boxplot_1.png")
boxplot(rawdata, 
        col=colours, 
        names=target$Name, 
        main="QC boxplot of microarray data",
        xlab="samples", 
        ylab="intensity",
        las=2)
dev.off()


## Step3 Data normalisation and Quality control after normalisation
# Normalise the data using RMA
RMA_data <- rma(rawdata)
# View(RMA_data)

# To obtain a matrix of the expression values, use exprs()
exp_data <- exprs(RMA_data)
View(exp_data)

# Boxplot to observe the results of normalisation
# Notice differences with the boxplot from the raw data
png("Nqc_boxplot1.png")
boxplot(exp_data, 
        col=colours,
        names=target$Name, 
        main="Normalized-QC boxplot of microarray data",
        xlab="samples", 
        ylab="intensity",
        las=2)
dev.off()

# MA plot of the normalised samples 
png("N_MVA_plot_normal.png", width=1000, height=600) 
par(cex=0.3) 
mva.pairs(exp_data[, c(1, 2, 3, 4, 5, 6)])
dev.off()

# The same plot for the non-normalised raw data
# Note that the mva.pairs call below only plots a few of the #samples – you may wish to plot them all but this is slow
png("MVA_plot_perfect.png", width=1000, height=600)
mva.pairs(pm(rawdata)[, c(1, 2, 3, 4, 5, 6)])
dev.off()



## Step4 Hierarchical clustering of normalised data and Principal Components Analysis (PCA) of normalised data
# To facilitate interpretation, replace the columns, header,currently
# displaying the filename, to show the name of each sample
pData(Target)
View(exp_data)
colnames(exp_data) <- rownames(pData(Target))

# Pearson’s Correlation Coefficient
cluster <-hclust(as.dist(1-cor(exp_data, method="pearson")), method="average")
png("Hier_cluster_plot1.png")
plot(cluster)
dev.off()

#load the library
library(scatterplot3d)
# Perform PCA
pca <- prcomp(t(exp_data), scale=T)
# Plot the PCA results
png("PCA_3d_plot1.png")
s3d<-scatterplot3d(pca$x[,1:3], type="h", pch=19, color=colours,main="PCA 3D Plot")
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(exp_data),pos = 3,offset = 0.5)
dev.off()



## Step5 Fold Filtering and Check Sample Name Order and Build a Fold Change Vector.
## obtaining a matrix of expression values
exprsvals <- exprs(RMA_data)
# RMA outputs log2 data while MAS5 outputs linear data
# To convert from log…
exprsvals10 <-2^exprsvals
# check conversion
exprsvals[1:100, ] 
# converted
exprsvals10[1:100, ]
View(exprsvals10)
# check order of sample names
mysamples <- sampleNames(RMA_data)
# display the list
mysamples
# it is useful to obtain a vector of ProbeIDs here
probesets <- probeNames(rawdata)
# display the first 100 ProbeSets
probesets[1:100]

# Calculate the means
# Note mean of the log is not the same as the log of the mean!!
Plac.mean <- apply(exprsvals10[,c("GSM3946100_01_Plac.CEL","GSM3946101_02_Plac.CEL","GSM3946102_03_Plac.CEL")],1,mean)
E2.mean <- apply(exprsvals10[,c("GSM3946103_04_E2.CEL","GSM3946104_05_E2.CEL","GSM3946105_06_E2.CEL")],1,mean)
View(Plac.mean)
View(E2.mean)
# calculate some fold changes
# ES_Plac <- E2.mean/Plac.mean 
Plac_E2 <- Plac.mean/E2.mean
Plac_E2_log2 <- abs(log(Plac.mean/E2.mean))

# build a summary table to hold all the data
# all.data= cbind(E2.mean, Plac.mean, ES_Plac)
all.foldchange.data=cbind(Plac.mean, E2.mean, Plac_E2, Plac_E2_log2)
View(all.foldchange.data)
rownames(all.foldchange.data)
# check the column names
colnames(all.foldchange.data)
nrow(all.foldchange.data)


## step6
## Tidy-Up &Rename Samples
## Annotate the Results with Gene Names
# Check original sample order
sampleNames(RMA_data)
# Rename the samples
sampleNames(RMA_data) <-c("Plac_1","Plac_2","Plac_3","E2_04","E2_05","E2_06")
# Check the samples have renamed
sampleNames(RMA_data)

#establish annotation for MOE430v2
#which annotation do we need
# modified from #http://gettinggeneticsdone.blogspot.co.uk/2012/01/annotating-limma-#results-with-gene.html
RMA_data@annotation
library(mouse4302.db)# load chip-specific annotation
# packages in the annotation package
ls("package:mouse4302.db")

# build an annotation table
ID <- featureNames(RMA_data)
ID
Symbol <- getSYMBOL(ID, "mouse4302.db")
View(Symbol)
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
Name
template <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
template[template =="NA"] <- NA #fix padding with NA characters

# assign as feature data of the current Eset
fData(RMA_data) <- template
fData(RMA_data) 



## step7 Statistical Analysis Using Limma
# Build the design matrix
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
design
colnames(design) <- c("Plac", "E2")
# Check it makes sense
sampleNames(RMA_data)
# output the design matrix
design

# Build the Contrasts Matrix
# This instructs Limma which comparisons to make
contrastmatrix <- makeContrasts(Plac-E2, levels=design)
contrastmatrix

# Fit the Linear Model
# issue these commands to fit the model
# and make the contrasts
fit1 <- lmFit(RMA_data, design)
fit2 <- contrasts.fit(fit1, contrastmatrix)
# this last part essentially moderates the t-statistic using
# the borrowed variance approach described in class
fit_eBays <- eBayes(fit2)
fit_treat <- treat(fit2)



## Step8 Writing a Summary Table
## make topTable analysis
topTable(fit_eBays,coef=1,adjust="fdr")
stat_top_results<- topTable(fit_eBays,coef=1, adjust="fdr", number=nrow(RMA_data))
nrow(stat_top_results)
View(stat_top_results)
write.table(stat_top_results,"Stat_top_results.txt")

# make topTreat analysis
topTreat(fit_treat,coef=1,adjust="fdr")
stat_treat_results<- topTreat(fit_treat,coef=1, adjust="fdr", number=nrow(RMA_data))
nrow(stat_treat_results)
View(stat_treat_results)
write.table(stat_treat_results, "Stat_treat_results.txt")

# make the VennDiagram
clas <- classifyTestsF(fit_eBays)
nrow(clas)
View(clas)
png("VennDiagram of differential sample gene.png")
vennDiagram(clas, main = "VennDiagram of differential sample gene")
dev.off

# Base on the stat_top_results,add the other items on foldchangefile
fold_select <- stat_top_results[rownames(all.foldchange.data),]
View(fold_select)
fold_item <- fold_select[, c("ID","Symbol","Name")]
View(fold_item)
foldchange_results <- cbind(fold_item,  all.foldchange.data)
View(foldchange_results)

# write the table of means as an output
write.table(foldchange_results,file="group_means_foldchange.txt", quote=F,
            sep = "\t", col.names = NA)



##step9
## Annotating the Expression Data with Entrez (and other useful) IDs
## Annotate the Expression Data
## Writing a Summary Table
#load the signatures library
Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds")
# View(Mm.H)

# Show the full contents of the annotation package
ls("package:mouse4302.db")
#Show the annotation keys in this database
keytypes(mouse4302.db)

# Mapping from
# http://genomicsclass.github.io/book/pages/mapping_features.ht#ml
sampleNames(RMA_data)

# Here we select from the annotation a number of keys with the primary key being PROBEID
res <- select(mouse4302.db, keys = rownames(RMA_data), columns = c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
# View the top of the table

# find the index of each row of the expression set in the #annotation object res
index <- match(rownames(RMA_data), res$PROBEID)
length(index)

# Use the index to set the phenotypic data in the ExpressionSet
fData(RMA_data) <- res[index, ]
head(fData(RMA_data), 10)
# View(fData(RMA_data))

# Find all rows that don’t have an EntrezID and remove then
No_NA_RMA_data<-RMA_data[is.na(fData(RMA_data)$ENTREZID)==0,]
#View(No_NA_RMA_data)
# No_NA_RMA_data



## step10
## Statistical Enrichment Analysis Using Limma Design & Contrast Matrices
## Run the Enrichment Analysis Using MRoast or Camera or Romer
#convert to indexes
H.indices <- ids2indices(Mm.H,fData(No_NA_RMA_data)$ENTREZID)
# View(H.indices)

# Pick the most suitable enrichment analysis tool to find #enrichment signatures in the data and run this tool So:
# irun mroast
mroast_results1 <-mroast(No_NA_RMA_data,
                         index=H.indices,
                         design=design,
                         contrast=contrastmatrix[,1],
                         adjust.method = "BH")


# run camera
camera_results2 <- camera(No_NA_RMA_data,
                          index=H.indices,
                          design=design,
                          contrast=contrastmatrix[,1])

#  run romer (takes 5mins or so)
romer_results3 <-romer(No_NA_RMA_data,
                       index=H.indices,
                       design=design,
                       contrast=contrastmatrix[,1] )

# View the results
View(mroast_results1)
View(camera_results2)
View(romer_results3)

# write the results to the table 
write.table(mroast_results1,"Mroast_enrichment_1.txt",sep="\t")
write.table(camera_results2,"Camera_enrichment_2.txt",sep="\t")
write.table(romer_results3,"Romer_enrichment_3.txt",sep="\t")
#You can then examine the results in “enrichment.txt”. It is a text file. It can be downloaded to view in a spreadsheet such as Excel.



## step11 draw the volcano plot and heatmap
# try to plot volcano plot with ggplot
# the codeing reference is from:
# https://cloud.tencent.com/developer/article/1486128

head(stat_top_results)
View(stat_top_results)
N_stat_top_results= na.omit(stat_top_results)
N_stat_treat_results= na.omit(stat_treat_results)
View(N_stat_top_results)
View(N_stat_treat_results)

## using |logFC| > 1, p-value < 0.05
foldChange = 1
p.value = 0.05
## generate the data frame
logFC <- N_stat_top_results$logFC
logFC

deg.p.value <- N_stat_top_results$P.Value
deg.p.value

Symbol <- N_stat_top_results$Symbol
Symbol

# build up the dataframe for volcano
vol_data <- data.frame(Symbol = Symbol, logFC = logFC, p.value = deg.p.value)
View(vol_data)

vol_data$group[(vol_data$p.value > 0.05 | vol_data$p.value == "NA") |
                 (abs(vol_data$logFC)< foldChange)]  <- "Not"
vol_data$group[(vol_data$p.value <= 0.05 & vol_data$logFC >= 1)] <-  "Up"
vol_data$group[(vol_data$p.value <= 0.05 & vol_data$logFC <= -1)] <- "Down"

View(vol_data)
sum(vol_data$group == "Up")
sum(vol_data$group == "Down")
sum(vol_data$group == "Not")

## draw the volcano plot
## the volcano plot codeing reference is from:
## https://www.jianshu.com/p/884415e37313
png("Volcano_plot_of_Genes.png")
ggplot(
  # data, mapping,color
  vol_data, aes(x = logFC, y =-log10(p.value)))+
  ggtitle("Volcano plot of Genes")+
  geom_point(aes(color=group), size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  #annotation
  geom_label_repel(
    data = subset(vol_data, p.value < 0.05 & abs(vol_data$logFC) >=1),
    aes(label = Symbol),
    size = 2, fill = "darkred", color = "white",
    box.padding = unit(0.2, "lines"),
    point.padding = unit(0.2, "lines"),
    max.overlaps = 30)+
  # Auxiliary lines
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(p.value),lty=4,col="black",lwd=0.5) +
  #axes 
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  # legends
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
dev.off()


## heatmap plotting
## screening out all differential gene results
All_top_diff <-N_stat_top_results[(N_stat_top_results$P.Value <= p.value & abs(N_stat_top_results$logFC)>=foldChange), ]
All_top_0.5_diff <-N_stat_top_results[(N_stat_top_results$P.Value <= p.value & abs(N_stat_top_results$logFC)>=0.5), ]
All_treat_diff <-N_stat_treat_results[(N_stat_treat_results$P.Value <= p.value & abs(N_stat_treat_results$logFC)>=foldChange), ]

View(All_top_diff)
View(All_treat_diff)
write.csv(All_top_diff, "all_top_diff.csv")
write.csv(All_treat_diff, "all_treat_diff.csv")

# get the heatmap value
heatmap_value <- exp_data[All_top_diff$ID, ]
heatmap_value_0.5 <- exp_data[All_top_0.5_diff$ID, ]

# pheatmap for logfoldchange >=1
png("heatmap_fc_1.png")
pheatmap(heatmap_value,
         scale='column',  
         legend = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "heatmap for differential genes(fold change>=1;p-value<=0.05)",
         fontsize_row = 4,
         fontsize_col = 12)
dev.off()


## pheatmap for logfoldchange >=0.5
png("heatmap_fc_0.5.png")
pheatmap(heatmap_value_0.5,
         scale='column', 
         legend = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "heatmap for differential genes(fold change>=0.5;p-value<=0.05)",
         fontsize_row = 4,
         fontsize_col = 12)
dev.off() 
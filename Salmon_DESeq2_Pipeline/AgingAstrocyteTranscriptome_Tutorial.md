RNA-sequencing Analysis Tutorial
================
Jillybeth Burgado
Last Updated: August 2022

# Introduction

The purpose of this tutorial is to go over how to import and summarize
transcript-level abundance estimates and complete differential
expression (DE) analysis using DESeq2. I also include some basic ideas
for downstream analysis and visualizing your data. The original article
describing these methods can be found below:

> Charlotte Soneson, Michael I. Love, Mark D. Robinson (2015):
> Differential analyses for RNA-seq: transcript-level estimates improve
> gene-level inferences. F1000Research
> <http://dx.doi.org/10.12688/f1000research.7563.1>

**Notes**

-   If you are primarily interested in how to use DESeq2 or
    visualizations of RNA-seq results, skip to Part 2: Differential
    Expression Analysis using DESeq2.
-   Another Bioconductor package, tximeta (Love et al. 2020), extends
    tximport, offering the same functionality, plus the additional
    benefit of automatic addition of annotation metadata for commonly
    used transcriptomes (GENCODE, Ensembl, RefSeq for human and mouse).
    tximeta also offers easy conversion to data objects used by edgeR
    and limma with the makeDGEList function.

# Installing and Loading packages

The first step is to install and load required packages. Visit
<https://www.datacamp.com/community/tutorials/r-packages-guide> for more
information about packages and package repositories. You only need to
install packages once - however, you need to keep track of versions and
you may need to update packages depending on new releases. Below is an
example of how to install a package from the Bioconductor repo.

``` r
#install packages 
if (!requireNamespace("BiocManager", quietly = FALSE))
  install.packages("BiocManager")
BiocManager::install("tximport")

install.packages("pheatmap")
```

``` r
#load necessary packages
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(ggplot2)
library(ggpubr)
library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(tximeta)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(EnhancedVolcano)
library(biomaRt)
library(enrichplot)
library(clusterProfiler)
library(scales)
library(readr)
library(ggnewscale)
library(ggridges)
library(M3C)
library(RColorBrewer)
library(viridis)
library(pathview)
library(rWikiPathways)
library(RCy3)
library(UpSetR)
library(ReactomePA)
library(ggupset)
library(knitr)
library(kableExtra)
```

# Transcript quantification for gene-level analysis

You will need a metadata file that includes your sample names as well as
any group/explanatory variables. Your sample names need to exactly match
your file names (from Salmon), 1 sample per row. I recommend exporting
the names.txt file you created when doing the Salmon analysis on the
server.

``` r
metadata <- read.csv("MattAgingAstrocyte_Metadata2.csv", header = TRUE)
metadata$AgeRegion <- as.character(paste(metadata$Age, metadata$RealRegion, sep= "_"))
```

## Import and Summarize transcript quantification files

The tximeta package has a single function for importing transcript-level
estimates. The type argument is used to specify what software was used
for estimation. A simple list with matrices, “abundance”, “counts”, and
“length”, is returned, where the transcript level information is
summarized to the gene-level. Typically, abundance is provided by the
quantification tools as TPM (transcripts-per-million), while the counts
are estimated counts (possibly fractional), and the “length” matrix
contains the effective gene lengths.

``` r
#First, create a list, called files, with the file path to the quant.sf (Salmon output) files for each sample
files <- file.path("AgingAstrocyteTranscriptome_mm39out", metadata$FullSampleName, "quant.sf")

#Check to make sure all the files are actually located in the right place
all(file.exists(files)) #This should output "TRUE"
```

    ## [1] TRUE

``` r
#create a dataframe containing those filepaths, the names of the samples, and any other metadata you would like to include. 
#I include all of the metadata (from the metadata file)
coldata <- data.frame(files, names = metadata$FullSampleName, metadata)
se <- tximeta(coldata, type="salmon")
```

    ## importing quantifications

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 
    ## found matching transcriptome:
    ## [ Ensembl - Mus musculus - release 105 ]
    ## loading existing EnsDb created: 2022-03-15 17:25:23
    ## loading existing transcript ranges created: 2022-03-15 17:25:27

``` r
#Code examples to view different outputs from the summarized experiment (se)
#assayNames(se)
#rowRanges(se)
#seqinfo(se)
#se.exons <- addExons(se)
#transcript <- as.data.frame(assay(se))

#Summarize to genes
gse <- summarizeToGene(se)
```

    ## loading existing EnsDb created: 2022-03-15 17:25:23
    ## obtaining transcript-to-gene mapping from database
    ## loading existing gene ranges created: 2022-03-15 17:35:44
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
#genes <- as.data.frame(assay(gse))
#You can add other metadata to the gene summarized experiment object (gse) using the following code
columns(org.Mm.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
    ##  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    ## [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MGI"         
    ## [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UNIPROT"

``` r
#gse <- addIds(gse, "SYMBOL", gene=TRUE)
#gse <- addIds(gse, "ENTREZID", gene=TRUE)
#gse <- addIds(gse, "ENSEMBL", gene=TRUE)
```

## Create tpm table and add gene annotations

``` r
annotation <- as.data.frame(mcols(gse))
```

``` r
annotation <- as.data.frame(mcols(gse))

mart <- biomaRt::useDataset("mmusculus_gene_ensembl", useMart("ensembl")) #for default mouse genome

gene_list <- getBM(filters= "ensembl_gene_id",
                   attributes=c("description","ensembl_gene_id"),
                   values= annotation$ENSEMBL,
                   mart= mart,
                   useCache = FALSE)

annotation <-merge(annotation, gene_list, by.x = c('ENSEMBL'), by.y = c('ensembl_gene_id'))

tpm <- as.data.frame(assay(gse, "abundance"))

#To change the column names based on sample names in metadata
#colnames(tpm) <- metadata$PrepID

tpm <- merge(tpm, annotation, by.x = 0, by.y = "gene_id")
tpm <- select(tpm, -1)
rownames(tpm) <- tpm$ENSEMBL

tpm_geneOnly <- filter(tpm, !ENTREZID == "NA")
```

# Differential Expression Analysis using DESeq2

Check out [this
vignette](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
for detailed information regarding the DESeq2 package. Check out [this
other
vignette](https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html)
for information about multi-group comparisons (e.g Sex + Genotype + Age
and interaction terms).

DESeq2 takes in a counts matrix - the number of sequence fragments that
have been assigned to each gene. To statistical determine which genes
are differentially expressed genes, we need to determine the changes
between our conditions compared to the within-condition variability.
DESeq2 uses a negative binomial generalized linear model. The normalized
counts provided by DESeq2 are called *median of ratios*. You should not
use these normalized counts for within sample comparison - They are
primarily used for DE analysis.

First you create a DESeqDataSet (dds) object (contains data + metadata)
and then run DE analysis. In my example data, astrocyte mRNA was
sequenced from 4-month and 2-year old male mice, from 5 different brain
regions, using the Ribotag system [original
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5783200/).

**Note** Check out the tutorial above for information regarding multiple
conditions and interaction terms. In this case, I would actually make a
new variable combining age and region, to keep things simple. Then I
would use that new variable in my design (instead of *Age* as shown
below; see the 2nd code chunk for how to do this).

``` r
dds <- DESeqDataSet(gse, design = ~ RealRegion + Age) #the ~ are your design variables. If you would like to test for an interaction you would do ~Age + RealRegion + Age:RealRegion

dds$Age <- relevel(dds$Age, ref = "4m") #Set a group as reference; default is alphabetical

dge <- DESeq(dds) #run Diff gene exp analysis

round(colSums(counts(dds))/1e6) #number of reads per sample
```

## DESeq2 with STAR Counts Table (or other table)

If you already have a counts table (MUST be raw counts, un-normalized),
then you can create the dds object with the following code:

*Note* You may need to clean up your counts matrix before running the
DESeq2 analysis

``` r
#Code if counts matrix (i.e. from STAR) is the input instead of txi from Salmon
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Group)
```

## Here I am creating a new varaible called *AgeRegion* that I will use as my variable of interest

``` r
dds <- DESeqDataSet(gse, design = ~ AgeRegion)
```

    ## using counts and average transcript lengths from tximeta

    ## Warning in DESeqDataSet(gse, design = ~AgeRegion): some variables in design
    ## formula are characters, converting to factors

``` r
#dds$AgeRegion <- relevel(dds$AgeRegion, ref = "4m_CB") #Set a group as reference; default is alphabetical

dds <- estimateSizeFactors(dds)
```

    ## using 'avgTxLength' from assays(dds), correcting for library size

``` r
counts <- as.data.frame(counts(dds, normalized = TRUE))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
```

    ## using pre-existing normalization factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

## To check for outliers:

``` r
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Extracting Results from DE Analysis

``` r
resultsNames(dds) #view names for the different combinations you tested for - but only for the ones against reference level. Since we have multiple groups (and no single baseline/control), use the contrast option (2nd option below)
```

    ##  [1] "Intercept"                 "AgeRegion_2y_HTH_vs_2y_CB"
    ##  [3] "AgeRegion_2y_MC_vs_2y_CB"  "AgeRegion_2y_SSC_vs_2y_CB"
    ##  [5] "AgeRegion_2y_VC_vs_2y_CB"  "AgeRegion_4m_CB_vs_2y_CB" 
    ##  [7] "AgeRegion_4m_HTH_vs_2y_CB" "AgeRegion_4m_MC_vs_2y_CB" 
    ##  [9] "AgeRegion_4m_SSC_vs_2y_CB" "AgeRegion_4m_VC_vs_2y_CB"

``` r
summary(results(dds, name = "AgeRegion_4m_VC_vs_2y_CB")) #this option won't give you all of your possible comparisons, just those against your reference level
```

    ## 
    ## out of 19363 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 3865, 20%
    ## LFC < 0 (down)     : 3875, 20%
    ## outliers [1]       : 17, 0.088%
    ## low counts [2]     : 2253, 12%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(results(dds, contrast =c("AgeRegion", "2y_VC","4m_VC"))) #use this to obtain all possible comparisons
```

    ## 
    ## out of 19363 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 189, 0.98%
    ## LFC < 0 (down)     : 166, 0.86%
    ## outliers [1]       : 17, 0.088%
    ## low counts [2]     : 3750, 19%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
res_VC <- results(dds, contrast =c("AgeRegion", "2y_VC","4m_VC")) #select results as table, needed for visualizations below
res_MC <- results(dds, contrast =c("AgeRegion", "2y_MC","4m_MC")) 
res_HTH <- results(dds, contrast =c("AgeRegion", "2y_HTH","4m_HTH")) 
res_SSC <- results(dds, contrast =c("AgeRegion", "2y_SSC","4m_SSC")) 
res_CB <- results(dds, contrast =c("AgeRegion", "2y_CB","4m_CB")) 

#make dataframes of all contrasts
res_VC.df <- as.data.frame(res_VC) 
res_MC.df <- as.data.frame(res_MC)
res_CB.df <- as.data.frame(res_CB)
res_HTH.df <- as.data.frame(res_HTH)
res_SSC.df <- as.data.frame(res_SSC)

#all in one step
#res_CB <- as.data.frame(results(dds, contrast =c("AgeRegion", "2y_CB","4m_CB")))

#write.csv(res_VC.df, "AgingAstroTranscriptome_Results_VC.csv", row.names = TRUE) #to export as csv file
```

## Visualizations

### MA Plot

the function plotMA shows the log2 fold changes attributable to a given
variable over the mean of normalized counts for all the samples in the
DESeqDataSet.

``` r
plotMA(res_CB, ylim=c(-2,2))
abline(h=c(-1,1), col="dodgerblue", lwd=2) #add horizontal line
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
par(mfrow=c(1,4), mar=c(2,2,2,2))
plotMA(res_VC, ylim=c(-2,2), main = "VC")
plotMA(res_MC, ylim=c(-2,2), main = "MC")
plotMA(res_HTH, ylim=c(-2,2), main = "HTH")
plotMA(res_CB, ylim=c(-2,2), main = "CB")
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->
\### PCA

There are several different options for making PCA plots in R. Here I
show an example with plotPCA and ggplot. The other popular option is
[PCAtools](https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html)

``` r
vsd <- vst(dds) #1st you need to transform data. I use variance stabilization transformation

data <- plotPCA(vsd, intgroup = "Age", returnData = TRUE)

percent_var <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(x = PC1, y = PC2, color = metadata$RealRegion, shape = metadata$Age)) + 
  geom_point(size = 5) +
  xlab(paste("PC1: ", percent_var[1], "%variance")) +
  ylab(paste("PC2: ", percent_var[2], "%variance")) + 
  theme_classic(base_size = 14) + theme(legend.title = element_blank()) 
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Here is another PCA plot, but with a few more modifications. The
commented lines are not included in the code to generate the plot
(included here as an example, if you wanted to add labels to each
point).

``` r
pca <- ggplot(data, aes(x = PC1, y = PC2, 
                        color = metadata$Age, shape = metadata$RealRegion))+ 
      geom_point(size = 3) + 
      # geom_text_repel(aes(label = 'labelname',
      #                 nudge_x = .5, nudge_y = 0.5,
      #                 fontface = 2,
      #                 check_overlap = TRUE,
      #                 size = 4) +
      xlab(paste("PC1: ", percent_var[1], "%variance")) +
      ylab(paste("PC2: ", percent_var[2], "%variance")) + 
      labs(shape = "Brain Region", color = "Age") +
      scale_shape_discrete(labels=c("Cerebellum", "Hypothalamus", "Motor Cortex",   "Somatosensory", "Visual Cortex")) +
      scale_color_discrete(labels=c("2 year old", "4 month old")) +
      scale_color_manual(values = c("deeppink1","deepskyblue1"))+
      theme_classic(base_size = 12)  +    
      theme(plot.title = element_text(hjust = 0.5, size=14, face="bold", margin = margin(0,0,10,0))) +
      ggtitle("PCA: Aging Astrocyte Transcriptome")
```

    ## Scale for 'colour' is already present. Adding another scale for 'colour',
    ## which will replace the existing scale.

``` r
pca
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
#Code to calculate ellipses around data - not possible with this dataset due to low sample size per group
# pca + stat_ellipse(aes(group = metadata$Region),
#                                   type = "t",
#                                   linetype = 2,
#                                   lwd = 1)
```

### Heatmaps

There are several options for heatmaps as well. I like to use
[pheatmap](https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/)
because the default options are pretty good. You can plot a heatmap with
all your genes, but it takes a long time and is not super informative.
Alternatively, try plotting a heatmap of the most variable genes or a
subset of genes of interest (examples shown below). The color scale is
based on Z-scores of the transformed counts matrix.

``` r
#Example 1: Selecting most variable genes
topVarGenes <- head(order(-rowVars(assay(vsd))),100) # select top 100 most variable genes (across all samples)

pheatmap(assay(vsd)[topVarGenes,]) #basic heatmap
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
sampleInfo <- as.data.frame(colData(dds)[,c("Age","RealRegion")]) #include annotations for plot based on metadata

map1 <- pheatmap(assay(vsd)[topVarGenes,], 
           annotation_col = sampleInfo,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           clustering_distance_rows = "euclidean",
           show_rownames = FALSE, 
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols =  FALSE)
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
#Example 2: Selecting significant DEGs for just Cerebellum 2y vs 4m 
sigGenes <- subset(as.data.frame(res_CB), padj <=0.01) 
sigGenes <- subset(sigGenes, abs(log2FoldChange) > 1)

#Create custom colors for your labels 
my_colour = list(
  Age = c("2y" = "steelblue2", "4m" = "brown2"),
  RealRegion = c( CB = "cyan", VC = "mediumpurple1", MC = "orange1", SSC = "orchid", HTH = "aquamarine"))

map2 <- pheatmap(assay(vsd)[row.names(sigGenes),], 
           annotation_col = sampleInfo, 
           annotation_colors = my_colour,
           scale = "row",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           clustering_distance_rows = "euclidean",
           show_rownames = FALSE, 
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols =  TRUE,
           cutree_cols = 2)
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->

``` r
cluster <- as.data.frame(sort(cutree(map2$tree_row, k=3))) #Export clusters based on heatmap.

#to export figures
#ggsave("fig1.pdf", height = 5, width = 7, unit = "in", dpi = 500)
```

Another type of visualization option is to plot a heatmap of
sample-to-sample distances. A heatmap of this distance matrix gives us
an overview over similarities and dissimilarities between samples. We
have to provide a hierarchical clustering hc to the heatmap function
based on the sample distances, or else the heatmap function would
calculate a clustering based on the distances between the rows/columns
of the distance matrix.

``` r
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Age, vsd$RealRegion, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### Volcano Plots

The plots below are made using [Enhanced Volcano
package](https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)

``` r
EnhancedVolcano(res_CB,
    lab = rownames(res_CB),
    x = 'log2FoldChange',
    y = 'pvalue')
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
EnhancedVolcano(res_CB,
    lab = rownames(res_CB),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = '4m vs 2y CB',
    pCutoff = 0.01,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 6.0)
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
res_CB.df <- merge(res_CB.df, annotation, by = 0)
res_CB.df <- res_CB.df[!is.na(res_CB.df$symbol),]

#Custom shape and color 
keyvals <- ifelse(
           res_CB.df$log2FoldChange < -2.5, 'royalblue',
           ifelse(res_CB.df$log2FoldChange > 2.5, 'gold','black'))

keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

Volcano1 <- EnhancedVolcano(res_CB.df,
            lab = res_CB.df$symbol,
            x = 'log2FoldChange',
            y = 'pvalue',
            selectLab = res_CB.df$symbol[which(names(keyvals) %in% c('high', 'low'))],
            xlab = bquote(~Log[2]~ 'fold change'),
            title = '4m vs 2y CB',
            subtitle = '',
            legendLabels=c('Not sig.',~Log[2]~ 'fold change >|1|','padj<0.01','Both'),
            pCutoff = 0.01,
            FCcutoff = 1.0,
            pointSize = 3.5,
            labSize = 4.5,
            shape = c(6, 4, 2, 11),
            colCustom = keyvals,
            colAlpha = 1,
            legendPosition = 'left',
            legendLabSize = 15,
            legendIconSize = 5.0,
            drawConnectors = TRUE,
            widthConnectors = 1.0,
            colConnectors = 'black',
            arrowheads = FALSE,
            gridlines.major = TRUE,
            gridlines.minor = FALSE,
            border = 'partial',
            borderWidth = 1.5,
            borderColour = 'black')


Volcano2 <- EnhancedVolcano(res_CB.df,
            lab = res_CB.df$symbol,
            x = "log2FoldChange",
            y = "pvalue",
            pCutoff = 0.01,
            FCcutoff = 1.5,
            pointSize = c(ifelse(res_CB.df$log2FoldChange>2, 2, 1)),
            labSize = 2.0,
            shape = c(6, 6, 19, 16),
            title = "DESeq2 results",
            subtitle = "Differential expression",
            caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 0.01"),
            legendPosition = "right",
            legendLabSize = 14,
            legendLabels=c('Not sig.',~Log[2]~ 'fold change >|1|','padj<0.01','Both'),
            col = c("grey30", "forestgreen", "royalblue", "red2"),
            colAlpha = 0.9,
            drawConnectors = TRUE,
            hline = c(10e-8),
            widthConnectors = 0.5)

Volcano1
```

    ## Warning: ggrepel: 536 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->

``` r
Volcano2
```

    ## Warning: It is deprecated to specify `guide = FALSE` to remove a guide. Please
    ## use `guide = "none"` instead.

    ## Warning: ggrepel: 192 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-17-4.png)<!-- -->

### Plot Individual Counts

You can modify these plots with the usual ggplot2. Below is a just a
basic example.

``` r
d <- plotCounts(dds, gene=which.min(res_CB$padj), intgroup="AgeRegion", 
                returnData=TRUE)

ggplot(d, aes(x=AgeRegion, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
#How to show only a subset of groups
d <- plotCounts(dds[ , dds$AgeRegion %in% c("2y_CB","4m_CB")], gene="ENSMUSG00000002985", intgroup="AgeRegion", returnData=TRUE)

#d$SampleNames <- metadata$simpleID


ggplot(d, aes(x = AgeRegion, y = count, fill = AgeRegion)) +
            geom_boxplot()+
            #geom_point(aes(color = SampleNames), size = 2,position=position_jitter(w = 0.1,h = 0)) +
            #geom_text_repel(aes(label = d$SampleNames)) + 
            ggtitle("RNA: APOE") +
            ylab("Normalized Count (Log10 Scale)")+
            expand_limits(y = 100000)+
            scale_y_continuous(trans = 'log10',
                               expand = expansion(mult = c(.1, 0.1)), labels = comma)+
            scale_fill_manual(values = c( "grey60","darkorange", "black",  "firebrick3"))+
            theme_classic(base_size = 16)+
            theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle = 45, hjust =1),
                  #axis.title.y=element_text(size=14),
                  axis.title.x=element_blank(),
                  legend.position="right") 
```

![](AgingAstrocyteTranscriptome_Tutorial_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### DESeq2 Report Tool

This
[function](https://www.rdocumentation.org/packages/regionReport/versions/1.6.5/topics/DESeq2Report)
creates example graphs and tables with your data.

``` r
dir.create("DESeq2Report-example", showWarnings = FALSE, recursive = TRUE)

## Generate the HTML report
report <- DESeq2Report(dge, "DESeq2-example", c("Age", "RealRegion"),
    outdir = "DESeq2Report-example"
)

if (interactive()) {
    ## Browse the report
    browseURL(report)
}
```

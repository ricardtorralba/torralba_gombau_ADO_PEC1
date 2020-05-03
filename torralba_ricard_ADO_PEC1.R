params <-
list(samples = 12L, groups = 4L, p.value = 0.05, log.fc = 1L)

## ----installBioC, message=FALSE, warning=FALSE, eval=FALSE-----------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install()


## ----installPackages, message=FALSE, warning=FALSE, eval=FALSE-------------------------------
## install.packages("knitr")
## install.packages("ggplot2")
## install.packages("ggrepel")
## install.packages("BiocManager")
## BiocManager::install("oligo")
## BiocManager::install("pd.mogene.2.0.st")
## BiocManager::install("arrayQualityMetrics")
## BiocManager::install("limma")
## BiocManager::install("genefilter")
## BiocManager::install("mogene20sttranscriptcluster.db")
## BiocManager::install("org.Mm.eg.db")
## BiocManager::install("ReactomePA")
## BiocManager::install("reactome.db")


## ----setup, include=FALSE--------------------------------------------------------------------
library(knitr)
opts_chunk$set(comment = NA,
               echo = T,
               prompt = T,
               message = F, 
               warning = F,
               cache = T)


## ----libraries, include = F, message = F, warning = F----------------------------------------
# Load packages
library(ggplot2)
library(ggrepel)
library(Biobase)
library(oligo)
library(genefilter)
library(limma)
library(ReactomePA)
library(reactome.db)
library(arrayQualityMetrics)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(org.Mm.eg.db)


## ---- echo = F-------------------------------------------------------------------------------
targets = read.csv2("data/targets_cor.csv",
                    header = T,
                    sep = ";")
kable(targets,
      caption = "Group at which each of the samples belongs")


## ---- echo = F, results='hide'---------------------------------------------------------------
celFiles = list.celfiles("./data",
                         full.names = T)
my.targets = read.AnnotatedDataFrame("data/targets_cor.csv",
                                     header = T,
                                     row.names = 1,
                                     sep = ";")
rawData = read.celfiles(celFiles,
                        phenoData = my.targets)
rownames(pData(rawData)) = my.targets@data$ShortName
colnames(rawData) = rownames(pData(rawData))


## ---- echo = F-------------------------------------------------------------------------------
arrayQualityMetrics(rawData,
                    outdir = "./results/QCDir.Raw",
                    force = T)


## ---- echo = F, fig.align = 'center',fig.cap = "Quality analysis produced by the arrayQualityMetrics package on the raw data."----
include_graphics("results/quality_raw.PNG")


## ---- echo = F, fig.width = 7.5, fig.height = 4.5, fig.align = 'center', fig.cap = "Intensity distribution of the arrays raw data."----
boxplot(rawData, 
        cex.axis = .8,
        las = 2,
        which = "all",
        col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
        main = "Distribution of raw instensity values")


## ----rawpca, echo = F, fig.width = 7, fig.height = 4.5, fig.align = 'center', fig.cap = "Principal Component Analysis for the raw data."----
plotPCA = function (datos, labels, factor, title, scale, colours, size = 1.5, glineas = 0.25) {
  data = prcomp(t(datos), scale = scale)
  dataDf = data.frame(data$x)
  Group = factor
  loads = round(data$sdev^2/sum(data$sdev^2)*100, 1)
  p1 = ggplot(dataDf,aes(x = PC1, y = PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5, max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels), segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1", loads[1], "%")), y = c(paste("PC2", loads[2], "%"))) +  
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colours)
  }

plotPCA(exprs(rawData),
        labels = targets$ShortName,
        factor = targets$Group,
        title = "Raw data PCA",
        scale = F,
        size = 3,
        colours = c("red", "blue", "green", "yellow"))


## ---- echo = F, results = 'hide'-------------------------------------------------------------
pca_raw = prcomp(t(exprs(rawData)))
loads_raw = round(pca_raw$sdev^2/sum(pca_raw$sdev^2)*100, 1)


## ---- echo = F, results = "hide"-------------------------------------------------------------
data_rma = rma(rawData)


## ---- echo = F-------------------------------------------------------------------------------
arrayQualityMetrics(data_rma,
                    outdir = "./results/QCDir.Norm",
                    force = T)


## ---- echo = F, fig.align = 'center',fig.cap = "Quality analysis for the normalized data produced by the arrayQualityMetrics package."----
include_graphics("results/quality_norm.PNG")


## ---- echo = F, fig.width = 7.5, fig.height = 4.5, fig.align = 'center', fig.cap = "Intensity distribution of the arrays normalized data."----
boxplot(data_rma, 
        cex.axis = .8,
        las = 2,
        which = "all",
        col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)),
        main = "Distribution of normalized instensity values")


## ---- echo = F, fig.width = 7, fig.height = 4.5, fig.align = 'center', fig.cap = "Principal Component Analysis for the normalized data."----
plotPCA(exprs(data_rma),
        labels = targets$ShortName,
        factor = targets$Group,
        title = "Normalized data PCA",
        scale = F,
        size = 3,
        colours = c("red", "blue", "green", "yellow"))


## ---- echo = F, results = "hide"-------------------------------------------------------------
pca_norm = prcomp(t(exprs(data_rma)))
loads_norm = round(pca_norm$sdev^2/sum(pca_norm$sdev^2)*100, 1)
norm_reduction = loads_raw[1] - loads_norm[1]


## ---- echo = F-------------------------------------------------------------------------------
annotation(data_rma) = "mogene20sttranscriptcluster.db"
filtered = nsFilter(data_rma,
                    require.entrez = T,
                    remove.dupEntrez = T,
                    var.filter = T,
                    var.func = IQR,
                    var.cutoff = .75,
                    filterByQuantile = T,
                    feature.exclude = "^AFFX")
data_filtered = filtered$eset
after_filtered_genes = dim(exprs(data_filtered))[1]


## ---- echo = F-------------------------------------------------------------------------------
designMat = model.matrix(~0 + Group,
                         pData(data_filtered))
colnames(designMat) = c("KO.NO", "KO.YES", "WT.NO", "WT.YES")
designMat[,]


## ---- echo = F-------------------------------------------------------------------------------
contrastMat = makeContrasts (KOvsWT.NO = KO.NO - WT.NO,
                             KOvsWT.YES = KO.YES - WT.YES,
                             INT = (KO.NO - WT.NO) - (KO.YES - WT.YES),
                             levels = designMat)
contrastMat


## ---- echo = F-------------------------------------------------------------------------------
fit = lmFit(data_filtered, designMat)
fit.main = contrasts.fit(fit, contrastMat)
fit.main = eBayes(fit.main)


## ---- echo = F-------------------------------------------------------------------------------
topTab_KOvsWT.NO = topTable(fit.main,
                            number = nrow(fit.main),
                            coef = "KOvsWT.NO",
                            adjust = "fdr")
kable(head(topTab_KOvsWT.NO, 10),
      caption = "Genes most differently expressed when PRDM16 transcription factor is deleted and Rosiglitazone is not present")


## ---- echo = F-------------------------------------------------------------------------------
topTab_KOvsWT.YES = topTable(fit.main,
                             number = nrow(fit.main),
                             coef = "KOvsWT.YES",
                             adjust = "fdr")
kable(head(topTab_KOvsWT.YES, 10),
      caption = "Genes most differently expressed when PRDM16 transcription factor is deleted and Rosiglitazone is present")


## ---- echo = F-------------------------------------------------------------------------------
topTab_INT = topTable(fit.main,
                      number = nrow(fit.main),
                      coef = "INT",
                      adjust = "fdr")
kable(head(topTab_INT, 10),
      caption = "Genes most differently expressed between the previous two comparisons")


## ---- echo = F-------------------------------------------------------------------------------
geneSymbols = select(mogene20sttranscriptcluster.db,
                     rownames(fit.main),
                     c("SYMBOL"))
SYMBOLS = geneSymbols$SYMBOL


## ---- echo = F, fig.align = 'center', fig.cap = "PRDM16 KO and Rosiglitazone not present"----
volcanoplot(fit.main,
            coef = 1,
            highlight = 5,
            names = SYMBOLS,
            main = paste("Differentially expressed genes",
                         colnames(contrastMat)[1],
                         sep = " "))
abline(v = c(-1,1))


## ---- echo = F, fig.align = 'center', fig.cap = "PRDM16 WT and Rosiglitazone present"--------
volcanoplot(fit.main,
            coef = 2,
            highlight = 5,
            names = SYMBOLS,
            main = paste("Differentially expressed genes",
                         colnames(contrastMat)[2],
                         sep = " "))
abline(v = c(-1,1))


## ---- echo = F, fig.align = 'center', fig.cap = "Interaction between factors"----------------
volcanoplot(fit.main,
            coef = 3,
            highlight = 5,
            names = SYMBOLS,
            main = paste("Differentially expressed genes",
                         colnames(contrastMat)[3],
                         sep = " "))
abline(v = c(-1,1))


## ---- echo = F-------------------------------------------------------------------------------
annotatedTopTable = function(topTab, anotPackage){
  topTab = cbind(PROBEID = rownames(topTab),
                 topTab)
  myProbes = rownames(topTab)
  thePackage = eval(parse(text = anotPackage))
  geneAnots = select(thePackage,
                     myProbes,
                     c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab = merge(x = geneAnots,
                          y = topTab,
                          by.x = "PROBEID",
                          by.y = "PROBEID")
  return(annotatedTopTab)
}


## ---- echo = F-------------------------------------------------------------------------------
topAnnotated_KOvsWT.NO = annotatedTopTable(topTab_KOvsWT.NO,
                                           anotPackage = "mogene20sttranscriptcluster.db")
topAnnotated_KOvsWT.YES = annotatedTopTable(topTab_KOvsWT.YES,
                                            anotPackage = "mogene20sttranscriptcluster.db")
topAnnotated_INT = annotatedTopTable(topTab_INT,
                                     anotPackage = "mogene20sttranscriptcluster.db") 


## ---- echo = F-------------------------------------------------------------------------------
low_adjusted_pvalue_KOvsWT.NO = head(order(topAnnotated_KOvsWT.NO$adj.P.Val), 10)
kable(topAnnotated_KOvsWT.NO[low_adjusted_pvalue_KOvsWT.NO, c(1:4,9)],
      caption = "Annotated most differently expressed genes when PRDM16 transcription factor is deleted and Rosiglitazone is not present")


## ---- echo = F-------------------------------------------------------------------------------
low_adjusted_pvalue_KOvsWT.YES = head(order(topAnnotated_KOvsWT.YES$adj.P.Val), 10)
kable(topAnnotated_KOvsWT.YES[low_adjusted_pvalue_KOvsWT.YES, c(1:4,9)],
      caption = "Annotated most differently expressed genes when PRDM16 transcription factor is deleted and Rosiglitazone is present")


## ---- echo = F-------------------------------------------------------------------------------
low_adjusted_pvalue_INT = head(order(topAnnotated_INT$adj.P.Val), 10)
kable(topAnnotated_INT[low_adjusted_pvalue_INT ,c(1:4,9)],
      caption = "Annotated most differently expressed genes between the previous two comparisons")


## ---- echo = F-------------------------------------------------------------------------------
res = decideTests(fit.main,
                  method = "separate",
                  adjust.method = "fdr",
                  p.value = params$p.value,
                  lfc = params$log.fc)
sum.res.rows = apply(abs(res),1,sum)
res.selected = res[sum.res.rows != 0,]
kable(summary(res), caption = "Differently expressed genes in each comparison.")


## ---- echo = F, fig.align = 'center', fig.cap = "Genes in common between the three comparisons"----
vennDiagram(res.selected[,1:3], cex = .9)


## ---- echo = F-------------------------------------------------------------------------------
listOfTables = list(KOvsWT.NO = topTab_KOvsWT.NO,
                    KOvsWT.YES = topTab_KOvsWT.YES,
                    INT = topTab_INT)
listOfSelected = list()
for (i in 1:length(listOfTables)){
  topTab = listOfTables[[i]]
  whichGenes = topTab["adj.P.Val"]< params$p.value
  selectedIDs = rownames(topTab)[whichGenes]
  EntrezIDs = select(mogene20sttranscriptcluster.db,
                     selectedIDs,
                     c("ENTREZID"))
  EntrezIDs = EntrezIDs$ENTREZID
  listOfSelected[[i]] = EntrezIDs
  names(listOfSelected)[i] = names(listOfTables)[i]
}


## ---- echo = F-------------------------------------------------------------------------------
kable(sapply(listOfSelected, length), caption = "Genes used in the enrichment analysis for each of the comparisons.")


## ---- echo = F-------------------------------------------------------------------------------
universe_genesGO = mappedkeys(org.Mm.egGO)
universe_genes2KEGG = mappedkeys(org.Mm.egPATH)
universe_genes = union(universe_genesGO, universe_genes2KEGG)


## ---- echo = F-------------------------------------------------------------------------------
universe = universe_genes
gene_enrichment = list()
results_enrich = list()

for(i in 1:length(listOfSelected)){
  genesIn = listOfSelected[[i]]
  enrich.result = enrichPathway(gene = genesIn,
                                pvalueCutoff = params$p.value,
                                readable = T,
                                pAdjustMethod = "BH",
                                organism = "mouse",
                                universe = universe)
  results_enrich[[i]] = enrich.result
  enrich.df = as.data.frame(enrich.result[,c(2:7,9)])
  gene_enrichment[[i]] = enrich.df
}
names(results_enrich) = names(listOfSelected)
names(gene_enrichment) = names(listOfSelected)
gene_enrichment


## ---- echo = F-------------------------------------------------------------------------------
barplot(results_enrich$KOvsWT.NO,
        showCategory = 15,
        font.size = 8,
        title = "Reactiome Pathway Analysis")


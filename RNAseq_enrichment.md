```
library("DESeq2")
library("dplyr")
library("tidyverse")
library("data.table")
library("ggplot2")
library("enrichR")
library(ComplexHeatmap)
```
## upload metadata file and feature counts table

```
colData <- read.csv("metadata_2groups.csv", sep=";") #metadata

countData <- fread(paste("featureCounts.txt", sep = ""), 
                   header = T, sep = "\t", check.names = FALSE) #feature counts
```

## DESeq2

```
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ Condition)

dds <- DESeq(dds)

res <- results(dds)

resultsNames(dds)

test <- as.data.frame((resOrdered <- res[order(res$pvalue), ]))

test.sign<- subset(test, padj <= 0.001)

vsd <- vst(dds)

vvt <-assay(vsd)

test.out<- test.sign
test.out$gene<- row.names(test.sign)

write.table(test.out,"results_deseq2_FDR0_001_1696genes_highVSneg.csv",quote=F, sep=",")
```

## Generating enrichment analysis

```
test.sign.down<- test.sign[which(test.sign$log2FoldChange>0),]

test.sign.up<- test.sign[which(test.sign$log2FoldChange<0),]

genenames<- row.names(test.sign.down)

genels = c()

strGenes <- as.character(genenames)
splitGenes = strsplit(strGenes, "[;]")

splitGenes <-matrix(unlist(splitGenes))

genesym <- unlist(splitGenes[!duplicated(splitGenes)])


dbs <- listEnrichrDbs()
dbs <- "KEGG_2019_Human"
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
enriched <-enrichr(genesym, databases = dbs)


KEGG_enrich<- as.data.frame(enriched[["KEGG_2019_Human"]])

KEGG_enrich.f<- subset(KEGG_enrich, Adjusted.P.value< 0.05)

gos <- KEGG_enrich.f

gos <- gos[order(-gos$P.value), ]

gos$Term <- factor(gos$Term, levels=gos$Term)

ggplot(gos, aes(x=Term, y=Adjusted.P.value , label=Adjusted.P.value)) + 
 geom_bar(stat='identity', width=.2,position="dodge")  +
 coord_flip()

dbs <- listEnrichrDbs()

dbs <- "GO_Biological_Process_2018"
options(enrichR.base.address="https://amp.pharm.mssm.edu/Enrichr/")
enriched <-enrichr(genesym, databases = dbs)


GO_enrich<- as.data.frame(enriched[["GO_Biological_Process_2018"]])

GO_enrich.f<- subset(GO_enrich, Adjusted.P.value< 0.05)

gos <- GO_enrich.f 

gos <- gos[order(-gos$Adjusted.P.value), ]
gos$Term <- factor(gos$Term, levels=gos$Term)

ggplot(gos, aes(x=Term, y=Adjusted.P.value , label=Adjusted.P.value)) + 
 geom_bar(stat='identity', width=.1,position="dodge")  +
 coord_flip()
```

## Hierarchical clustering
```
data<-as.data.frame(vvt)

data$ID<- row.names(data)

test.sign$ID<- row.names(test.sign)

data.all<-merge(data, test.sign, by.x = "ID", by.y = "ID")

heat.data<- data.all[,c(2:21)]
row.names(heat.data)<- data.all$ID

z.mat <- t(scale(t(heat.data), center=TRUE, scale=TRUE))

column_ha3 = HeatmapAnnotation(Groups = colData$Condition, IGHV= colData$IGHV,FISH= colData$FISH, #Patient= colData$Patient_ID,
                               col = list(Groups = c("PD-1 high" = "green", "PD-1 Intermediate" = "gold", "PD-1 negative" = "red"),
                                          IGHV = c("Mutated" = "gray", "Unmutated" = "black", "N/A"="white"),
                                          FISH= c("del13q"='lightblue3', "del13q, del11q"="lightgreen", "del13q, trisomy 12"= "purple", "Normal"= "darkorange2", "Trisomy 12"="blue")),na_col="white"
)

Heatmap(as.matrix(z.mat), clustering_distance_columns = "euclidean",
clustering_method_columns = "complete", 
        top_annotation = column_ha3, row_names_gp = gpar(fontsize = 4), show_row_names = FALSE) #column_km = 2)

```

## Hierarchical clustering based on gene list

```
genes<- fread(paste("geneList.txt", sep = ""), 
            header = T, sep = "\t", check.names = FALSE) #select specific genes



countData<-left_join(genes, countData,  by="Geneid")
head(countData)
row.names(countData)<- countData$Geneid
names(countData)  <- colnames(countData)  %>% str_replace_all("_sorted.bam", "")
colnames(countData)



colData<-as.data.frame(colData)
countData<-as.data.frame(countData)
row.names(countData)<- countData$Geneid
list<-as.character(colData$file_id)
countData<- countData[,list]




#deseq2
dds <- DESeqDataSetFromMatrix(countData = countData, 
                                       colData = colData, design = ~ Condition)

dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)


test <- as.data.frame((resOrdered <- res[order(res$pvalue), ]))
test.sign<- subset(test, padj <= 0.05)

vsd <- varianceStabilizingTransformation(dds)
vvt <-assay(vsd)


test.out<- test.sign
test.out$gene<- row.names(test.sign)

##Heatmap
data<-as.data.frame(vvt)
head(data)
data$ID<- row.names(data)

test.sign$ID<- row.names(test.sign)

data.all<-data
data.all<-merge(data, test.sign, by.x = "ID", by.y = "ID")
colnames(data.all)
heat.data<- data.all[,c(2:21)]
row.names(heat.data)<- data.all$ID

z.mat <- t(scale(t(heat.data), center=TRUE, scale=TRUE))



colData
library(ComplexHeatmap)
column_ha3 = HeatmapAnnotation(Groups = colData$Condition, IGHV= colData$IGHV,FISH= colData$FISH, #Patient= colData$Patient_ID,
                               col = list(Groups = c("PD-1 high" = "green", "PD-1 Intermediate" = "gold", "PD-1 negative" = "red"),
                                          IGHV = c("Mutated" = "gray", "Unmutated" = "black", "N/A"="white"),
                                          FISH= c("del13q"='lightblue3', "del13q, del11q"="lightgreen", "del13q, trisomy 12"= "purple", "Normal"= "darkorange2", "Trisomy 12"="blue")
                                          
                               ),
                               na_col="white"
)

Heatmap(as.matrix(z.mat), clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        top_annotation = column_ha3, row_names_gp = gpar(fontsize = 4), show_row_names = TRUE) #column_km = 2)

```

## Volcano
```
res_tableOE<-res

threshold_OE <- res_tableOE$padj <= 0.001
res_tableOE$threshold <- threshold_OE 
res_tableOE<-as.data.frame(res_tableOE)

ggplot(res_tableOE) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  



## Sort by ordered padj
res_tableOE_ordered <- res_tableOE[order(res_tableOE$log2FoldChange), ] 

## Create a column to indicate which genes to label
res_tableOE_ordered$genelabels <- ""
sign<-c("CD5", "CXCR4", "PDCD1", "MKI67", "CD27", "CCL3", "CCL4", "BACH2", 
"RPTOR", "MTORC", "FOXO1")

for(i in 1:nrow(res_tableOE_ordered)){
       if(row.names(res_tableOE_ordered)[i] %in% sign){
         res_tableOE_ordered$genelabels[i] <- rownames(res_tableOE_ordered)[i]
       } else { 
         res_tableOE_ordered$genelabels[i] <- ""
         }
       }

volc<-ggplot(res_tableOE_ordered) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  


volc+geom_text_repel(data=res_tableOE_ordered, max.overlaps = Inf,
                     aes(x = log2FoldChange, y = -log10(padj),label=genelabels),
                     point.padding = 0.2,
                     nudge_x = .15,
                     nudge_y = .5,
                     segment.linetype = 6,
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc")))

volc<-ggplot(res_tableOE_ordered) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=ifelse(genelabels %in% sign, "red", "black"))) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

volc2<-ggplot(res_tableOE_ordered) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), color=threshold),
             shape=ifelse(res_tableOE_ordered$genelabels %in% sign, 24, 19), 
             fill="red", size=4, alpha=0.5) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


volc2+geom_text_repel(data=res_tableOE_ordered, max.overlaps = Inf,
                     aes(x = log2FoldChange, y = -log10(padj),label=genelabels),
                     point.padding = 0.2,
                     nudge_x = .15,
                     nudge_y = .5,
                     segment.linetype = 6,
                     segment.curvature = -1e-20,
                     arrow = arrow(length = unit(0.015, "npc")))
```

### GSE
```
test <- as.data.frame((resOrdered <- res[order(res$pvalue), ]))
test.sign<- subset(test, padj <= 0.001)

vsd <- vst(dds)
vvt <-assay(vsd)


test.sign.down<- test.sign[which(test.sign$log2FoldChange>0),]
test.sign.up<- test.sign[which(test.sign$log2FoldChange<0),]


library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# we want the log2 fold change 
#original_gene_list <- df$log2FoldChange
original_gene_list <-test.sign.up$log2FoldChange
# name the vector
names(original_gene_list) <- row.names(test.sign.up)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


#gene.df <- bitr(gene, fromType = "ENTREZID",
#                toType = c("ENSEMBL", "SYMBOL"),
#                OrgDb = org.Hs.eg.db)

ego3 <- gseGO(geneList     = gene_list,
              keyType = "SYMBOL", 
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

require(DOSE)
dotplot(ego3, showCategory=10, split=".sign")

ridgeplot(ego3) + labs(x = "enrichment distribution")


gseaplot(ego3, by = "all", title = ego3$Description[1], geneSetID = 1)
```

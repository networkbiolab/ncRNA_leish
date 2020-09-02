
#----------------------------ReadME use DESeq2----------------------------------------
## try http:// if https:// URLs are not supported

# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")


# cts -> is the matrix counts
# coldata -> is the sample information

# NOTE: "cts" in the first column must have "gene_id" and not be empty!
# NOTE2: When import "coldata", you must specified separation by \t

# pasCts <- read.csv("a_n0_m1_first_run.csv",sep="\t",row.names="gene_id")

#---------------------------Run DESeq2 ------------------------------------------------

library("DESeq2")
#args = commandArgs(trailingOnly=TRUE)
cts <- as.matrix(read.csv("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/LmajFriedlin_htseq_data/LmajFreidlin_Tritryp_final2.htseq",
                          row.names="gene_id", 
                          header = T,
                          sep = "\t"))
coldata <- read.csv("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/experiments/metada_LmajFriedlin_cut.csv", 
                    header = T,
                    sep = "\t",
                    row.names = 1) 

# With the count matrix, cts, and the sample information, coldata, 
# we can construct a DESeqDataSet:

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~Class)



# CONDITIONS FROM SAMPLE CONDITION FILE
#dds$Class <- factor(dds$Class, levels = c("amastigote","metacyclic_promastigote","procyclic_promastigote"))
#dds$Class <- factor(dds$Class, levels = c("amastigote","promastigote"))
dds$Class <- factor(dds$Class, levels = c("Metacyclic_promastigote","Procyclic_promastigote"))

#run Deseq
dds <- DESeq(dds)

#---------------------------Results all vs all -----------------------------------------------
#obtain Table
res <- results(dds)
write.table(as.data.frame(res), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LmajFriedlin_tritryp/LmajFriedlin_tritryp.csv",
            sep = "\t")
head(results(dds, tidy=TRUE))


#Summary of differential gene expression
summary(res) 

#Sort summary list by p-value

res <- res[order(res$padj,res$log2FoldChange),]
#q=head(res, n=500)
write.table(as.data.frame(res), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LmajFriedlin_tritryp/LmajFriedlin_tritryp_Sorted.csv",
            sep = "\t")
head(results(dds, tidy=TRUE))



#PCA plot all conditions 
library(ggplot2)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Class")
vsdata
dist1<-dist(t(cts))
hclust1 = hclust(dist1)
plot(hclust1)

# --------------------------COMPARATIVE VS lifestage----------------------------
res2 <- results(dds, contrast=c("Class", "amastigote", "promastigote"))
write.table(as.data.frame(res2), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LdonBPK282_tritryp/development_comparison_LdonBPK282_amastigote_promastigote.csv",
            sep = "\t")
topgenes_2 <- head(rownames(subset(res2, padj < 0.05 & log2FoldChange <= -1.5 |padj < 0.05 &
                                   log2FoldChange >= 1.5)),1000)

write.table(as.data.frame(topgenes_2), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LdonBPK282_tritryp/topgenes_development_comparison_LdonBPK282_amastigote_promastigote.csv",
            sep = "\t",row.names = F)



res3 <- results(dds, contrast=c("Class","amastigote","procyclic_promastigote"))
write.table(as.data.frame(res3), 
          file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_ncRNA/development_comparison_LbraM2903_amastigote_procyclic.csv",
          sep = "\t")

topgenes_3 <- head(rownames(subset(res3, padj < 0.05 & log2FoldChange <= -1.5 |padj < 0.05 &
                                   log2FoldChange >= 1.5)),1000)

write.table(as.data.frame(topgenes_3), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LbraM2903_ncRNA/topgenes_development_comparison_LbraM2903_amastigote_procyclic.csv",
            sep = "\t",row.names = F)


res4 <- results(dds, contrast=c("Class","Metacyclic_promastigote","Procyclic_promastigote"))

write.table(as.data.frame(res4), 
          file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LmajFriedlin_tritryp/development_comparison_LmajFriedlin_metacyclic_procyclic.csv",
          sep="\t")

topgenes_4 <- head(rownames(subset(res4, padj < 0.05 & log2FoldChange <= -1.5 |padj < 0.05 &
                                     log2FoldChange >= 1.5)),1000)

write.table(as.data.frame(topgenes_4), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LmajFriedlin_tritryp/topgenes_development_comparison_LmajFriedlin_metacyclic_procyclic.csv",
            sep = "\t",row.names = F)


#---------------------------Normalized count -----------------------------------------------


#dds1 <- estimateSizeFactors(dds)
#dds2<-counts(dds1, normalized=TRUE)
#write.table(as.data.frame(dds2), file = "norm_reads_LbraM2903.csv", sep = "\t")

resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds,normalized =TRUE)), 
                 by = 'row.names', sort = FALSE)

write.table(resdata,file = "/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LmajFriedlin_tritryp/norm_reads_LmajFriedlin_tritryp.csv",
            sep = "\t", row.names = F)











#-----------------------------Volcano Plot ----------------------------------------------------------------

# volcano plot 
library(EnhancedVolcano)
volcano<-EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = T,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                xlim = c(-5, 5),
                ylim = c(0,100),
                title =NULL)




#-----------------------------PLOT Heatmap of DE genes-------------------------------------------------


#metadata
my_sample_col2<- data.frame(rownames(coldata), coldata$Class)
colnames(my_sample_col2) <- c("run","Class")
rownames(my_sample_col2)<- my_sample_col2$run
my_sample_col2$run<- NULL
#my_colour = list(Class=c(amastigote="#1a5276",promastigote="#ff7f00"))
#my_colour = list(Class=c(amastigote="#1a5276",metacyclic_promastigote="#ff7f00",procyclic_promastigote ="#d4ac0d"))
my_colour = list(Class=c(Metacyclic_promastigote="#ff7f00",Procyclic_promastigote ="#d4ac0d"))

#heatmap
library(pheatmap)
topgenes <- head(rownames(subset(res, padj < 0.05 & log2FoldChange <= -1.5 |padj < 0.05 &
                                   log2FoldChange >= 1.5)),1000)
write.table(as.data.frame(topgenes), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/DE_analysis/LmajFriedlin_tritryp/topgenes_LmajFriedlin_tritryp.csv",
            sep = "\t",row.names = F)

rld <- rlog(dds, blind=FALSE)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
heatMap1=pheatmap(mat, 
         annotation_colors = my_colour,
         annotation_col=my_sample_col2,
         fontsize = 6,
         show_rownames = F)
         #, 
         #filename = "~/Desktop/01_htseq_outs/DE_analysis/Heatmap_DE_ncRNA_LdonBPK282.pdf" )



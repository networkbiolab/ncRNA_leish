library(fgsea)
library(data.table)
library(ggplot2)



#cargo los datos
paths_data<-read.delim("~/Escritorio/WGCNA/KEGG_annotation/GMT_Files/LbrM2903_parsed.gmt",
                                 header = F, 
                                 sep = "\t")

#selecciono los nombres
nams <- paths_data[, 1]
#eliminar la columna de los nombres y la que no ocupo
paths_data <- paths_data[, -c(1:2)]
#hacer la lista
pathways <- split(paths_data, seq_len(nrow(paths_data)))
#eliminar los elementos en blanco 
pathways <- lapply(pathways, function(x) x[x != ""])
#agragar los nombres de a la lista
names(pathways) <- nams
str(head(pathways))

ranks <- read.table('~/Escritorio/WGCNA/LbraM2903_median.htseq',
                    header=TRUE, 
                    sep = '\t',
                    colClasses = c("character", "numeric"))
ranks <- setNames(ranks$median, ranks$gene_id)
str(ranks)
fgseaRes <- fgsea(pathways = pathways, 
                  stats = ranks,
                  nperm=10000)

head(fgseaRes[order(pval), ])

sum(fgseaRes[, padj < 0.05])



topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam = 0.05)



barplot(sort(ranks, decreasing = T))

plotEnrichment(pathways, ranks)


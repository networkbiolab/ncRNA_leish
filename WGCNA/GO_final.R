library(goseq)
library(GO.db)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(grid)


##########PLOTS#############
#crea un tema para los plots

theme_set(
  theme_pubr() +
    theme(legend.position = "right")
)

#########LbraM2903#########
#read and parse GAF file
gaf_Lbra<-read.delim("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/seq_Leishmania_Companion_GeneID/AA_Tritryp/gaf_files_Release35/TriTrypDB-35_LbraziliensisMHOMBR75M2903_GO.gaf", 
                     header = F, sep="\t")

final_Lbra<-data.frame(gaf_Lbra$V2, gaf_Lbra$V5)
final_Lbra<-final_Lbra[-1,]

# read Gene Length file and convert to numeric vector
length_Lbra= read.delim("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/seq_Leishmania_Companion_GeneID/CDS_Tritryp/GC_content/TriTrypDB-35_LbraziliensisMHOMBR75M2903_AnnotatedCDSs.fasta.txt",
                        header = T, sep = "\t")
flLbra<-data.frame(length_Lbra$ID,length_Lbra$Total.Count)
flLbra$length_Lbra.ID<-gsub("\\.1","", flLbra$length_Lbra.ID)
flbra<-as.numeric(as.vector(flLbra$length_Lbra.Total.Count))
names(flbra)<-flLbra$length_Lbra.ID


# make vector of all genes
agLbra<- data.frame(flLbra$length_Lbra.ID)
avLbra<-c(t(agLbra))

# make vector of DE genes or subset of interest genes
#list of genes
Lbra_AMAvsMETA<- read.delim("~/Escritorio/info_paperncRNA/LbraM2903/module_comparison/META_module_comparison.list", 
                            header = F)
Lbra_AMAvsMETA$V1<-gsub("\\.1","", Lbra_AMAvsMETA$V1)
vLbra_AMAvsMETA<-c(t(Lbra_AMAvsMETA))


#Construct a new vector that adds a 0 next to every gene that is not in our DEG list and a 1 next to every gene that is in our DEG list.
gv.METAvsAMALbra=as.integer(avLbra%in%vLbra_AMAvsMETA)
names(gv.METAvsAMALbra)=avLbra 

pwf.LbraMETAvsAMA=nullp(gv.METAvsAMALbra, bias.data = flbra)
head(pwf.LbraMETAvsAMA)

GO.pwf.LbraAMAvsMETA<- goseq(pwf.LbraMETAvsAMA, gene2cat = final_Lbra, use_genes_without_cat=F)
head(GO.pwf.LbraAMAvsMETA)
GO.pwf.LbraAMAvsMETA$padjust<-p.adjust(GO.pwf.LbraAMAvsMETA$over_represented_pvalue, method = "fdr",  n = length(GO.pwf.LbraAMAvsMETA$over_represented_pvalue))

enriched_GO_LbraAMAvsMETA=GO.pwf.LbraAMAvsMETA$category[GO.pwf.LbraAMAvsMETA$padjust<0.05]



#obtain a file with all information about your enriched terms

SubsetGO<-subset(GO.pwf.LbraAMAvsMETA,  padjust<0.05)
write.table(as.data.frame(SubsetGO), 
            file="~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/GO_enrichment_PRO_modules.csv",
            sep = "\t",row.names = F)



capture.output(for(go in enriched_GO_LbraAMAvsMETA) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/GO_enrichment_PRO_modules.txt")
  

g_LbraAMAvsMETA<- SubsetGO%>% 
  top_n(20, wt=-padjust) %>% 
  mutate(hitsPerc=numDEInCat/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  #  expand_limits(x=0) +
  labs(x="Gene ratio", y="GO term", colour="p value", size="Count")
g_LbraAMAvsMETA
#########LdonBPK282###########
gaf_Ldon<-read.delim("~/Escritorio/WGCNA/GAF/TriTrypDB-35_LdonovaniBPK282A1_GO.gaf",
                     header = F,sep = '\t')
final_Ldon<-data.frame(gaf_Ldon$V2, gaf_Ldon$V5)
final_Ldon<-final_Ldon[-1,]

length_Ldon= read.delim("~/Escritorio/WGCNA/CDS_GCcont/TriTrypDB-35_LdonovaniBPK282A1_AnnotatedCDSs.fasta.txt",
                        header = T, sep = "\t")
flLdon<-data.frame(length_Ldon$ID,length_Ldon$Total.Count)
fldon<-as.numeric(as.vector(flLdon$length_Ldon.Total.Count))
names(fldon)<-flLdon$length_Ldon.ID


agLdon<- data.frame(flLdon$length_Ldon.ID)
avLdon<-c(t(agLdon))



mbLdon<- read.delim("~/Escritorio/info_paperncRNA/LdonBPK282/networks/CytoscapeInput-edges-all_modules.list", 
                    header = F)
mbLdon$V1<-gsub("\\.1\\.1","\\.1", mbLdon$V1)
mbvLdon<-c(t(mbLdon))

gv.mbLdon=as.integer(avLdon%in%mbvLdon)
names(gv.mbLdon)=avLdon 

pwf.mbLdon=nullp(gv.mbLdon, bias.data = fldon)
head(pwf.mbLdon)

GO.pwf.mbLdon<- goseq(pwf.mbLdon, gene2cat = final_Ldon, use_genes_without_cat=F)
head(GO.pwf.mbLdon)
GO.pwf.mbLdon$padjust<-p.adjust(GO.pwf.mbLdon$over_represented_pvalue, method = "fdr",  n = length(GO.pwf.mbLdon$over_represented_pvalue))

enriched_GO_Ldon=GO.pwf.mbLdon$category[GO.pwf.mbLdon$padjust<0.05]


SubsetGO<-subset(GO.pwf.mbLdon, padjust <0.05)
write.table(as.data.frame(SubsetGO), 
            file="~/Escritorio/info_paperncRNA/LdonBPK282/GOenrinchment/GO_enrichment_all_modules.csv",
            sep = "\t",row.names = F)

capture.output(for(go in enriched_GO_Ldon) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="~/Escritorio/info_paperncRNA/LdonBPK282/GOenrinchment/GO_enrichment_all_modules.txt")



g_mbLdon<- GO.pwf.mbLdon%>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=category, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  #  expand_limits(x=0) +
  labs(x="Gene ratio", y="GO category", colour="p value", size="Count")



#########LmajFriedlin##########

gaf_ann<- read.delim("~/Escritorio/WGCNA/GAF/TriTrypDB-35_LmajorFriedlin_GO.gaf", 
                     header = F, sep="\t")

final_ann<-data.frame(gaf_ann$V2, gaf_ann$V5)
final_ann<-final_ann[-1,]

# read Gene Length file and convert to numeric vector
gene_length= read.delim("~/Escritorio/WGCNA/CDS_GCcont/TriTrypDB-35_LmajorFriedlin_AnnotatedCDSs.fasta.txt",
                        header = T, sep = "\t")
final_length<-data.frame(gene_length$ID,gene_length$Total.Count)
#final_length$gene_length.ID<-gsub("\\.1","", final_length$gene_length.ID)
fl<-as.numeric(as.vector(final_length$gene_length.Total.Count))
names(fl)<-final_length$gene_length.ID


# make vector of all genes
assayed_genes<- data.frame(final_length$gene_length.ID)
assayed_vector<-c(t(assayed_genes))


############## META LmajFriedlin############
DE_genes<- read.delim("~/Escritorio/info_paperncRNA/LmajFriedlin/networks/CytoscapeInput-edges-allmodules.list", 
                      header = F)
DEG_vector<-c(t(DE_genes))
#Construct a new vector that adds a 0 next to every gene that is not in our DEG list and a 1 next to every gene that is in our DEG list.
gene.vector=as.integer(assayed_vector%in%DEG_vector)
names(gene.vector)=assayed_vector 

pwf=nullp(gene.vector, bias.data = fl)
head(pwf)

go_pwf<- goseq(pwf, gene2cat = final_ann, use_genes_without_cat=F)
head(go_pwf)
go_pwf$padjust<-p.adjust(go_pwf$over_represented_pvalue, method = "fdr",  n = length(go_pwf$over_represented_pvalue))

enriched_GO=go_pwf$category[go_pwf$padjust<0.05]




SubsetGO<-subset(go_pwf, padjust <0.05)
write.table(as.data.frame(SubsetGO), 
            file="~/Escritorio/info_paperncRNA/LmajFriedlin/GOenrichment/GO_enrichment_all_modules2.csv",
            sep = "\t",row.names = F)



capture.output(for(go in enriched_GO) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="~/Escritorio/info_paperncRNA/LmajFriedlin/GOenrichment/GO_enrichment_all_modules.txt")



g_mbLmaj<-go_pwf %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=category, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  #  expand_limits(x=0) +
  labs(x="Gene ratio", y="GO category", colour="p value", size="Count")

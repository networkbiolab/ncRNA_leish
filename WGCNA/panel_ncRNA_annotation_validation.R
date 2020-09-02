library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)

classes<-read.table("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/ncRNAs_Leishmania_spp/clustering_nc/class_Leishmania_filtered", header=F,
                    sep = "\t",
                    as.is=TRUE)
colnames(classes)<-c( "Species", "Class")

new_clasess <-classes %>%
  group_by(Species,Class) %>%
  summarise(n=n())%>%
  mutate(percent = (n / sum(n)), cumsum = cumsum(percent), label=ifelse(Class=="Unclassified", paste0("N=", sum(n)),""))

jerusalem = c("#CBC5B8", "#4F779F", "#729E93", "#99985D", "#699CB6", "#485178",
              "#A47443", "#FAE3A0", "#CEAE88")

theme_set(
  theme_pubr() +
    theme(legend.position = "right")
)


freq_plot<- ggplot(new_clasess,aes(x=Species, y=percent, fill=forcats::fct_rev(Class))) +
  scale_y_continuous(labels = scales::percent) +
  labs ( y = "Percentage", fill= "ncRNA annotation") +
  geom_bar(position = 'fill',color = "black", stat="identity") +
  geom_text(aes(y=cumsum, label=label),hjust = -0.1, size=3,
            inherit.aes = TRUE) +
  guides(fill=guide_legend(reverse=T)) +
  scale_fill_manual(values = jerusalem) +
  coord_flip()



validation<- read.table("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/ncRNAs_Leishmania_spp/clustering_nc/validation_data", header=F,
                        sep = "\t",
                        as.is=TRUE)

colnames(validation)<-c( "Species", "total","w_reads","porcentage_1","CPM_cutoff", "porcentage_2")
validation$no_CPM=validation$total- validation$CPM_cutoff
validation$porcentage_2=(validation$CPM_cutoff *100)/validation$total
validation$porcentage_3=(validation$no_CPM *100)/validation$total
validation$sum = validation$porcentage_2 + validation$porcentage_3

FinVal<- data.frame(validation$Species, validation$porcentage_2, validation$porcentage_3)
colnames(FinVal)<-c("Specie", "TE", "non_TE")


mVal<-melt(FinVal)
val_plot<- ggplot(mVal, aes(x=Specie, y=value, fill=forcats::fct_rev(variable)))+
  scale_y_continuous(labels = scales::percent) +
  labs ( y = "Percentage", fill= "ncRNA transcriptional evidence") +
  geom_bar(position = 'fill',color = "black", stat="identity") +
    guides(fill=guide_legend(reverse=T)) +
  scale_fill_manual(values = jerusalem)

setEPS()
postscript("~/Desktop/Leishmania_info_06-01-20/Imagenes_ncRNAs_paper/panel_Annotation_Validation.eps", width = 20.0, height = 17.0)
ggarrange(freq_plot, val_plot, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()
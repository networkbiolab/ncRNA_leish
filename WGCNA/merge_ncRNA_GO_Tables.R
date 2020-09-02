df1<-read.delim("~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/ncRNA_GOenrichment/tabla_blue_module.tab", 
                header = T)
df2<-read.delim("~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/ncRNA_GOenrichment/tabla_yellow_module.tab", 
                header = T)
df3<-read.delim("~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/ncRNA_GOenrichment/tabla_turquoise_module.tab", 
                header = T)

df4<- merge(x=df1, y=df2, by=colnames(df1),all=T)
df5<- merge(x=df3, y=df4, by=colnames(df1),all=T)
df5<-df5[order(df5$ncRNA_id),]


write.table(df5, 
            file = "~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/ncRNA_GOenrichment/ncRNA_GO_annotation.tsv",
            sep = "\t",
            row.names = F)

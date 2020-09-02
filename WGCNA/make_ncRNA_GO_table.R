df1<-read.delim("~/Escritorio/info_paperncRNA/LbraM2903/module_comparison/genes2.tmp", 
                header = F)
df2<-read.delim("~/Escritorio/info_paperncRNA/LbraM2903/module_comparison/genes3.tmp", 
                header = F)
df3<-merge(x= df1 ,y=df2 , by = "V1", all = TRUE) 
colnames(df3)[4]<- "category"
df4<-read.delim("~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/ModuleComparison/GO_enrichment_PRO_modules.csv",
                header = T)
df5<-merge(x= df3 ,y=df4 , by = "category", all = TRUE) 
colnames(df5)[2]<- "CodGene"
colnames(df5)[3]<-"ncRNA_id"
colnames(df5)[4]<-"ncRNA_Class"
colnames(df5)[5]<-"pvalue"
df5<-na.omit(df5)
df5<-df5[order(df5$ncRNA_id),]
df5 <- df5[,c(1,2,3,4,5,9,10,11,12)]
df5<-df5[c(3,4,1,5,8,6,7,9,2)]
df5<-aggregate(df5[9], df5[-9], 
          FUN = function(X) paste(unique(X), collapse=", "))
df5<-df5[order(df5$ncRNA_id),]
write.table(df5, 
            file = "~/Escritorio/info_paperncRNA/LbraM2903/GOenrichment/ModuleComparison/tabla_PRO_sharedmodules.csv",
            sep = "\t",
            row.names = F)
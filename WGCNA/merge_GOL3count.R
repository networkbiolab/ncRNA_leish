df1<-read.delim("~/Escritorio/info_paperncRNA/LdonBPK282/GOenrinchment/goL3.list", 
                header = T)
df2<-read.delim("~/Escritorio/info_paperncRNA/LdonBPK282/GOenrinchment/count_go.list", 
                header = T)

df4<- merge(x=df1, y=df2, by="category",all=T)

write.table(df4, 
            file = "~/Escritorio/info_paperncRNA/LdonBPK282/GOenrinchment/goL3count.csv",
            sep = "\t",
            row.names = F)


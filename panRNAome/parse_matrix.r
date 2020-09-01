#!/usr/bin/env Rscript

#usage Rscript	--vanilla parse_matrix.r input_matrix.file  out_folder
# input file count matrix 


args = commandArgs(trailingOnly=TRUE)


df <- read.delim(args[1],header = T,row.names ="Cluster")

df$X.conserved<- NULL
df$Classes<-NULL
df2<-data.matrix(df)
df2[df2 > 0] = 1

name="binary_matrix.list"

write.table(df2,
            file= paste(args[2],name, sep= "/"), 
            sep="\t",
            row.names = F, 
            col.names = T) 
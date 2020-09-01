args = commandArgs(trailingOnly=TRUE)



# get a list of files
my.file.list <- list.files(pattern = "htseq", path =args[1],full.names = T)
  # for each file, run read.table and select only the first column
my.list <- lapply(X = my.file.list, FUN = function(x){read.table(x, sep = "\t", header=F)[,2]})


# merge columns that are in a list into one data.frame
my_row<- read.table(my.file.list[1], sep = '\t', header = F)
my.df <- do.call("cbind", my.list)

# extract row names and put name of file as column name 
rownames(my.df)<-(my_row$V1)
path=my.file.list
strsplit(path, "/")->a
b<-sapply(a, tail, 1)
c<-strsplit(b, "[.]")
header <- do.call( rbind, c)[,1]
colnames(my.df)<- header

# save file as table
write.table(my.df, sep = "\t", row.names = T ,col.names = T, file = args[2])



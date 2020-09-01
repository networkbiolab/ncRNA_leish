library(edgeR)

cts <- as.matrix(read.csv("/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/LmajFriedlin_htseq_data/LmajFreidlin_Tritryp_final2.htseq",
                          row.names="gene_id", 
                          header = T,
                          sep = "\t"))
nCts<-cpm(cts, log = T)

write.table(as.data.frame(nCts), 
            file="/media/lalo/Seagate Backup + BK/Leishmania_info_06-01-20/01_htseq_outs/LmajFriedlin_htseq_data/LmajFriedlin_tritryp_CPMLog2.csv",
            sep = "\t")



####flowerfunction #####
flower_plot <- function(sample, value, start, a, b,  
                        ellipse_col = rgb(212, 128, 28, 150, max = 255), 
                        circle_col = rgb(239, 195, 93,102, max = 255),
                        circle_col2 = rgb(239, 195, 93,190, max = 255),
                        circle_text_cex = 1, labels=labels, labels2=labels) {
  par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type="n")
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                          y = 5 + sin((start + deg * (t - 1.5)) * pi / 180), 
                          col = ellipse_col,
                          border = ellipse_col,
                          a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.9 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.9 * sin((start + deg * (t - 1)) * pi / 180),font = 2,
         value[t]
    )
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = circle_text_cex,font = 2
      )
      
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = circle_text_cex,font = 2
      )
    }           
  })
  plotrix::draw.circle(x = 5, y = 5, r = 1.5, col = circle_col, border =circle_col)
  text(x =5, y = 4, labels=labels2,cex = .8,font = 2)
  plotrix::draw.circle(x = 5, y = 5, r = 1, col = circle_col2, border = circle_col2)
  text(x = 5, y = 5, labels=labels,font = 2)
}



#####plot_your_flower########

df<-read.delim("~/Escritorio/PsalminisPangenome/fasta_2/pangenome_Psalmonis/total_counts2.csv", header = T, sep = '\t')
df1<-subset(df, df$strain!="Accessory" & df$strain!="Core")
df1<-df1[order(df1$New_genogroup),]
df2<-subset(df, df$strain=="Accessory" | df$strain=="Core")
df2$New_genogroup=NULL

#ca<-read.delim("~/Escritorio/datasets_tablas_paperncRNA/CoreAccessoryncRNA.csv", header = T)
strian<-as.vector(df1$strain)
count<-as.vector(df1$count)
flower_plot(sample =strian ,
            value = count,
            90, 0.5, 2, 
            labels=paste("Core\n",df2$count[2]),
            labels2 = paste("\n\n\nAccessory",df2$count[1]),
            #ellipse_col =rgb(152, 90, 113,90,max= 255) ,
            #circle_col2 =rgb(152, 51, 102,max= 255), 
            #circle_col = rgb(152, 90, 113,180,max = 255)
)


#blue
#ellipse_col =rgb(69,131,249,90,max= 255) ,
#circle_col2 =rgb(69,131,249,max= 255), 
#circle_col = rgb(179, 193, 220,180,max = 255)
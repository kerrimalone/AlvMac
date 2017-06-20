###############################
# Load required packages     #
##############################
library("ggplot2")
library("reshape2")
library("RColorBrewer")
library("grid")
library("gridExtra")

###############################
# Read in and manipulate data #
##############################

setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/Common DE/")
my_data<-read.csv("common_DE.csv", sep=",", header=TRUE)
head(my_data)
names(my_data)

dim(my_data)
names(my_data)
my_data_bov<-my_data[,c(1,4,6,10,14)]
my_data_bov<-as.data.frame(my_data_bov)
head(my_data_bov)
dim(my_data_bov)
row.names(my_data_bov)<-my_data_bov[,1]
my_data_bov<-my_data_bov[,2:5]
matrix_bov<-as.matrix(my_data_bov)
matrix_bov
longData_bov <- melt(matrix_bov)
head(longData_bov, 10)
melt_lab<-c("Gene","Timepoint.log2FC", "Log2FC")
colnames(longData_bov)<-melt_lab
head(longData_bov, 10)
dim(longData_bov)


names(my_data)
my_data_bov_FDR<-my_data[,c(1,5,7,11,15)]
my_data_bov_FDR<-as.data.frame(my_data_bov_FDR)
names(my_data_bov_FDR)
dim(my_data_bov_FDR)
row.names(my_data_bov_FDR)<-my_data_bov_FDR[,1]
my_data_bov_FDR<-my_data_bov_FDR[,2:5]
matrix_bov_FDR<-as.matrix(my_data_bov_FDR)
head(matrix_bov_FDR,10)
longData_bov_FDR <- melt(matrix_bov_FDR)
head(longData_bov_FDR, 10)
melt_lab<-c("Gene","Timepoint.FDR", "FDR")
colnames(longData_bov_FDR)<-melt_lab
dim(longData_bov)
dim(longData_bov_FDR)
head(longData_bov)
head(longData_bov_FDR)

all_bov<-cbind(longData_bov,longData_bov_FDR)
head(all_bov)
all_bov<-all_bov[,2:6]
all_bov$bool <- (all_bov$FDR < 0.05)
View(all_bov)

dim(my_data)
names(my_data)
my_data_TB<-my_data[,c(1,2,8,12,16)]
my_data_TB<-as.data.frame(my_data_TB)
head(my_data_TB)
dim(my_data_TB)
row.names(my_data_TB)<-my_data_TB[,1]
my_data_TB<-my_data_TB[,2:5]
matrix_TB<-as.matrix(my_data_TB)
head(matrix_TB)
longData_TB <- melt(matrix_TB)
head(longData_TB, 10)
melt_lab<-c("Gene","Timepoint.log2FC", "Log2FC")
colnames(longData_TB)<-melt_lab
head(longData_TB, 10)
names(my_data)

my_data_TB_FDR<-my_data[,c(1,3,9,13,17)]
my_data_TB_FDR<-as.data.frame(my_data_TB_FDR)
dim(my_data_TB_FDR)
row.names(my_data_TB_FDR)<-my_data_TB_FDR[,1]
my_data_TB_FDR<-my_data_TB_FDR[,2:5]
matrix_TB_FDR<-as.matrix(my_data_TB_FDR)
head(matrix_TB_FDR)
longData_TB_FDR <- melt(matrix_TB_FDR)
head(longData_TB_FDR, 10)
melt_lab<-c("Gene","Timepoint.FDR", "FDR")
colnames(longData_TB_FDR)<-melt_lab
dim(longData_TB)
dim(longData_TB_FDR)
head(longData_TB)
head(longData_TB_FDR)

all_TB<-cbind(longData_TB,longData_TB_FDR)
head(all_TB)
all_TB<-all_TB[,2:6]
all_TB$bool <- (all_TB$FDR < 0.05)
View(all_TB)

#########
# Plot #
#########
col1 = colorRampPalette(c("darkslategray","cadetblue4","powderblue"))(4) #changing the number in the bracket alters the gradient
col2 <- rep("snow", 1) #can add in diff sections on gradient
col3 = colorRampPalette(c("orange1","chocolate3","brown","brown4"))(4)
colors2 <- c(col1,col2, col3)

plot_bov <- ggplot(all_bov,
                   aes(x = Timepoint.log2FC, y = (Gene), fill =Log2FC))

plot_bov<- plot_bov + theme(panel.border=element_rect(fill = NA, colour= 'black',size=10))
plot_bov <- plot_bov + geom_tile(colour="black") +
  ylim(rev(levels(all_bov$Gene))) + xlim(levels(all_bov$Timepoint.log2FC))

plot_bov<-plot_bov + geom_point(data=all_bov[all_bov$bool,], aes(x=Timepoint.log2FC, 
                                                                 y=Gene, size= as.numeric(bool)),colour="black",shape=42,size=2)  
plot_bov <- plot_bov + scale_fill_gradientn(colours = colors2,limits=c(-8,8)) 
plot_bov <- plot_bov + coord_fixed()
plot_bov <- plot_bov + theme_bw(base_size = 8, base_family = "")
plot_bov <- plot_bov + theme(text = element_text(size=6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=0.5,face="bold"))
plot_bov <- plot_bov + theme(plot.margin=unit(c(0.5, 0.5, 0.5 , 0.5), "cm")) 
plot_bov <- plot_bov + theme(axis.text.y = element_text(size=4,face="bold.italic",colour="black"))
plot_bov <- plot_bov + theme(axis.text.x = element_text(size=5,face="bold",colour = "black"))
plot_bov <- plot_bov + xlab("") + ylab("")
plot_bov <- plot_bov + theme(legend.key = element_rect(colour="black", size=8)) 
plot_bov <- plot_bov + guides(fill=FALSE)
plot_bov <- plot_bov + scale_x_discrete("", labels = c("logFC_MB_2hr" = "2hr","logFC_MB_6hr" = "6hr","logFC_MB_24hr" = "24hr",
                                                       "logFC_MB_48hr" = "48hr"))

print(plot_bov)

plot_TB <- ggplot(all_TB,
                  aes(x = Timepoint.log2FC, y = (Gene), fill =Log2FC))

plot_TB<- plot_TB + theme(panel.border=element_rect(fill = NA, colour= 'black',size=10))
plot_TB <- plot_TB + geom_tile(colour="black") +
  ylim(rev(levels(all_TB$Gene))) + xlim(levels(all_TB$Timepoint.log2FC))

plot_TB<-plot_TB + geom_point(data=all_TB[all_TB$bool,], aes(x=Timepoint.log2FC, 
                                                             y=Gene, size= as.numeric(bool)),colour="black",shape=42,size=2)  
plot_TB <- plot_TB + scale_fill_gradientn(colours = colors2,limits=c(-8,8)) 
plot_TB <- plot_TB + coord_fixed()
plot_TB <- plot_TB + theme_bw(base_size = 8, base_family = "")
plot_TB <- plot_TB + theme(text = element_text(size=6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=0.5,face="bold"))
plot_TB <- plot_TB + theme(plot.margin=unit(c(0.5, 0.5, 0.5 , -18.25), "cm")) 
plot_TB <- plot_TB + theme(axis.text.x = element_text(size=5,face="bold",colour = "black"))
plot_TB <- plot_TB + theme(axis.text.y=element_blank())
plot_TB <- plot_TB + xlab("") + ylab("")
plot_TB <- plot_TB + theme(legend.key = element_rect(colour="black", size=8)) 
plot_TB <- plot_TB + theme(legend.key = element_rect(colour="black", size=10)) 
leg=expression(paste(Log[2],"",FC))
plot_TB<- plot_TB + labs(fill=leg) 
plot_TB<- plot_TB + theme(legend.title=element_text(size=7)) 
plot_TB<- plot_TB + theme(legend.text=element_text(size=6)) 
plot_TB <- plot_TB + scale_x_discrete("", labels = c("logFC_TB_2hr" = "2hr","logFC_TB_6hr" = "6hr","logFC_TB_24hr" = "24hr",
                                                     "logFC_TB_48hr" = "48hr"))
print(plot_TB)

# Add custom lines and labels
annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
  {
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
  }

Text_star = textGrob("* FDR < 0.05", gp = gpar(fontsize=5.75))
Text1 = textGrob("M. tuberculosis", gp = gpar(fontsize=5,fontface=4))
Text2 = textGrob("M. bovis", gp = gpar(fontsize=5,fontface=4))
Text3 = textGrob("Chemokines", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text4 = textGrob("Interleukins", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text5 = textGrob("CD molecules", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text6 = textGrob("Leukocyte \nimmunoglobulin-\nlike receptors", rot=90, gp=gpar(fontsize=5.5,fontface="bold"))
Text7 = textGrob("TLRs", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text8 = textGrob("Immuno-\nregulatory \nmolecules", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text9 = textGrob("Immuno-\nregulatory \nreceptors", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text10 = textGrob("Immune related genes",gp=gpar(fontsize=8,fontface="bold"))


p1 <- plot_TB + annotation_custom2(Text1,  xmin = 3.5, xmax = 3.5, ymin = -3.5, ymax = -3.5, data = longData_TB) +
  annotation_custom2(linesGrob(), xmin = 0.5, xmax = 4.7, ymin = -2.7, ymax = -2.7, data = longData_TB)  +
  annotation_custom2(Text_star, xmin = 10, xmax = 10, ymin = 27.5, ymax = 27.5, data = longData_TB) 

plot_TB_anno <- ggplot_gtable(ggplot_build(p1))
plot_TB_anno$layout[grepl("panel", plot_TB_anno$layout$name), ]$clip <- "off"
grid.newpage()
grid.draw(plot_TB_anno)


p2 <- plot_bov + annotation_custom2(Text2,  xmin = 2.5, xmax = 2.5, ymin = -3.5, ymax = -3.5, data = longData_bov) +
  annotation_custom2(linesGrob(), xmin = 0.5, xmax = 4.7, ymin = -2.7, ymax = -2.7, data = longData_bov)  +
  annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 60.75, ymax = 77.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 60.75, ymax = 60.75, data = longData_bov)    +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 77.5, ymax = 77.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 46.75, ymax = 60.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 46.75, ymax = 46.75, data = longData_bov)    +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 60.5, ymax = 60.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 32.75, ymax = 46.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 32.75, ymax = 32.75, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 46.5, ymax = 46.5, data = longData_bov) +  
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 24.3, ymax = 32.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 24.3, ymax = 24.3, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 32.5, ymax = 32.5, data = longData_bov) + 
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 20.8, ymax = 24.05, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 24.05, ymax = 24.05, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 20.8, ymax = 20.8, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 11.75, ymax = 20.55, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 11.75, ymax = 11.75, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 20.55, ymax = 20.55, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -6.3, ymin = 0.5, ymax = 11.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 0.5, ymax = 0.5, data = longData_bov)  +
   annotation_custom2(linesGrob(), xmin = -6.3, xmax = -5, ymin = 11.5, ymax = 11.5, data = longData_bov)  +
   annotation_custom2(Text3,  xmin = -7.2, xmax = -7.2, ymin = 68.5, ymax = 68.5, data = longData_bov) + 
   annotation_custom2(Text4,  xmin = -7.2, xmax = -7.2, ymin = 54, ymax = 54, data = longData_bov) +
   annotation_custom2(Text5,  xmin = -7.2, xmax = -7.2, ymin = 39.5, ymax = 39.5, data = longData_bov) +
   annotation_custom2(Text6,  xmin = -8.75, xmax = -8.75, ymin = 28.75, ymax = 28.75, data = longData_bov) +
   annotation_custom2(Text7,  xmin = -7.2, xmax = -7.2, ymin = 22.5, ymax = 22.5, data = longData_bov) +
   annotation_custom2(Text8,  xmin = -8.75, xmax = -8.75, ymin = 16, ymax = 16, data = longData_bov) +
   annotation_custom2(Text9,  xmin = -8.75, xmax = -8.75, ymin = 5.5, ymax = 5.5, data = longData_bov) +
   annotation_custom2(Text10,  xmin = 5.5, xmax = 5.5, ymin = 79.25, ymax = 79.25, data = longData_bov)  



plot_bov_anno <- ggplot_gtable(ggplot_build(p2))
plot_bov_anno$layout[grepl("panel", plot_bov_anno$layout$name), ]$clip <- "off"
grid.newpage()
grid.draw(plot_bov_anno)

grid.arrange(plot_bov_anno,plot_TB_anno,ncol=2)


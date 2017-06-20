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

setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/Warburg")
my_data<-read.csv("warburg_gene_list.csv", sep=",", header=TRUE)
head(my_data)
names(my_data)

dim(my_data)
names(my_data)
my_data_bov<-my_data[,1:5]
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
names(my_data)

my_data_bov_FDR<-my_data[,c(1,11,12,14,16)]
my_data_bov_FDR<-as.data.frame(my_data_bov_FDR)
dim(my_data_bov_FDR)
row.names(my_data_bov_FDR)<-my_data_bov_FDR[,1]
my_data_bov_FDR<-my_data_bov_FDR[,2:5]
matrix_bov_FDR<-as.matrix(my_data_bov_FDR)
matrix_bov_FDR
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
my_data_TB<-my_data[,c(1,6:9)]
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

my_data_TB_FDR<-my_data[,c(1,10,13,15,17)]
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
col1 = colorRampPalette(c("cadetblue4","powderblue"))(4) #changing the number in the bracket alters the gradient
col2 <- rep("snow", 1) #can add in diff sections on gradient
col3 = colorRampPalette(c("orange1","chocolate3","brown"))(4)
colors2 <- c(col1,col2, col3)

plot_bov <- ggplot(all_bov,
                   aes(x = Timepoint.log2FC, y = (Gene), fill =Log2FC))

plot_bov<- plot_bov + theme(panel.border=element_rect(fill = NA, colour= 'black',size=10))
plot_bov <- plot_bov + geom_tile(colour="black") +
  ylim(rev(levels(all_bov$Gene))) + xlim(levels(all_bov$Timepoint.log2FC))

plot_bov<-plot_bov + geom_point(data=all_bov[all_bov$bool,], aes(x=Timepoint.log2FC, 
                                                                 y=Gene, size= as.numeric(bool)),colour="black",shape=42,size=3)  
plot_bov <- plot_bov + scale_fill_gradientn(colours = colors2,limits=c(-4,4)) 
plot_bov <- plot_bov + coord_fixed()
plot_bov <- plot_bov + theme_bw(base_size = 8, base_family = "")
plot_bov <- plot_bov + theme(text = element_text(size=6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=0.5,face="bold"))
plot_bov <- plot_bov + theme(plot.margin=unit(c(1, 1, 1 , 1), "cm")) 
plot_bov <- plot_bov + theme(axis.text.y = element_text(face="bold.italic",colour="black"))
plot_bov <- plot_bov + theme(axis.text.x = element_text(face="bold",colour = "black"))
plot_bov <- plot_bov + xlab("") + ylab("")
plot_bov <- plot_bov + theme(legend.key = element_rect(colour="black", size=8)) 
plot_bov <- plot_bov + guides(fill=FALSE)
plot_bov <- plot_bov + scale_x_discrete("", labels = c("MB.logfc.2" = "2hr","MB.logfc.6" = "6hr","MB.logfc.24" = "24hr",
                                                       "MB.logfc.48" = "48hr"))
print(plot_bov)

plot_TB <- ggplot(longData_TB,
                  aes(x = Timepoint.log2FC, y = (Gene), fill =Log2FC))
plot_TB<- plot_TB + theme(panel.border=element_rect(fill = NA, colour= 'black',size=10))
plot_TB <- plot_TB + geom_tile(colour="black") +
  ylim(rev(levels(all_TB$Gene))) + xlim(levels(all_TB$Timepoint.log2FC))

plot_TB<-plot_TB + geom_point(data=all_TB[all_TB$bool,], aes(x=Timepoint.log2FC, 
                                                                 y=Gene, size= as.numeric(bool)),colour="black",shape=42,size=3)  
plot_TB <- plot_TB + scale_fill_gradientn(colours = colors2,limits=c(-4,4)) 
plot_TB <- plot_TB + coord_fixed()
plot_TB <- plot_TB + theme_bw(base_size = 8, base_family = "")
plot_TB <- plot_TB + theme(text = element_text(size=6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=0.5,face="bold"))
plot_TB <- plot_TB + theme(plot.margin=unit(c(1, 1, 1 , -13.25), "cm")) 
plot_TB <- plot_TB + theme(axis.text.x = element_text(face="bold",colour = "black"))
plot_TB <- plot_TB + theme(axis.text.y=element_blank())
plot_TB <- plot_TB + xlab("") + ylab("")
plot_TB <- plot_TB + theme(legend.key = element_rect(colour="black", size=10)) 
leg=expression(paste(Log[2],"",FC))
plot_TB<- plot_TB + labs(fill=leg) 
plot_TB<- plot_TB + theme(legend.title=element_text(size=7)) 
plot_TB<- plot_TB + theme(legend.text=element_text(size=6)) 
plot_TB <- plot_TB + scale_x_discrete("", labels = c("TB.logfc.2" = "2hr","TB.logfc.6" = "6hr","TB.logfc.24" = "24hr",
                                                     "TB.logfc.48" = "48hr"))
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
Text1 = textGrob("M. tuberculosis", gp = gpar(fontsize=6,fontface=4))
Text2 = textGrob("M. bovis", gp = gpar(fontsize=6,fontface=4))
Text3 = textGrob(label=expression(bold(paste("HIF1", alpha, " regulation"))), rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text4 = textGrob("Glucose transport & metabolism", rot=90, gp=gpar(fontsize=6,fontface="bold"))
Text5 = textGrob("Monocarboxylate transport \n& V-ATPase", rot=90, gp=gpar(fontsize=5.5,fontface="bold"))
Text6 = textGrob("Glucose Metabolism",gp=gpar(fontsize=10,fontface="bold"))


p1 <- plot_TB + annotation_custom2(Text1,  xmin = 2, xmax = 3, ymin = -1.75, ymax = -1.75, data = longData_TB) +
  annotation_custom2(linesGrob(), xmin = 0.5, xmax = 4.6, ymin = -1.25, ymax = -1.25, data = longData_TB)  +
  annotation_custom2(Text_star, xmin = 7, xmax = 7, ymin = 12, ymax = 12, data = longData_TB) 
plot_TB_anno <- ggplot_gtable(ggplot_build(p1))
plot_TB_anno$layout[grepl("panel", plot_TB_anno$layout$name), ]$clip <- "off"
grid.newpage()
grid.draw(plot_TB_anno)

p2 <- plot_bov + annotation_custom2(Text2,  xmin = 2, xmax = 3, ymin = -1.75, ymax = -1.75, data = longData_bov) +
  annotation_custom2(linesGrob(), xmin = 0.5, xmax = 4.6, ymin = -1.25, ymax = -1.25, data = longData_bov)   +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -2.5, ymin = 21.65, ymax = 34.25, data = longData_bov)   +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -1.75, ymin = 34.25, ymax = 34.25, data = longData_bov)    +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -1.75, ymin = 21.65, ymax = 21.65, data = longData_bov)  +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -2.5, ymin = 9.7, ymax = 21.2, data = longData_bov)  + 
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -1.75, ymin = 9.7, ymax = 9.7, data = longData_bov)    +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -1.75, ymin = 21.2, ymax = 21.2, data = longData_bov)  +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -2.5, ymin = 0.5, ymax = 9.3, data = longData_bov)  +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -1.75, ymin = 0.5, ymax = 0.5, data = longData_bov)    +
  annotation_custom2(linesGrob(), xmin = -2.5, xmax = -1.75, ymin = 9.3, ymax = 9.3, data = longData_bov) + 
  annotation_custom2(Text3,  xmin = -2.9, xmax = -2.9, ymin = 28, ymax = 28, data = longData_bov)  + 
  annotation_custom2(Text4,  xmin = -2.9, xmax = -2.9, ymin = 15, ymax = 15, data = longData_bov) + 
  annotation_custom2(Text5,  xmin = -3.3, xmax = -3.3, ymin = 4.75, ymax = 4.75, data = longData_bov) +
  annotation_custom2(Text6,  xmin = 5.05, xmax = 5.05, ymin = 35.5, ymax = 35.5, data = longData_bov)  
  
plot_bov_anno <- ggplot_gtable(ggplot_build(p2))
plot_bov_anno$layout[grepl("panel", plot_bov_anno$layout$name), ]$clip <- "off"
grid.newpage()
grid.draw(plot_bov_anno)

grid.arrange(plot_bov_anno,plot_TB_anno,ncol=2)


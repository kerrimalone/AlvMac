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

setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/ABC transporters/")
#Read in all data for heatmap
my_data<-read.csv("ABC_transporters.csv", sep=",", header=TRUE)
head(my_data)
names(my_data)

#Subset data for logFC & TB
logFC_data<-my_data[,c(1,20,22)]
head(logFC_data)
logFC_data<-as.data.frame(logFC_data)
head(logFC_data)
dim(logFC_data)
row.names(logFC_data)<-logFC_data[,1]
logFC_data<-logFC_data[,2:3]
logFC_data<-na.omit(logFC_data)

#Create matrix and melt
logFC_matrix<-as.matrix(logFC_data)
logFC_matrix
logFC_longData <- melt(logFC_matrix)
head(logFC_longData, 10)
melt_lab<-c("Gene","Timepoint.log2FC", "Log2FC")
colnames(logFC_longData)<-melt_lab
head(logFC_longData, 10)
dim(logFC_longData)

#Subset data for FDR MB & TB
names(my_data)
FDR_data<-my_data[,c(1,21,23)]
head(FDR_data)
FDR_data<-as.data.frame(FDR_data)
head(FDR_data)
dim(FDR_data)
row.names(FDR_data)<-FDR_data[,1]
FDR_data<-FDR_data[,2:3]
FDR_data<-na.omit(FDR_data)

#Create matrix and melt
FDR_matrix<-as.matrix(FDR_data)
FDR_matrix
FDR_longData <- melt(FDR_matrix)
head(FDR_longData, 10)
melt_lab<-c("Gene","Timepoint.FDR", "FDR")
colnames(FDR_longData)<-melt_lab
head(FDR_longData, 10)
dim(FDR_longData)

#Create df for logFC and FDR 
all_logFC<-cbind(logFC_longData,FDR_longData)
head(all_logFC)
all_logFC<-all_logFC[,2:6]
all_logFC$bool <- (all_logFC$FDR < 0.05) #Create boolean for adding FDR * to heatmap
View(all_logFC)

#Subset data for logFC delta
names(my_data)
delta_logFCdata<-my_data[,c(1,24)]
head(delta_logFCdata)
delta_logFCdata<-as.data.frame(delta_logFCdata)
head(delta_logFCdata)
dim(delta_logFCdata)
row.names(delta_logFCdata)<-delta_logFCdata[,1]
delta_logFCdata<-delta_logFCdata[2]
delta_logFCdata<-na.omit(delta_logFCdata)

#Create matrix and melt
delta_logFCmatrix<-as.matrix(delta_logFCdata)
delta_logFCmatrix
delta_logFClongData <- melt(delta_logFCmatrix)
head(delta_logFClongData, 10)
melt_lab<-c("Gene","Timepoint.log2FC", "Log2FC")
colnames(delta_logFClongData)<-melt_lab
head(delta_logFClongData, 10)
dim(delta_logFClongData)


#Subset data for FDR delta
names(my_data)
delta_FDR_data<-my_data[,c(1,25)]
head(delta_FDR_data)
delta_FDR_data<-as.data.frame(delta_FDR_data)
head(delta_FDR_data)
dim(delta_FDR_data)
row.names(delta_FDR_data)<-delta_FDR_data[,1]

#Create matrix and melt
delta_FDR_data<-delta_FDR_data[2]
delta_FDR_data<-na.omit(delta_FDR_data)
delta_FDR_matrix<-as.matrix(delta_FDR_data)
delta_FDR_matrix
delta_FDR_longData <- melt(delta_FDR_matrix)
head(delta_FDR_longData, 10)
melt_lab<-c("Gene","Timepoint.delta_FDR", "delta_FDR")
colnames(delta_FDR_longData)<-melt_lab
head(delta_FDR_longData, 10)
dim(delta_FDR_longData)

#Create df for delta logFC and FDR 
all_delta<-cbind(delta_logFClongData,delta_FDR_longData)
head(all_delta)
all_delta<-all_delta[,2:6]
all_delta$bool <- (all_delta$delta_FDR < 0.05) #Create boolean for adding FDR * to heatmap
View(all_delta)

#########
# Plot #
#########

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

col1 = colorRampPalette(c("darkslategray","cadetblue4","powderblue"))(4) #changing the number in the bracket alters the gradient, keep ratios same
col2 <- rep("snow", 2) 
col3 = colorRampPalette(c("orange1","chocolate3","brown","brown4"))(4)
colors2 <- c(col1,col2, col3)

mylabels_logFC <- c(expression(paste(bolditalic("M. bovis"))), 
                    expression(paste(bolditalic("   MTB")))) 

plot_logFC <- ggplot(all_logFC,
                     aes(x = Timepoint.log2FC, y = (Gene), fill =Log2FC))
plot_logFC<- plot_logFC + theme(panel.border=element_rect(fill = NA, colour= 'black',size=10))
plot_logFC <- plot_logFC + geom_tile(colour="black") +
  ylim(rev(levels(all_logFC$Gene))) + xlim(levels(all_logFC$Timepoint.log2FC))
plot_logFC<-plot_logFC + geom_point(data=all_logFC[all_logFC$bool,], aes(x=Timepoint.log2FC, 
                                                                         y=Gene, size= as.numeric(bool)),colour="black",shape=42,size=4) #Add FDR stars 
plot_logFC <- plot_logFC + scale_fill_gradientn(colours = colors2,limits=c(-5,5),breaks=c(-6,-4,-2,0,2,4,6), 
                                                labels=c("-6","-4","-2"," 0"," 2"," 4"," 6")) 
plot_logFC <- plot_logFC + coord_fixed()
plot_logFC <- plot_logFC + theme_bw(base_size = 8, base_family = "")
plot_logFC <- plot_logFC + theme(text = element_text(size=6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=0.5,face="bold"))
plot_logFC <- plot_logFC + theme(plot.margin=unit(c(2, 2, 2 , -0.5), "cm")) #Moving for custom placement later on
plot_logFC <- plot_logFC + theme(axis.text.y = element_text(size=7,face="bold.italic",colour="black"))
plot_logFC <- plot_logFC + theme(axis.text.x = element_text(size=8,face="bold.italic",colour = "black"))
plot_logFC <- plot_logFC + xlab("") + ylab("")
plot_logFC <- plot_logFC + theme(legend.key = element_rect(colour="black", size=8)) 
plot_logFC <- plot_logFC + theme(legend.key = element_rect(colour="black", size=10)) 
leg=expression(paste(Log[2],"",FC))
plot_logFC<- plot_logFC + labs(fill=leg) 
plot_logFC<- plot_logFC + theme(legend.title=element_text(size=9)) 
plot_logFC<- plot_logFC + theme(legend.text=element_text(size=7)) 
plot_logFC<- plot_logFC + theme(legend.key.size =  unit(0.6, "cm"))
plot_logFC<- plot_logFC + theme(legend.position=c(0.38, 0.65))
plot_logFC <- plot_logFC + scale_x_discrete("", labels = mylabels_logFC)
legend <- g_legend(plot_logFC)
plot_logFC <- plot_logFC + guides(fill=FALSE)
print(plot_logFC)


col1 = colorRampPalette(c("Red3", 'Lightsalmon'))(8) #changing the number in the bracket alters the gradient
col2 <- rep("snow", 2) #can add in diff sections on gradient
col3 = colorRampPalette(c("skyblue2", "navyblue"))(8)
colors2 <- c(col1,col2, col3)

mylabels_delta <- (expression(paste(bold("Delta"))))


plot_delta <- ggplot(all_delta,
                     aes(x = Timepoint.log2FC, y = (Gene), fill =Log2FC))
plot_delta<- plot_delta + theme(panel.border=element_rect(fill = NA, colour= 'black',size=10))
plot_delta <- plot_delta + geom_tile(colour="black") +
  ylim(rev(levels(all_delta$Gene))) + xlim(levels(all_delta$Timepoint.log2FC))
plot_delta<-plot_delta + geom_point(data=all_delta[all_delta$bool,], aes(x=Timepoint.log2FC, 
                                                                         y=Gene, size= as.numeric(bool)),colour="black",shape=42,size=4) #Add FDR stars  
plot_delta <- plot_delta + scale_fill_gradientn(colours = colors2,limits=c(-3,3), breaks=c(-6,-4,-2,0,2,4,6), 
                                                labels=c(" 6"," 4"," 2"," 0"," 2"," 4"," 6"))
plot_delta <- plot_delta + coord_fixed()
plot_delta <- plot_delta + theme_bw(base_size = 8, base_family = "")
plot_delta <- plot_delta + theme(text = element_text(size=6),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=0.5,face="bold"))
plot_delta <- plot_delta + theme(plot.margin=unit(c(2, 2, 2.4 , -15.8), "cm")) #changing to move graph to left for custom placement later on
plot_delta <- plot_delta + theme(axis.text.y = element_text(size=7,face="bold.italic",colour="black"))
plot_delta <- plot_delta + theme(axis.text.x = element_text(size=8,face="bold",colour = "black"))
plot_delta <- plot_delta + xlab("") + ylab("")
plot_delta <- plot_delta + theme(legend.key = element_rect(colour="black", size=8)) 
plot_delta <- plot_delta + theme(legend.key = element_rect(colour="black", size=10)) 
plot_delta <- plot_delta + theme(axis.text.y=element_blank())
leg=expression(paste(Log[2],"",FC))
plot_delta<- plot_delta + labs(fill=leg) 
plot_delta<- plot_delta + theme(legend.title=element_text(size=9)) 
plot_delta<- plot_delta + theme(legend.text=element_text(size=7)) 
plot_delta<- plot_delta + theme(legend.key.size =  unit(0.6, "cm")) #change legend size
plot_delta <- plot_delta + scale_x_discrete("", labels = mylabels_delta)
plot_delta<- plot_delta + theme(legend.position=c(0.38, 0.4))
legend2 <- g_legend(plot_delta)
plot_delta <- plot_delta + guides(fill=FALSE)
print(plot_delta)

#Plot the basic heatmap without annotation
grid.arrange(plot_logFC,plot_delta,ncol=2)
grid.draw(legend)
grid.draw(legend2)

#Annotating the heatmap
#Custom annotation function to allow placement of text and line Grobs outside of grid.arrange margins
annotation_custom2 <- 
  function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
  {
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax)) }

#Can add as many objects as you need here
Text1 = textGrob("48hr post infection", gp = gpar(fontsize=9,fontface="bold"))
Text2 = textGrob("ABC transporters, cholesterol associated",rot=90,gp=gpar(fontsize=14,fontface="bold"))

#Run the following lines together, else you get an error saying "don't know how to add X to plot"
p1 <- plot_delta +
  annotation_custom2(Text2, xmin = -4.3, xmax = -4.3, ymin = 7, ymax = 7, data = all_delta) 

plot_delta_anno <- ggplot_gtable(ggplot_build(p1))
plot_delta_anno$layout[grepl("panel", plot_delta_anno$layout$name), ]$clip <- "off"
grid.newpage()
grid.draw(plot_delta_anno)

p2 <- plot_logFC + annotation_custom2(Text1,  xmin = 2.3, xmax = 2.3 , ymin = -2.0, ymax = -2.0, data = all_logFC) +
  annotation_custom2(linesGrob(), xmin = 0.5, xmax = 4.1, ymin = -1.6, ymax = -1.6, data = all_logFC) 

plot_logFC_anno <- ggplot_gtable(ggplot_build(p2))
plot_logFC_anno$layout[grepl("panel", plot_logFC_anno$layout$name), ]$clip <- "off"

grid.arrange(plot_logFC_anno,plot_delta_anno,ncol=2)
grid.draw(legend)
grid.draw(legend2)


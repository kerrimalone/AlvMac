###############################
# Load required packages     #
##############################
library("ggplot2")

###############################
# Read in and manipulate data #
##############################
setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")
#List of delta 48hr genes logFC > 1, FDR < 0.05
delta_48hr<-read.csv("FDR_0.05_logFC_delta_48H.txt", sep="\t", header=TRUE) 
#Take out rows with NA as these correspond to ncRNA
delta_48hr<-na.omit(delta_48hr)
dim(delta_48hr)
#Pull out gene symbols only and rename.Make as df for merge ("vlookup")
delta_48hr<-(delta_48hr[,3])
str(delta_48hr)
delta_48hr<-as.data.frame(delta_48hr)
head(delta_48hr)
colnames(delta_48hr)[1]<-("gene.symbol")

#Pulling in all data for all genes across timepoint regardless of DE status
all_genes<-read.csv("Full_DE_alv_mac.txt", sep="\t", header=TRUE)
names(all_genes)
#Take out rows with NA as these correspond to ncRNA
all_genes<-na.omit(all_genes)
head(all_genes)
all_genes<-as.data.frame(all_genes)
all_genes<-na.omit(all_genes)
#Pull out columns corresponding to delta model
all_genes_delta<-all_genes[,c(3,5,9,30,34,45,49,60,64)]
names(all_genes_delta)
head(all_genes_delta)
dim(all_genes_delta)

#Pull out records for all genes in delta_48hr from all_genes_delta record,
#using gene.symbol as lookup value
DEgene_merge<-merge(delta_48hr,all_genes_delta, by="gene.symbol")
head(DEgene_merge)
dim(DEgene_merge)
dim(delta_48hr)
head(delta_48hr)

#inserting new rows into df so that I can manually melt the df for graphing
new_col_2hr<-rep("2hr",nrow(DEgene_merge))
new_col_6hr<-rep("6hr",nrow(DEgene_merge))
new_col_24hr<-rep("24hr",nrow(DEgene_merge))
new_col_48hr<-rep("48hr",nrow(DEgene_merge))
length(new_col_6hr)

#inserting the new columns into df
DEgene_merge <- data.frame(DEgene_merge[,1], new_col_2hr, DEgene_merge[,2:length(DEgene_merge)])
head(DEgene_merge)
DEgene_merge <- data.frame(DEgene_merge[,1:4], new_col_6hr, DEgene_merge[,5:length(DEgene_merge)])
head(DEgene_merge)
DEgene_merge <- data.frame(DEgene_merge[,1:7], new_col_24hr, DEgene_merge[,8:length(DEgene_merge)])
head(DEgene_merge)
DEgene_merge <- data.frame(DEgene_merge[,1:10], new_col_48hr, DEgene_merge[,11:length(DEgene_merge)])
head(DEgene_merge)
names(DEgene_merge)[1]<-"gene.symbol"
names(DEgene_merge)


#Pull out columns corresponding to each time point including the new inputted columns
DE_gene_merge_2hr<-DEgene_merge[,c(1,2:4)]
DE_gene_merge_6hr<-DEgene_merge[,c(1,5:7)]
DE_gene_merge_24hr<-DEgene_merge[,c(1,8:10)]
DE_gene_merge_48hr<-DEgene_merge[,c(1,11:13)]
#column names need to be identical for binding to work
names(DE_gene_merge_2hr)<-c("gene.symbol","time","logFC","FDR")
names(DE_gene_merge_6hr)<-c("gene.symbol","time","logFC","FDR")
names(DE_gene_merge_24hr)<-c("gene.symbol","time","logFC","FDR")
names(DE_gene_merge_48hr)<-c("gene.symbol","time","logFC","FDR")

#Bind the sub-column dfs so that the original dataframe is now melted with
#corresponding time value
DEgene_melt<-rbind(DE_gene_merge_2hr,DE_gene_merge_6hr,DE_gene_merge_24hr,
                   DE_gene_merge_48hr)

#########
# Plot #
#########
col1 = colorRampPalette(c("Red3", 'Lightsalmon'))() #changing the number in the bracket alters the gradient
col2 <- rep("white","snow", 1) #can add in diff sections on gradient
col3 = colorRampPalette(c("skyblue2", "navyblue"))(8)
colors2 <- c(col1,col2, col3)
min<-min(DEgene_melt$logFC)
max<-max(DEgene_melt$logFC)
breaks<-c(min[1],-2.5,-1,-0.8,0,0.8,1,3,max[1])
y_lab<-expression("Log"[2]*"FC") #need * to separate sections, [ ]=subscript

p<-ggplot(data=DEgene_melt, aes(x=time, y=logFC, group=gene.symbol, colour=logFC),
          ylim(-6:6)) +
  geom_line() +
  geom_point(aes(fill=logFC))+
  guides(fill=FALSE) +
  #theme(axis.title.y=element_blank()) +
  theme_bw() +
  theme(axis.title.x = element_text(size=10)) +
  theme(axis.title.y = element_text(size=10)) +
  scale_y_continuous(limits = c(-5, 5),breaks=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                     labels=c("5","4","3","2","1","0","1","2","3","4","5"),name=y_lab) +
  scale_x_discrete("Time post-infection") +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=-1,linetype="dotted") +
  geom_hline(yintercept=1,linetype="dotted") +
  geom_point(shape=18,size=0.005) +
  scale_colour_gradientn(breaks = as.vector(breaks), colours = 
                         c("brown4","darkred","firebrick3","palevioletred3","snow1","steelblue1","dodgerblue","blue1","darkblue"), values = as.vector(breaks), 
                       oob = identity, rescaler = function(x,...) x) +
  theme(legend.position="none") +
  annotate("text", x=0.65, y=1.25, label="M. bovis",fontface="italic",size=3) +
  annotate("text", x=0.8, y=-1.25, label="M. tuberculosis", fontface="italic",size=3) +
  annotate("text", x=0.7, y=-5, label="*FDR < 0.05", size=2.75) +
  annotate("text", x=4.025, y=5, label="* ", size=5) +
  annotate("text", x=4.025, y=-3.7, label="* ", size=5) +
  theme(axis.line = element_line(colour = "black",size=0.2),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  
p




########################################
# Put y axis on right of graph
# Rotate graph 90 to the left so that it looks like a Y
# Put a straight line in at y=0
# Replace y axis labels so that the values are positive in both directions for a TB/bovis split
#Colour genes based on logFC using my old red-blue colour ramp

          
# # axis tweaks http://stackoverflow.com/questions/15334494/how-to-change-positions-of-x-and-y-axis-in-ggplot2
# g <- ggplot_gtable(ggplot_build(p))
# ia <- which(g$layout$name == "axis-l")
# ax <- g$grobs[[ia]]$children[[2]]
# ax$widths <- rev(ax$widths)
# ax$grobs <- rev(ax$grobs)
# ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
# pp <- c(subset(g$layout, name == "panel", select = t:r))
# g <- gtable_add_cols(g, g$widths[g$layout[ia, ]$l], length(g$widths) - 1)
# g <-  gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
# g$grobs[[ia]]$children[[2]] <- NULL
# ##############################
# ia <- which(g$layout$name == "ylab")
# ylab <- g$grobs[[ia]]
# g <- gtable_add_cols(g, g$widths[g$layout[ia, ]$l], length(g$widths) - 1)
# g <-  gtable_add_grob(g, ylab, pp$t, length(g$widths) - 1, pp$b)
# g$grobs[[ia]]$label = ''
# grid.draw(g)
# g + theme(axis.title.x = element_text(vjust=-1))



# ##########EXAMPLE
# # Basic line graph with points
# dat1 <- data.frame(
#   sex = factor(c("Female","Female","Male","Male")),
#   time = factor(c("Lunch","Dinner","Lunch","Dinner"), levels=c("Lunch","Dinner")),
#   total_bill = c(13.53, 16.81, 16.24, 17.42)
# )
# dat1
# ggplot(data=dat1, aes(x=time, y=total_bill, group=sex)) +
#   geom_line() +
#   geom_point()
# 
# # Map sex to color
# ggplot(data=dat1, aes(x=time, y=total_bill, group=sex, colour=sex)) +
#   geom_line() +
#   geom_point()
# 
# 
# test.df<-data.frame(
#   sex = factor(c("Female","Female","Male","Male")),
#   time = factor(c("Lunch","Dinner","Lunch","Dinner"), levels=c("Lunch","Dinner")),
#   total_bill = c(13.53, 16.81, 16.24, 17.42)
# )
# 
# test2.df<-data.frame(
#   sex = factor(c("Female"))
# )
# 
# 
# test.merge<-(merge(test2.df,test.df, by="sex"))



# #######################################################################
# basic plot of what I needed
# 
# setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")
# #List of delta 48hr genes logFC > 1, FDR < 0.05
# delta_48hr<-read.csv("FDR_0.05_logFC_delta_48H.txt", sep="\t", header=TRUE) 
# #Take out rows with NA as these correspond to ncRNA
# delta_48hr<-na.omit(delta_48hr)
# dim(delta_48hr)
# #Pull out gene symbols only and rename.Make as df for merge ("vlookup")
# delta_48hr<-(delta_48hr[,3])
# str(delta_48hr)
# delta_48hr<-as.data.frame(delta_48hr)
# head(delta_48hr)
# colnames(delta_48hr)[1]<-("gene.symbol")
# 
# #Pulling in all data for all genes across timepoint regardless of DE status
# all_genes<-read.csv("Full_DE_alv_mac.txt", sep="\t", header=TRUE)
# names(all_genes)
# #Take out rows with NA as these correspond to ncRNA
# all_genes<-na.omit(all_genes)
# head(all_genes)
# all_genes<-as.data.frame(all_genes)
# all_genes<-na.omit(all_genes)
# #Pull out columns corresponding to delta model
# all_genes_delta<-all_genes[,c(3,5,9,30,34,45,49,60,64)]
# names(all_genes_delta)
# colnames(all_genes_delta)[2]<-("logFC_delta_2hr")
# colnames(all_genes_delta)[3]<-("FDR_delta_2hr")
# head(all_genes_delta)
# dim(all_genes_delta)
# 
# #Pull out records for all genes in delta_48hr from all_genes_delta record,
# #using gene.symbol as lookup value
# DEgene_merge<-merge(delta_48hr,all_genes_delta, by="gene.symbol")
# head(DEgene_merge)
# dim(DEgene_merge)
# dim(delta_48hr)
# head(delta_48hr)
# 
# #inserting new rows into df so that I can manually melt the df for graphing
# new_col_2hr<-rep("2hr",nrow(DEgene_merge))
# new_col_6hr<-rep("6hr",nrow(DEgene_merge))
# new_col_24hr<-rep("24hr",nrow(DEgene_merge))
# new_col_48hr<-rep("48hr",nrow(DEgene_merge))
# length(new_col_6hr)
# 
# #inserting the new columns into df
# DEgene_merge <- data.frame(DEgene_merge[,1], new_col_2hr, DEgene_merge[,2:length(DEgene_merge)])
# head(DEgene_merge)
# DEgene_merge <- data.frame(DEgene_merge[,1:4], new_col_6hr, DEgene_merge[,5:length(DEgene_merge)])
# head(DEgene_merge)
# DEgene_merge <- data.frame(DEgene_merge[,1:7], new_col_24hr, DEgene_merge[,8:length(DEgene_merge)])
# head(DEgene_merge)
# DEgene_merge <- data.frame(DEgene_merge[,1:10], new_col_48hr, DEgene_merge[,11:length(DEgene_merge)])
# head(DEgene_merge)
# names(DEgene_merge)[1]<-"gene.symbol"
# names(DEgene_merge)
# 
# 
# #Pull out columns corresponding to each time point including the new inputted columns
# DE_gene_merge_2hr<-DEgene_merge[,c(1,2:4)]
# DE_gene_merge_6hr<-DEgene_merge[,c(1,5:7)]
# DE_gene_merge_24hr<-DEgene_merge[,c(1,8:10)]
# DE_gene_merge_48hr<-DEgene_merge[,c(1,11:13)]
# #column names need to be identical for binding to work
# names(DE_gene_merge_2hr)<-c("gene.symbol","time","logFC","FDR")
# names(DE_gene_merge_6hr)<-c("gene.symbol","time","logFC","FDR")
# names(DE_gene_merge_24hr)<-c("gene.symbol","time","logFC","FDR")
# names(DE_gene_merge_48hr)<-c("gene.symbol","time","logFC","FDR")
# 
# #Bind the sub-column dfs so that the original dataframe is now melted with
# #corresponding time value
# DEgene_melt<-rbind(DE_gene_merge_2hr,DE_gene_merge_6hr,DE_gene_merge_24hr,
#                    DE_gene_merge_48hr)
# 
# 
# p<-ggplot(data=DEgene_melt, aes(x=time, y=logFC, group=gene.symbol, colour=logFC)) +
#   geom_line() +
#   geom_point()
# 
# g <- ggplot_gtable(ggplot_build(p))
# 
# 
# ########################################################

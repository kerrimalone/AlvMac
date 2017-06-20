###############################
# Load required packages     #
##############################

library("ggplot2")
library("plyr")
#Pulling out common DE genes based on logFC FDR lists for each time point.
# Also minusing the delta gene lists from the common DE so that they can be appropriately graphed

#http://stackoverflow.com/questions/38268741/geom-bar-ggplot2-stacked-grouped-bar-plot-with-positive-and-negative-values-p


#Subtract elements in one smaller vector from another using setdiff. Note small vector goes second.
# x<-as.vector(c("1","2","3"))
# y<-as.vector(c("1","2","3","4","5"))
# setdiff(y,x)

###############################
# Read in and manipulate data #
##############################
common_DE_genes_counts.vec<-as.vector(rep("x",8))

setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")
TB_2hr<-read.csv("FDR_0.05_logFC_DE_TB_2H.txt",sep="\t",header=TRUE)
head(TB_2hr)
dim(TB_2hr)
TB_2hr<-na.omit(TB_2hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_2hr)

bovis_2hr<-read.csv("FDR_0.05_logFC_DE_MB_2H.txt",sep="\t",header=TRUE)
head(bovis_2hr)
dim(bovis_2hr)
bovis_2hr<-na.omit(bovis_2hr)
dim(bovis_2hr)

TB_2hr_up<-TB_2hr[TB_2hr$logFC > 1,]
TB_2hr_down<-TB_2hr[TB_2hr$logFC < 1,]
bovis_2hr_up<-bovis_2hr[bovis_2hr$logFC > 1,]
bovis_2hr_down<-bovis_2hr[bovis_2hr$logFC < 1,]

delta_2hr<-read.csv("delta_2hr.txt",sep="\t",header=TRUE)
delta_2hr<-subset(delta_2hr, abs(logFC) > 1)
head(delta_2hr)
#write.csv(delta_2hr, file="logFC_only_delta_2hr.csv")
dim(delta_2hr)
delta_2hr<-na.omit(delta_2hr)
dim(delta_2hr)
delta_2hr_genes<-as.vector(delta_2hr[,"gene.symbol"])

bovis_2hr_up_genes<-as.vector(bovis_2hr_up[,"gene.symbol"])
bovis_2hr_down_genes<-as.vector(bovis_2hr_down[,"gene.symbol"])
TB_2hr_up_genes<-as.vector(TB_2hr_up[,"gene.symbol"])
TB_2hr_down_genes<-as.vector(TB_2hr_down[,"gene.symbol"])

#Taking out delta genes from common DE gene lists
length(bovis_2hr_up_genes)
bovis_2hr_up_genes<-setdiff(bovis_2hr_up_genes,delta_2hr_genes)
length(bovis_2hr_up_genes)
length(bovis_2hr_down_genes)
bovis_2hr_down_genes<-setdiff(bovis_2hr_down_genes,delta_2hr_genes)
length(bovis_2hr_down_genes)
length(TB_2hr_up_genes)
TB_2hr_up_genes<-setdiff(TB_2hr_up_genes,delta_2hr_genes)
length(TB_2hr_up_genes)
length(TB_2hr_down_genes)
TB_2hr_down_genes<-setdiff(TB_2hr_down_genes,delta_2hr_genes)
length(TB_2hr_down_genes)


intersect_2hr_up.vec<-intersect(x =TB_2hr_up_genes, y=bovis_2hr_up_genes)
intersect_2hr_down.vec<-intersect(x =TB_2hr_down_genes, y=bovis_2hr_down_genes)


common_DE_genes_counts.vec[1]<-length(intersect_2hr_up.vec)
common_DE_genes_counts.vec[2]<-length(intersect_2hr_down.vec)
common_DE_genes_counts.vec

TB_6hr<-read.csv("FDR_0.05_logFC_DE_TB_6H.txt",sep="\t",header=TRUE)
head(TB_6hr)
dim(TB_6hr)
TB_6hr<-na.omit(TB_6hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_6hr)

bovis_6hr<-read.csv("FDR_0.05_logFC_DE_MB_6H.txt",sep="\t",header=TRUE)
head(bovis_6hr)
dim(bovis_6hr)
bovis_6hr<-na.omit(bovis_6hr)
dim(bovis_6hr)

delta_6hr<-read.csv("delta_6hr.txt",sep="\t",header=TRUE)
delta_6hr<-subset(delta_6hr, abs(logFC) > 1)
head(delta_6hr)
#write.csv(delta_6hr, file="logFC_only_delta_6hr.csv")
dim(delta_6hr)
delta_6hr<-na.omit(delta_6hr)
dim(delta_6hr)
delta_6hr_genes<-as.vector(delta_6hr[,"gene.symbol"])

TB_6hr_up<-TB_6hr[TB_6hr$logFC > 1,]
TB_6hr_down<-TB_6hr[TB_6hr$logFC < 1,]
bovis_6hr_up<-bovis_6hr[bovis_6hr$logFC > 1,]
bovis_6hr_down<-bovis_6hr[bovis_6hr$logFC < 1,]

bovis_6hr_up_genes<-as.vector(bovis_6hr_up[,"gene.symbol"])
bovis_6hr_down_genes<-as.vector(bovis_6hr_down[,"gene.symbol"])
TB_6hr_up_genes<-as.vector(TB_6hr_up[,"gene.symbol"])
TB_6hr_down_genes<-as.vector(TB_6hr_down[,"gene.symbol"])

length(bovis_6hr_up_genes)
bovis_6hr_up_genes<-setdiff(bovis_6hr_up_genes,delta_6hr_genes)
length(bovis_6hr_up_genes)
length(bovis_6hr_down_genes)
bovis_6hr_down_genes<-setdiff(bovis_6hr_down_genes,delta_6hr_genes)
length(bovis_6hr_down_genes)
length(TB_6hr_up_genes)
TB_6hr_up_genes<-setdiff(TB_6hr_up_genes,delta_6hr_genes)
length(TB_6hr_up_genes)
length(TB_6hr_down_genes)
TB_6hr_down_genes<-setdiff(TB_6hr_down_genes,delta_6hr_genes)
length(TB_6hr_down_genes)


intersect_6hr_up.vec<-intersect(x =TB_6hr_up_genes, y=bovis_6hr_up_genes)
intersect_6hr_down.vec<-intersect(x =TB_6hr_down_genes, y=bovis_6hr_down_genes)
length(intersect_6hr_up.vec)
length(intersect_6hr_down.vec)

common_DE_genes_counts.vec[3]<-length(intersect_6hr_up.vec)
common_DE_genes_counts.vec[4]<-length(intersect_6hr_down.vec)
common_DE_genes_counts.vec

TB_24hr<-read.csv("FDR_0.05_logFC_DE_TB_24H.txt",sep="\t",header=TRUE)
head(TB_24hr)
dim(TB_24hr)
TB_24hr<-na.omit(TB_24hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_24hr)

bovis_24hr<-read.csv("FDR_0.05_logFC_DE_MB_24H.txt",sep="\t",header=TRUE)
head(bovis_24hr)
dim(bovis_24hr)
bovis_24hr<-na.omit(bovis_24hr)
dim(bovis_24hr)

delta_24hr<-read.csv("delta_24hr.txt",sep="\t",header=TRUE)
delta_24hr<-subset(delta_24hr, abs(logFC) > 1)
head(delta_24hr)
#write.csv(delta_24hr, file="logFC_only_delta_24hr.csv")
dim(delta_24hr)
delta_24hr<-na.omit(delta_24hr)
dim(delta_24hr)
delta_24hr_genes<-as.vector(delta_24hr[,"gene.symbol"])

TB_24hr_up<-TB_24hr[TB_24hr$logFC > 1,]
TB_24hr_down<-TB_24hr[TB_24hr$logFC < 1,]
bovis_24hr_up<-bovis_24hr[bovis_24hr$logFC > 1,]
bovis_24hr_down<-bovis_24hr[bovis_24hr$logFC < 1,]

bovis_24hr_up_genes<-as.vector(bovis_24hr_up[,"gene.symbol"])
bovis_24hr_down_genes<-as.vector(bovis_24hr_down[,"gene.symbol"])
TB_24hr_up_genes<-as.vector(TB_24hr_up[,"gene.symbol"])
TB_24hr_down_genes<-as.vector(TB_24hr_down[,"gene.symbol"])

length(bovis_24hr_up_genes)
bovis_24hr_up_genes<-setdiff(bovis_24hr_up_genes,delta_24hr_genes)
length(bovis_24hr_up_genes)
length(bovis_24hr_down_genes)
bovis_24hr_down_genes<-setdiff(bovis_24hr_down_genes,delta_24hr_genes)
length(bovis_24hr_down_genes)
length(TB_24hr_up_genes)
TB_24hr_up_genes<-setdiff(TB_24hr_up_genes,delta_24hr_genes)
length(TB_24hr_up_genes)
length(TB_24hr_down_genes)
TB_24hr_down_genes<-setdiff(TB_24hr_down_genes,delta_24hr_genes)
length(TB_24hr_down_genes)

intersect_24hr_up.vec<-intersect(x =TB_24hr_up_genes, y=bovis_24hr_up_genes)
intersect_24hr_down.vec<-intersect(x =TB_24hr_down_genes, y=bovis_24hr_down_genes)
length(intersect_24hr_up.vec)
length(intersect_24hr_down.vec)

common_DE_genes_counts.vec[5]<-length(intersect_24hr_up.vec)
common_DE_genes_counts.vec[6]<-length(intersect_24hr_down.vec)
common_DE_genes_counts.vec

TB_48hr<-read.csv("FDR_0.05_logFC_DE_TB_48H.txt",sep="\t",header=TRUE)
head(TB_48hr)
dim(TB_48hr)
TB_48hr<-na.omit(TB_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_48hr)

bovis_48hr<-read.csv("FDR_0.05_logFC_DE_MB_48H.txt",sep="\t",header=TRUE)
head(bovis_48hr)
dim(bovis_48hr)
bovis_48hr<-na.omit(bovis_48hr)
dim(bovis_48hr)

delta_48hr<-read.csv("FDR_0.05_logFC_delta_48H.txt",sep="\t",header=TRUE)
delta_48hr<-subset(delta_48hr, abs(logFC) > 1)
head(delta_48hr)
#write.csv(delta_48hr, file="logFC_only_delta_48hr.csv")
dim(delta_48hr)
delta_48hr<-na.omit(delta_48hr)
dim(delta_48hr)
delta_48hr_genes<-as.vector(delta_48hr[,"gene.symbol"])

TB_48hr_up<-TB_48hr[TB_48hr$logFC > 1,]
TB_48hr_down<-TB_48hr[TB_48hr$logFC < 1,]
bovis_48hr_up<-bovis_48hr[bovis_48hr$logFC > 1,]
bovis_48hr_down<-bovis_48hr[bovis_48hr$logFC < 1,]

bovis_48hr_up_genes<-as.vector(bovis_48hr_up[,"gene.symbol"])
bovis_48hr_down_genes<-as.vector(bovis_48hr_down[,"gene.symbol"])
TB_48hr_up_genes<-as.vector(TB_48hr_up[,"gene.symbol"])
TB_48hr_down_genes<-as.vector(TB_48hr_down[,"gene.symbol"])

length(bovis_48hr_up_genes)
bovis_48hr_up_genes<-setdiff(bovis_48hr_up_genes,delta_48hr_genes)
length(bovis_48hr_up_genes)
length(bovis_48hr_down_genes)
bovis_48hr_down_genes<-setdiff(bovis_48hr_down_genes,delta_48hr_genes)
length(bovis_48hr_down_genes)
length(TB_48hr_up_genes)
TB_48hr_up_genes<-setdiff(TB_48hr_up_genes,delta_48hr_genes)
length(TB_48hr_up_genes)
length(TB_48hr_down_genes)
TB_48hr_down_genes<-setdiff(TB_48hr_down_genes,delta_48hr_genes)
length(TB_48hr_down_genes)


intersect_48hr_up.vec<-intersect(x =TB_48hr_up_genes, y=bovis_48hr_up_genes)
intersect_48hr_down.vec<-intersect(x =TB_48hr_down_genes, y=bovis_48hr_down_genes)
length(intersect_48hr_up.vec)
length(intersect_48hr_down.vec)


common_DE_genes_counts.vec[7]<-length(intersect_48hr_up.vec)
common_DE_genes_counts.vec[8]<-length(intersect_48hr_down.vec)
common_DE_genes_counts.vec

Time.vec<-c("02hr","02hr","06hr","06hr","24hr","24hr","48hr","48hr")
Variable.vec<-c("Up","Down","Up","Down","Up","Down","Up","Down")

bar_data.raw<-data.frame(a=character(),b=character(),d=numeric())
bar_data<-rbind(bar_data.raw, data.frame(a=Time.vec, b=Variable.vec, c=as.numeric(common_DE_genes_counts.vec)))
colnames(bar_data)<-c("Time","Variable","Value")

#########
# Plot #
########
label_1<-expression(paste("Up"))
label_2<-expression(paste("Down")) 

q<-ggplot(bar_data, aes(Time), ylim(-400:500)) + 
  geom_bar(data = subset(bar_data, Variable == "Up"), 
           aes(y = Value, fill = Variable), stat = "identity", position = "dodge",width=0.85,colour="black",size=0.4) +
  scale_fill_manual(values=c("#53aaba","#c97e2a"), 
                    name=" ") +
  geom_bar(data = subset(bar_data, Variable == "Down"), #colours are bovis up, tb up, bovis down, tb down
           aes(y = -Value, fill = Variable), stat = "identity", position = "dodge",width=0.85, colour="black",size=0.4) + 
  geom_hline(yintercept = 0,colour = "black") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text.align = 0)  #aligning the legend labels to legend boxes


q + 
  geom_text(data = subset(bar_data, Variable == "Up"), 
            aes(Time, Value, group=Time, label=Value),
            position = position_dodge(width=0.9), vjust = -0.25, size=4) +
  geom_text(data = subset(bar_data, Variable == "Down"), 
            aes(Time, -Value, group=Time, label=Value),
            position = position_dodge(width=0.9), vjust = 1.25, size=4) +
  coord_cartesian(ylim = c(-400, 500)) +
  scale_x_discrete(name="Time post-infection", breaks=c("02hr","06hr","24hr","48hr"),
                   labels=c("2hr","6hr","24hr","48hr")) + #getting rid of the 0 in 02hr and 06hr
  scale_y_continuous("Number of differentially expressed genes") +
  theme(legend.text=element_text(size=9),legend.key.size=unit(0.4,"cm")) + #changing size of legend
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.title.y=element_text(size=11)) +
  theme(aspect.ratio = 1.3) +
  theme(legend.position="bottom", legend.box = "horizontal") + #horizontal legend at bottom of graph
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.2))


common_genes_all<-list(intersect_2hr_up.vec,intersect_2hr_down.vec,intersect_6hr_up.vec,
                       intersect_6hr_down.vec,intersect_24hr_up.vec,intersect_24hr_down.vec,
                       intersect_48hr_up.vec,intersect_48hr_down.vec)
names(common_genes_all)<-c("2hr_up","2hr_down","6hr_up","6hr_down","24hr_up","24hr_down","48hr_up","48hr_down")
names(common_genes_all)

#Write out screen print
sink("all_common_DE_genes_FDR_logFC_final.txt")
print(common_genes_all)
sink() #closes sink file
###############################
# Load required packages     #
##############################

library("ggplot2")

#http://stackoverflow.com/questions/38268741/geom-bar-ggplot2-stacked-grouped-bar-plot-with-positive-and-negative-values-p

###############################
# Read in and manipulate data #
##############################
# Set working directory and load any previously saved data
setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")

#Create vectors with desirable variables for graphing
Time.vec<-c(rep("02hr",4),rep("06hr",4),rep("24hr",4),rep("48hr",4))
Condition.vec<-c("MB","TB","MB","TB","MB","TB","MB","TB","MB","TB","MB","TB","MB","TB","MB","TB")
Variable.vec<-c(rep("Up",2),rep("Down",2),rep("Up",2),rep("Down",2),rep("Up",2),rep("Down",2),
                rep("Up",2),rep("Down",2))
Variable.condition.vec<-c("MB_up","TB_up","MB_down","TB_down","MB_up","TB_up","MB_down","TB_down",
                          "MB_up","TB_up","MB_down","TB_down","MB_up","TB_up","MB_down","TB_down")

#Going to count how many genes are up and down so,
#set up empty vector with blank entries to store the gene counts
values.vec<-rep("x",16)

#Read in DE gene data for each timepoint and treatment and subset
#FDR < 0.05 and log2FC > 1
TB_2hr<-read.csv("FDR_0.05_logFC_DE_TB_2H.txt",sep="\t",header=TRUE)
head(TB_2hr)
dim(TB_2hr)
TB_2hr<-na.omit(TB_2hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_2hr)
TB_2hr_logFC<-as.vector(TB_2hr["logFC"])
#counting up and down regged genes based on logFC values 
#save result in a particular entry in the blank vector values.vec
values.vec[2] <-sum(TB_2hr_logFC > 1)
values.vec[4] <-sum(TB_2hr_logFC < 1)

#Repeat for all timepoints and treatments
bovis_2hr<-read.csv("FDR_0.05_logFC_DE_MB_2H.txt",sep="\t",header=TRUE)
head(bovis_2hr)
dim(bovis_2hr)
bovis_2hr<-na.omit(bovis_2hr)
dim(bovis_2hr)
bovis_2hr_logFC<-as.vector(bovis_2hr["logFC"])
values.vec[1] <-sum(bovis_2hr_logFC > 1)
values.vec[3] <-sum(bovis_2hr_logFC < 1)

bovis_6hr<-read.csv("FDR_0.05_logFC_DE_MB_6H.txt",sep="\t",header=TRUE)
head(bovis_6hr)
dim(bovis_6hr)
bovis_6hr<-na.omit(bovis_6hr)
dim(bovis_6hr)
bovis_6hr_logFC<-as.vector(bovis_6hr["logFC"])
values.vec[5] <-sum(bovis_6hr_logFC > 1)
values.vec[7] <-sum(bovis_6hr_logFC < 1)

TB_6hr<-read.csv("FDR_0.05_logFC_DE_TB_6H.txt",sep="\t",header=TRUE)
head(TB_6hr)
dim(TB_6hr)
TB_6hr<-na.omit(TB_6hr)
dim(TB_6hr)
TB_6hr_logFC<-as.vector(TB_6hr["logFC"])
values.vec[6] <-sum(TB_6hr_logFC > 1)
values.vec[8] <-sum(TB_6hr_logFC < 1)

bovis_24hr<-read.csv("FDR_0.05_logFC_DE_MB_24H.txt",sep="\t",header=TRUE)
head(bovis_24hr)
dim(bovis_24hr)
bovis_24hr<-na.omit(bovis_24hr)
dim(bovis_24hr)
bovis_24hr_logFC<-as.vector(bovis_24hr["logFC"])
values.vec[9] <-sum(bovis_24hr_logFC > 1)
values.vec[11] <-sum(bovis_24hr_logFC < 1)

TB_24hr<-read.csv("FDR_0.05_logFC_DE_TB_24H.txt",sep="\t",header=TRUE)
head(TB_24hr)
dim(TB_24hr)
TB_24hr<-na.omit(TB_24hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_24hr)<-TB_24hr[,1]
dim(TB_24hr)
TB_24hr_logFC<-as.vector(TB_24hr["logFC"])
values.vec[10] <-sum(TB_24hr_logFC > 1)
values.vec[12] <-sum(TB_24hr_logFC < 1)

bovis_48hr<-read.csv("FDR_0.05_logFC_DE_MB_48H.txt",sep="\t",header=TRUE)
head(bovis_48hr)
dim(bovis_48hr)
bovis_48hr<-na.omit(bovis_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(bovis_48hr)
bovis_48hr_logFC<-as.vector(bovis_48hr["logFC"])
values.vec[13] <-sum(bovis_48hr_logFC > 1)
values.vec[15] <-sum(bovis_48hr_logFC < 1)

TB_48hr<-read.csv("FDR_0.05_logFC_DE_TB_48H.txt",sep="\t",header=TRUE)
head(TB_48hr)
dim(TB_48hr)
TB_48hr<-na.omit(TB_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_48hr)
TB_48hr_logFC<-as.vector(TB_48hr["logFC"])
values.vec[14] <-sum(TB_48hr_logFC > 1)
values.vec[16] <-sum(TB_48hr_logFC < 1)
values.vec


#create a new df to store all of the above info with desired variables
bar_data.raw<-data.frame(a=character(),b=character(),c=character(),d=numeric(), e=character())
bar_data<-rbind(bar_data.raw, data.frame(a=Time.vec, b=Condition.vec, c=Variable.vec, d=as.numeric(values.vec), e=Variable.condition.vec))
colnames(bar_data)<-c("Time","Condition","Variable","Value","Variable.condition")

#Make custom labels for legend of graph to include both italicised and plain text
label_1<-expression(paste(italic("M. bovis")," up"))
label_2<-expression(paste(italic("M. tuberculosis")," up")) 
label_3<-expression(paste(italic("M. bovis")," down")) 
label_4<-expression(paste(italic("M. tuberculosis")," down"))


#########
# Plot #
#########
q<-ggplot(bar_data, aes(Time), ylim(-1300:1300)) + 
  geom_bar(data = subset(bar_data, Variable == "Up"), 
           aes(y = Value, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4) +
  scale_fill_manual(values=c("#75a5e5","#323cd3","#f7c0cb","#bc2944","#b5b1b2","#605e5f"), 
                    name=" ",
                    breaks=c("MB_up", "TB_up", "MB_down", "TB_down"), #define the
                    #breaks so that you can relabel
                    labels=c(label_1,label_2,label_3,label_4)) +
  geom_bar(data = subset(bar_data, Variable == "Down"), #colours are bovis up, tb up, bovis down, tb down
           aes(y = -Value, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4) + 
  geom_hline(yintercept = 0,colour = "black") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text.align = 0)  #aligning the legend labels to legend boxes


q + 
  geom_text(data = subset(bar_data, Variable == "Up"), 
            aes(Time, Value, group=Condition, label=Value),
            position = position_dodge(width=0.9), vjust = -0.25, size=4) +
  geom_text(data = subset(bar_data, Variable == "Down"), 
            aes(Time, -Value, group=Condition, label=Value),
            position = position_dodge(width=0.9), vjust = 1.25, size=4) +
  coord_cartesian(ylim = c(-1300, 1300)) +
  scale_x_discrete(name="Time post-infection", breaks=c("02hr","06hr","24hr","48hr"),
                   labels=c("2hr","6hr","24hr","48hr")) + #getting rid of the 0 in 02hr and 06hr
  scale_y_continuous("Number of differentially expressed genes") +
  theme(legend.text=element_text(size=9),legend.key.size=unit(0.4,"cm")) + #changing size of legend
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.title.y=element_text(size=11)) +
  theme(legend.position="bottom", legend.box = "horizontal") + #horizontal legend at bottom of graph
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.2))



#######################################
# Read in and manipulate data for FDR #
#######################################
setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")

#Create vectors with desirable variables for graphing
Time.vec<-c(rep("02hr",4),rep("06hr",4),rep("24hr",4),rep("48hr",4))
Condition.vec<-c("MB","TB","MB","TB","MB","TB","MB","TB","MB","TB","MB","TB","MB","TB","MB","TB")
Variable.vec<-c(rep("Up",2),rep("Down",2),rep("Up",2),rep("Down",2),rep("Up",2),rep("Down",2),
                rep("Up",2),rep("Down",2))
Variable.condition.vec<-c("MB_up","TB_up","MB_down","TB_down","MB_up","TB_up","MB_down","TB_down",
                          "MB_up","TB_up","MB_down","TB_down","MB_up","TB_up","MB_down","TB_down")

#Going to count how many genes are up and down so,
#set up empty vector with blank entries to store the gene counts
values.vec<-rep("x",16)

#Read in DE gene data for each timepoint and treatment and subset
TB_2hr<-read.csv("FDR_0.05_DE_TB_2H.txt",sep="\t",header=TRUE)
head(TB_2hr)
dim(TB_2hr)
TB_2hr<-na.omit(TB_2hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_2hr)
TB_2hr<-as.vector(TB_2hr["logFC"])
#counting up and down regged genes based on FDR values 
#save result in a particular entry in the blank vector values.vec
values.vec[2] <-sum(TB_2hr > 0)
values.vec[4] <-sum(TB_2hr < 0)

#Repeat for all timepoints and treatments
bovis_2hr<-read.csv("FDR_0.05_DE_MB_2H.txt",sep="\t",header=TRUE)
head(bovis_2hr)
dim(bovis_2hr)
bovis_2hr<-na.omit(bovis_2hr)
dim(bovis_2hr)
bovis_2hr<-as.vector(bovis_2hr["logFC"])
values.vec[1] <-sum(bovis_2hr > 0)
values.vec[3] <-sum(bovis_2hr < 0)

bovis_6hr<-read.csv("FDR_0.05_DE_MB_6H.txt",sep="\t",header=TRUE)
head(bovis_6hr)
dim(bovis_6hr)
bovis_6hr<-na.omit(bovis_6hr)
dim(bovis_6hr)
bovis_6hr<-as.vector(bovis_6hr["logFC"])
values.vec[5] <-sum(bovis_6hr > 0)
values.vec[7] <-sum(bovis_6hr < 0)

TB_6hr<-read.csv("FDR_0.05_DE_TB_6H.txt",sep="\t",header=TRUE)
head(TB_6hr)
dim(TB_6hr)
TB_6hr<-na.omit(TB_6hr)
dim(TB_6hr)
TB_6hr<-as.vector(TB_6hr["logFC"])
values.vec[6] <-sum(TB_6hr > 0)
values.vec[8] <-sum(TB_6hr < 0)

bovis_24hr<-read.csv("FDR_0.05_DE_MB_24H.txt",sep="\t",header=TRUE)
head(bovis_24hr)
dim(bovis_24hr)
bovis_24hr<-na.omit(bovis_24hr)
dim(bovis_24hr)
bovis_24hr<-as.vector(bovis_24hr["logFC"])
values.vec[9] <-sum(bovis_24hr > 0)
values.vec[11] <-sum(bovis_24hr < 0)

TB_24hr<-read.csv("FDR_0.05_DE_TB_24H.txt",sep="\t",header=TRUE)
head(TB_24hr)
dim(TB_24hr)
TB_24hr<-na.omit(TB_24hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_24hr)<-TB_24hr[,1]
dim(TB_24hr)
TB_24hr<-as.vector(TB_24hr["logFC"])
values.vec[10] <-sum(TB_24hr > 0)
values.vec[12] <-sum(TB_24hr < 0)

bovis_48hr<-read.csv("FDR_0.05_DE_MB_48H.txt",sep="\t",header=TRUE)
head(bovis_48hr)
dim(bovis_48hr)
bovis_48hr<-na.omit(bovis_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(bovis_48hr)
bovis_48hr<-as.vector(bovis_48hr["logFC"])
values.vec[13] <-sum(bovis_48hr > 0)
values.vec[15] <-sum(bovis_48hr < 0)

TB_48hr<-read.csv("FDR_0.05_DE_TB_48H.txt",sep="\t",header=TRUE)
head(TB_48hr)
dim(TB_48hr)
TB_48hr<-na.omit(TB_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
dim(TB_48hr)
TB_48hr<-as.vector(TB_48hr["logFC"])
values.vec[14] <-sum(TB_48hr > 0)
values.vec[16] <-sum(TB_48hr < 0)
values.vec


#create a new df to store all of the above info with desired variables
bar_data.raw<-data.frame(a=character(),b=character(),c=character(),d=numeric(), e=character())
bar_data<-rbind(bar_data.raw, data.frame(a=Time.vec, b=Condition.vec, c=Variable.vec, d=as.numeric(values.vec), e=Variable.condition.vec))
colnames(bar_data)<-c("Time","Condition","Variable","Value","Variable.condition")

#########
# Plot #
#########
#Make custom labels for legend of graph to include both italicised and plain text
label_1<-expression(paste(italic("M. bovis")," up"))
label_2<-expression(paste(italic("M. tuberculosis")," up")) 
label_3<-expression(paste(italic("M. bovis")," down")) 
label_4<-expression(paste(italic("M. tuberculosis")," down"))

q<-ggplot(bar_data, aes(Time), ylim(-4000:4000)) + 
  geom_bar(data = subset(bar_data, Variable == "Up"), 
           aes(y = Value, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4) +
  scale_fill_manual(values=c("#75a5e5","#323cd3","#f7c0cb","#bc2944","#b5b1b2","#605e5f"), 
                    name=" ",
                    breaks=c("MB_up", "TB_up", "MB_down", "TB_down"), #define the
                    #breaks so that you can relabel
                    labels=c(label_1,label_2,label_3,label_4)) +
  geom_bar(data = subset(bar_data, Variable == "Down"), #colours are bovis up, tb up, bovis down, tb down
           aes(y = -Value, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4) + 
  geom_hline(yintercept = 0,colour = "black") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text.align = 0)  #aligning the legend labels to legend boxes


q + 
  geom_text(data = subset(bar_data, Variable == "Up"), 
            aes(Time, Value, group=Condition, label=Value),
            position = position_dodge(width=0.9), vjust = -0.25, size=4) +
  geom_text(data = subset(bar_data, Variable == "Down"), 
            aes(Time, -Value, group=Condition, label=Value),
            position = position_dodge(width=0.9), vjust = 1.25, size=4) +
  coord_cartesian(ylim = c(-4000, 4000)) +
  scale_x_discrete(name="Time post-infection", breaks=c("02hr","06hr","24hr","48hr"),
                   labels=c("2hr","6hr","24hr","48hr")) + #getting rid of the 0 in 02hr and 06hr
  scale_y_continuous("Number of differentially expressed genes") +
  theme(legend.text=element_text(size=9),legend.key.size=unit(0.4,"cm")) + #changing size of legend
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.title.y=element_text(size=11)) +
  theme(legend.position="bottom", legend.box = "horizontal") + #horizontal legend at bottom of graph
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.2))






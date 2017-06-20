###############################
# Load required packages     #
##############################
library("VennDiagram")     
library("ggplot2")         
library("gridExtra")
library("grid")
library("lattice")

###############################
# Read in and manipulate data #
##############################
setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")

#Read in DE gene data for each timepoint and treatment
#FDR < 0.05 and log2FC > 1
TB_2hr<-read.csv("FDR_0.05_logFC_DE_TB_2H.txt",sep="\t",header=TRUE)
head(TB_2hr)
dim(TB_2hr)
TB_2hr<-na.omit(TB_2hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_2hr)<-TB_2hr[,1]
DE_TB_2hr.vector <- c(rownames(TB_2hr))

bovis_2hr<-read.csv("FDR_0.05_logFC_DE_MB_2H.txt",sep="\t",header=TRUE)
head(bovis_2hr)
dim(bovis_2hr)
bovis_2hr<-na.omit(bovis_2hr)
rownames(bovis_2hr)<-bovis_2hr[,1]
DE_MB_2hr.vector <- c(rownames(bovis_2hr))

bovis_6hr<-read.csv("FDR_0.05_logFC_DE_MB_6H.txt",sep="\t",header=TRUE)
head(bovis_6hr)
dim(bovis_6hr)
bovis_6hr<-na.omit(bovis_6hr)
rownames(bovis_6hr)<-bovis_6hr[,1]
DE_MB_6hr.vector <- c(rownames(bovis_6hr))


TB_6hr<-read.csv("FDR_0.05_logFC_DE_TB_6H.txt",sep="\t",header=TRUE)
head(TB_6hr)
dim(TB_6hr)
TB_6hr<-na.omit(TB_6hr)
rownames(TB_6hr)<-TB_6hr[,1]
DE_TB_6hr.vector <- c(rownames(TB_6hr))

bovis_24hr<-read.csv("FDR_0.05_logFC_DE_MB_24H.txt",sep="\t",header=TRUE)
head(bovis_24hr)
dim(bovis_24hr)
bovis_24hr<-na.omit(bovis_24hr)
rownames(bovis_24hr)<-bovis_24hr[,1]
DE_MB_24hr.vector <- c(rownames(bovis_24hr))


TB_24hr<-read.csv("FDR_0.05_logFC_DE_TB_24H.txt",sep="\t",header=TRUE)
head(TB_24hr)
dim(TB_24hr)
TB_24hr<-na.omit(TB_24hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_24hr)<-TB_24hr[,1]
DE_TB_24hr.vector <- c(rownames(TB_24hr))

bovis_48hr<-read.csv("FDR_0.05_logFC_DE_MB_48H.txt",sep="\t",header=TRUE)
head(bovis_48hr)
dim(bovis_48hr)
bovis_48hr<-na.omit(bovis_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(bovis_48hr)<-bovis_48hr[,1]
DE_MB_48hr.vector <- c(rownames(bovis_48hr))


TB_48hr<-read.csv("FDR_0.05_logFC_DE_TB_48H.txt",sep="\t",header=TRUE)
head(TB_48hr)
dim(TB_48hr)
TB_48hr<-na.omit(TB_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_48hr)<-TB_48hr[,1]
DE_TB_48hr.vector <- c(rownames(TB_48hr))

#########
# Plot #
#########
# Plot all time points per treatment
venn.plot.bovis <- venn.diagram(list(A = DE_MB_2hr.vector,
                                    B = DE_MB_6hr.vector,
                                    C = DE_MB_24hr.vector,
                                    D = DE_MB_48hr.vector),
                               filename        = NULL,
                               fill            = c("#c6e1ef",
                                                   "#a1d7f4",
                                                   "#2aa7ea",
                                                   "#0c15cc"),
                               euler.d=FALSE,
                               scaled=FALSE,
                               sub=substitute( paste(bolditalic('M. bovis'))),
                               sub.fontfamily = "Arial",
                               sub.cex=1.2,
                               sub.pos=c(0.5,1),
                               lwd=1,
                               alpha           = rep(0.50,4),
                               label.col       = "#003333",
                               cex             = 0.85,
                               fontfamily      = "Arial",
                               category.names  = c("2hr",
                                                   "6hr",
                                                   "24hr",
                                                   "48hr"),
                               cat.pos=c(-5,10,10,0),
                               cat.col         = "black",
                               cat.cex         = 0.95,
                               cat.fontfamily  = "Arial",
                               cat.fontface    = 2,
                               rotation.degree = 360,
                               margin          = 0,
                               height          = 10,
                               width           = 4,
                               units           = 'cm',
                               compression     = 'lzw',
                               resolution      = 1200)

venn.plot.TB <- venn.diagram(list(A = DE_TB_2hr.vector,
                                     B = DE_TB_6hr.vector,
                                     C = DE_TB_24hr.vector,
                                     D = DE_TB_48hr.vector),
                                filename        = NULL,
                                fill            = c("#eabbc4",
                                                    "#ed8e9f",
                                                    "#ef3456",
                                                    "#96031d"),
                                euler.d=FALSE,
                                scaled=FALSE,
                                lwd=1,
                                alpha           = rep(0.50,4),
                                sub=substitute( paste(bolditalic('M. tuberculosis'))),
                                sub.fontfamily = "Arial",
                                sub.cex=1.2,
                                sub.pos=c(0.5,1),
                                label.col       = "#003333",
                                cex             = 0.85,
                                fontfamily      = "Arial",
                                category.names  = c("2hr",
                                                    "6hr",
                                                    "24hr",
                                                    "48hr"),
                                cat.pos=c(-5,10,10,0),
                                cat.col         = "black",
                                cat.cex         = 0.9,
                                cat.fontfamily  = "Arial",
                                cat.fontface    = 2,
                                rotation.degree = 360,
                                margin          = 0,
                                height          = 10,
                                width           = 4,
                                units           = 'cm',
                                compression     = 'lzw',
                                resolution      = 1200)


#title1=textGrob(expression("FDR < 0.05, -1 > Log"[2]*"FC > 1"), gp=gpar(fontface="bold", fontsize = 10))
quartz()
venn.plot<-grid.arrange(gTree(children=venn.plot.bovis),
             gTree(children=venn.plot.TB),
             ncol = 2,
             widths = c(1,1),
             heights = c(1,1))
             

#top=title1)



################FDR only################

###############################
# Read in and manipulate data #
##############################

TB_2hr<-read.csv("FDR_0.05_DE_TB_2H.txt",sep="\t",header=TRUE)
head(TB_2hr)
dim(TB_2hr)
TB_2hr<-na.omit(TB_2hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_2hr)<-TB_2hr[,1]
DE_TB_2hr.vector <- c(rownames(TB_2hr))

bovis_2hr<-read.csv("FDR_0.05_DE_MB_2H.txt",sep="\t",header=TRUE)
head(bovis_2hr)
dim(bovis_2hr)
bovis_2hr<-na.omit(bovis_2hr)
rownames(bovis_2hr)<-bovis_2hr[,1]
DE_MB_2hr.vector <- c(rownames(bovis_2hr))

bovis_6hr<-read.csv("FDR_0.05_DE_MB_6H.txt",sep="\t",header=TRUE)
head(bovis_6hr)
dim(bovis_6hr)
bovis_6hr<-na.omit(bovis_6hr)
rownames(bovis_6hr)<-bovis_6hr[,1]
DE_MB_6hr.vector <- c(rownames(bovis_6hr))


TB_6hr<-read.csv("FDR_0.05_DE_TB_6H.txt",sep="\t",header=TRUE)
head(TB_6hr)
dim(TB_6hr)
TB_6hr<-na.omit(TB_6hr)
rownames(TB_6hr)<-TB_6hr[,1]
DE_TB_6hr.vector <- c(rownames(TB_6hr))

bovis_24hr<-read.csv("FDR_0.05_DE_MB_24H.txt",sep="\t",header=TRUE)
head(bovis_24hr)
dim(bovis_24hr)
bovis_24hr<-na.omit(bovis_24hr)
rownames(bovis_24hr)<-bovis_24hr[,1]
DE_MB_24hr.vector <- c(rownames(bovis_24hr))


TB_24hr<-read.csv("FDR_0.05_DE_TB_24H.txt",sep="\t",header=TRUE)
head(TB_24hr)
dim(TB_24hr)
TB_24hr<-na.omit(TB_24hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_24hr)<-TB_24hr[,1]
DE_TB_24hr.vector <- c(rownames(TB_24hr))

bovis_48hr<-read.csv("FDR_0.05_DE_MB_48H.txt",sep="\t",header=TRUE)
head(bovis_48hr)
dim(bovis_48hr)
bovis_48hr<-na.omit(bovis_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(bovis_48hr)<-bovis_48hr[,1]
DE_MB_48hr.vector <- c(rownames(bovis_48hr))


TB_48hr<-read.csv("FDR_0.05_DE_TB_48H.txt",sep="\t",header=TRUE)
head(TB_48hr)
dim(TB_48hr)
TB_48hr<-na.omit(TB_48hr) #NA rows correspond to ncRNAs and thus do not have gene symbols. Need to be removed for venn.
rownames(TB_48hr)<-TB_48hr[,1]
DE_TB_48hr.vector <- c(rownames(TB_48hr))

#########
# Plot #
#########
# Plot all time points per treatment for FDR only 
venn.plot.bovis <- venn.diagram(list(A = DE_MB_2hr.vector,
                                     B = DE_MB_6hr.vector,
                                     C = DE_MB_24hr.vector,
                                     D = DE_MB_48hr.vector),
                                filename        = NULL,
                                fill            = c("#c6e1ef",
                                                    "#a1d7f4",
                                                    "#2aa7ea",
                                                    "#0c15cc"),
                                euler.d=FALSE,
                                scaled=FALSE,
                                sub=substitute( paste(bolditalic('M. bovis'))),
                                sub.fontfamily = "Arial",
                                sub.cex=1.2,
                                sub.pos=c(0.5,1),
                                lwd=1,
                                alpha           = rep(0.50,4),
                                label.col       = "#003333",
                                cex             = 0.85,
                                fontfamily      = "Arial",
                                category.names  = c("2hr",
                                                    "6hr",
                                                    "24hr",
                                                    "48hr"),
                                cat.pos=c(-5,10,10,0),
                                cat.col         = "black",
                                cat.cex         = 0.95,
                                cat.fontfamily  = "Arial",
                                cat.fontface    = 2,
                                rotation.degree = 360,
                                margin          = 0,
                                height          = 10,
                                width           = 4,
                                units           = 'cm',
                                compression     = 'lzw',
                                resolution      = 1200)

venn.plot.TB <- venn.diagram(list(A = DE_TB_2hr.vector,
                                  B = DE_TB_6hr.vector,
                                  C = DE_TB_24hr.vector,
                                  D = DE_TB_48hr.vector),
                             filename        = NULL,
                             fill            = c("#eabbc4",
                                                 "#ed8e9f",
                                                 "#ef3456",
                                                 "#96031d"),
                             euler.d=FALSE,
                             scaled=FALSE,
                             lwd=1,
                             alpha           = rep(0.50,4),
                             sub=substitute( paste(bolditalic('M. tuberculosis'))),
                             sub.fontfamily = "Arial",
                             sub.cex=1.2,
                             sub.pos=c(0.5,1),
                             label.col       = "#003333",
                             cex             = 0.85,
                             fontfamily      = "Arial",
                             category.names  = c("2hr",
                                                 "6hr",
                                                 "24hr",
                                                 "48hr"),
                             cat.pos=c(-5,10,10,0),
                             cat.col         = "black",
                             cat.cex         = 0.9,
                             cat.fontfamily  = "Arial",
                             cat.fontface    = 2,
                             rotation.degree = 360,
                             margin          = 0,
                             height          = 10,
                             width           = 4,
                             units           = 'cm',
                             compression     = 'lzw',
                             resolution      = 1200)


#title1=textGrob(expression("FDR < 0.05, -1 > Log"[2]*"FC > 1"), gp=gpar(fontface="bold", fontsize = 10))
quartz()
venn.plot<-grid.arrange(gTree(children=venn.plot.bovis),
                        gTree(children=venn.plot.TB),
                        ncol = 2,
                        widths = c(1,1),
                        heights = c(1,1))


#############################
# Install required packages # ----
#############################

source("https://bioconductor.org/biocLite.R")
biocLite("AnnotationFuncs")
biocLite("Biobase")
biocLite("edgeR")
biocLite("geneLenDataBase")
biocLite("GO.db")
biocLite("goseq")
biocLite("limma")
biocLite("org.Bt.eg.db")

install.packages("dplyr")
install.packages("extrafont")
install.packages("gdata")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("magrittr")
install.packages("MASS")
install.packages("plyr")
install.packages("svglite")
install.packages("devtools")
install.packages("VennDiagram")
##########################
# Load required packages # ----
##########################

library(AnnotationFuncs)  # Version 1.22.0
library(Biobase)          # Version 2.32.0
library(plyr)             # Version 1.8.4
library(dplyr)            # Version 0.4.3
library(edgeR)            # Version 3.2.4_1
library(extrafont)        # Version 0.17
library(ggplot2)          # Version 2.1.0
library(GO.db)            # Version 3.3.0
library(goseq)            # Version 1.24.0
library(grDevices)        # Version 3.3.0
library(grid)             # Version 3.3.0
library(gridExtra)        # Version 2.2.1
library(limma)            # Version 3.28.5
library(magrittr)         # Version 1.5
library(MASS)             # Version 7.3-45
library(org.Bt.eg.db)     # Version 3.3.0
library(RColorBrewer)     # Version 1.1-2
library(svglite)          # Version 1.1.0, not working?
library(tools)            # Version 3.3.0
library(VennDiagram)      # Version 1.6.17

# Set working directory and load any previously saved data
setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR")
getwd()
workDir <- getwd()
workDir
# Create a vector of all files names
fileDir="/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR/counts"
raw_files <- list.files(path = fileDir,
                    pattern         = "*txt",
                    all.files       = TRUE,
                    full.names      = FALSE,
                    recursive       = FALSE,
                    ignore.case     = FALSE)
raw_files

# Reads and merges a set of files containing counts
Counts <- readDGE(files = raw_files, path=fileDir, header = TRUE,columns = c(1, 7), comment.char = "#")
names(Counts)
head(Counts$samples)
head(Counts$counts)
dim(Counts)

# Output data
write.table(x = Counts$samples,file = "alv_mac_samples.txt", sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)
write.table(x = Counts$counts,file = "alv_mac_raw_counts.txt", sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)

#Clean input files
raw.counts <- read.table(file   = "alv_mac_raw_counts.txt",
                         header = TRUE)
colnames(raw.counts)
head(raw.counts)
#Remove extraneous file extenstion info from column sample names to leave animal,timepoint and treatment info
colnames(raw.counts) <- gsub("_alignment_sense.counts","",colnames(raw.counts))
head(raw.counts)
#Keep ENTREZID numbers only for gene column
head(rownames(raw.counts))
rownames(raw.counts)<-gsub("BGD.*,","",rownames(raw.counts)) #.* means any no. of characters
head(rownames(raw.counts))
rownames(raw.counts)<-gsub("GeneID:","",rownames(raw.counts))
head(rownames(raw.counts))
rownames(raw.counts)<-gsub(",miRBase.*","",rownames(raw.counts))
head(rownames(raw.counts))
#Check outputs and make sure all is okay
head(raw.counts)
rownames(raw.counts)
colnames(raw.counts)
write.table(x = raw.counts,file = "alv_mac_raw_counts_clean.txt", sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)

#Sample file needs to include Animal, time and condition data and file names formatted as above
samples <- read.table(file   = "alv_mac_samples.txt",
                      header = TRUE)
names(samples)
head(samples)
rownames(samples)<-c()
samples$files <- gsub("_alignment_sense-counts.txt","",samples$files)
head(samples)
#Copy the files column three times and append to samples df
samples<-cbind(samples,samples[,1])
samples<-cbind(samples,samples[,1])
samples<-cbind(samples,samples[,1])
dim(samples)
head(samples)
colnames(samples)[5]<-"Animal"
colnames(samples)[6]<-"Time"
colnames(samples)[7]<-"Condition"
head(samples)
samples$Animal <- gsub("_.*","",samples$Animal)
samples$Time <- gsub(".*_.._","",samples$Time)
samples$Condition <- gsub(".*CN.*","CN",samples$Condition)
samples$Condition <- gsub(".*TB.*","TB",samples$Condition)
samples$Condition <- gsub(".*MB.*","MB",samples$Condition)
head(samples)
#Remove group column
samples <- subset(samples, select = -c(group) )
#Reorder columns for tidyness
dim(samples)
samples <-samples[c(1,4,5,6,2,3)]
head(samples)
write.table(x = samples,file = "alv_mac_samples_clean.txt", sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)


##################################################################
# Get gene information using the org.Bt.eg.db annotation package # ----
##################################################################

# Create annotation table with counts information
annotated.counts <- raw.counts
head(annotated.counts)
columns(org.Bt.eg.db)

# Get gene names from NCBI gene identifiers
annotated.counts$gene.name <- mapIds(org.Bt.eg.db,keys = rownames(annotated.counts),column = "GENENAME",keytype = "ENTREZID",multiVals = "first")

# Get gene symbols from NCBI gene identifiers
annotated.counts$gene.symbol <- mapIds(org.Bt.eg.db,
                                       keys      = rownames(annotated.counts),
                                       column    = "SYMBOL",
                                       keytype   = "ENTREZID",
                                       multiVals = "first")

# Get ENSEMBL gene ids from NCBI gene identifiers
annotated.counts$ENSEMBL.tag <- mapIds(org.Bt.eg.db,
                                       keys      = rownames(annotated.counts),
                                       column    = "ENSEMBL",
                                       keytype   = "ENTREZID",
                                       multiVals = "first")

head(annotated.counts)
dim(annotated.counts)

# Output data
write.table(x         = annotated.counts,
            file      = "alv_mac_raw_counts_clean_annotated.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

###########rough below works to a certain point
head(annotated.counts)
gene.annotation <- dplyr::select(annotated.counts,
                                 gene.name,
                                 gene.symbol,
                                 ENSEMBL.tag)
head(gene.annotation)
dim(gene.annotation)
rownames(gene.annotation)

# Create DGElist containing information about condition and annotation
head(raw.counts)
head(samples)
#Check classes of variables
str(samples)
samples$Condition<-as.factor(samples$Condition)
condition <- factor(relevel(samples$Condition, ref = "CN"))
head(condition)

Alv_mac_dgelist <- DGEList(counts       = raw.counts,
                        group        = condition,
                        genes        = gene.annotation,
                        lib.size     = NULL,
                        norm.factors = NULL,
                        remove.zeros = FALSE)

names(Alv_mac_dgelist)
dim(Alv_mac_dgelist)
head(Alv_mac_dgelist$counts)
head(Alv_mac_dgelist$samples)
head(Alv_mac_dgelist$genes)
rownames(Alv_mac_dgelist$samples)

# Add animal, batch and time point information to DGElist
samples.info <- dplyr::select(samples, Animal, Time, Condition)
head(samples.info)
rownames(samples.info)<-rownames(Alv_mac_dgelist$samples)
#dim(samples)
#samples<-samples[,2:6]
#head(samples)
Alv_mac_dgelist$samples <- merge(x  = Alv_mac_dgelist$samples,
                                 y  = samples.info,
                                 by = "row.names") #0 means rownames, Carol had "row.names" but it wouldn't work for me

head(Alv_mac_dgelist$samples)
rownames(samples.info)
rownames(Alv_mac_dgelist$samples)

# Correct row names
rownames(Alv_mac_dgelist$samples) <- Alv_mac_dgelist$samples[, 1]
Alv_mac_dgelist$samples <- Alv_mac_dgelist$samples[, -1]
head(Alv_mac_dgelist$samples)

#####################################################################
# Quality check of libraries by plotting density of raw gene counts # ----
#####################################################################

# Log10 transform the count data for better visualization
count_log10 <- log10(x = (Alv_mac_dgelist$counts[ , 1 : ncol(Alv_mac_dgelist$counts)]
                          + 1))
summary(count_log10)

# Plot density of raw counts for all libraries
png(filename = "Density_raw_alv_mac.png", width = 1366, height = 768, units = "px")

plot(x    = density(count_log10[, 1]),
     main = "Density plot of raw counts per gene",
     lty  = 1,
     xlab = "Log10 of raw counts per gene",
     ylab = "Density",
     col  = "black",
     ylim = c(0.0,1.25))


for (i in 2 : ncol(count_log10)) {
  lines(density(count_log10[, i]),
        lty = 1,
        col = "black")
}

dev.off()
##############################################
# Filtering of zero and lowly expressed tags # ----
##############################################

# Filter non expressed tags (all genes that have zero counts in all samples)
alv_mac_no_zeros <- Alv_mac_dgelist[rowSums(Alv_mac_dgelist$counts) > 1, ]
dim(alv_mac_no_zeros$counts)
head(alv_mac_no_zeros$counts)
rownames(alv_mac_no_zeros$counts)

alv_mac_filt <- alv_mac_no_zeros[rowSums(cpm(alv_mac_no_zeros$counts) > 1) >= 10, ]
head(alv_mac_filt$counts)


# Output file of filtered counts
write.table(x         = alv_mac_filt$counts,
            file      = "alv_mac_filtered_rawcounts.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

# Recompute the library size
colSums(alv_mac_filt$counts)
alv_mac_filt$samples<-Alv_mac_dgelist$samples
head(alv_mac_filt$samples)

alv_mac_filt$samples$lib.size <- base::colSums(alv_mac_filt$counts)
head(alv_mac_filt$samples)
head(Alv_mac_dgelist$samples)

########################################################
# Normalization of data using Trimmed Mean of M-values # ----
#    (based on RNA composition between libraries)      #
########################################################

# Calculate normalisation factor for our DGElist.
# With edgeR, counts are not transformed in any way after normalization,
# instead normalization will modify library size.
alv_mac_norm <- calcNormFactors(alv_mac_filt, method = "TMM")
head(alv_mac_norm$samples) #check norm.factors before and after normalisation

#####################################################################
# Quality check of filtered libraries by plotting density of counts # ----
#####################################################################

# Log10 transform the filtered count data for better visualization
count_filt_log10 <- log10(alv_mac_norm$counts[, 1 : ncol(alv_mac_norm$counts)] + 1)

# Plot density of count for all libraries
png(filename = "Density_alv_mac_post_filter.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plot(density(count_filt_log10[, 1]),
     main = "Density plot of count per gene post filtering",
     lty  = 1,
     xlab = "Log10 of count per gene",
     ylab = "Density",
     col  = "black",
     ylim = c(0.0,0.6))

for (i in 2 : ncol(count_filt_log10)) {
  lines(density(count_filt_log10[, i]),
        lty = 1,
        col = "black")
}

dev.off()

#############################################################
# Exploratory data analysis: Multidimensional scaling plots # ----
#############################################################

# Plot MDS of all samples
png(filename = "MDS_all_samples.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(alv_mac_norm)

dev.off()

# Plot MDS of 0hr time point
png(filename = "MDS_0hr.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = alv_mac_norm[, grep(pattern = "_0H", x = colnames(alv_mac_norm))])

dev.off()

# Plot MDS of 2hr time point
png(filename = "MDS_2hr.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = alv_mac_norm[, grep(pattern = "_2H", x = colnames(alv_mac_norm))])

dev.off()

# Plot MDS of Week 24Hr time point
png(filename = "MDS_24Hr.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = alv_mac_norm[, grep(pattern = "_24H", x = colnames(alv_mac_norm))])

dev.off()

# Plot MDS of Week 24hr time point
png(filename = "MDS_24hr.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = alv_mac_norm[, grep(pattern = "_24H", x = colnames(alv_mac_norm))])

dev.off()

# Plot MDS of Week 48hr time point
png(filename = "MDS_48hr.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(x = alv_mac_norm[, grep(pattern = "_48H", x = colnames(alv_mac_norm))])

dev.off()

###############################
# Define experimental factors # ---- 
###############################

#View(alv_mac_norm$samples)
head(alv_mac_norm$samples$Condition)
head(alv_mac_norm$samples$Time)


######################################
#cond.time interaction model, paired analysis
#######################################
Animal <-factor(alv_mac_norm$samples$Animal)
head(Animal)

Condition<-factor(alv_mac_norm$samples$Condition, levels=c("CN","MB","TB"))
head(Condition)
Condition<-relevel(Condition,ref="CN")

Time<-factor(alv_mac_norm$samples$Time)
head(Time)
Time<-relevel(Time,ref="2H")

Cond.Time.levels<-rownames(alv_mac_norm$samples)
head(Cond.Time.levels)

Cond.Time.levels<-c(rownames(alv_mac_norm$samples))
head(Cond.Time.levels)
as.data.frame(Cond.Time.levels)
Cond.Time.levels<-gsub("N.*_C","C",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_T","T",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_M","M",Cond.Time.levels) 
head(Cond.Time.levels)
Cond.Time.levels<-unique(Cond.Time.levels)
Cond.Time.levels
Cond.Time <- factor(paste(alv_mac_norm$samples$group,
                          alv_mac_norm$samples$Time,
                          sep="_"),
                    levels = c(Cond.Time.levels))
head(Cond.Time)
Cond.Time<-relevel(Cond.Time,ref="CN_2H")

test_matrix_2H <- model.matrix(~Animal + Cond.Time,
                             data = alv_mac_norm$samples)

colnames(test_matrix_2H)
dim(test_matrix_2H)
#View(test_matrix_2H)

write.table(x         = test_matrix_2H,
            file      = "alv_mac_test_matrix_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

############################################################################
# Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
#                      (for multi-factor data)                             #
############################################################################

alv_mac_disp_2H <- estimateGLMCommonDisp(y    = alv_mac_norm,
                                      design  = test_matrix_2H,
                                      verbose = TRUE)

alv_mac_disp_2H <- estimateGLMTrendedDisp(y   = alv_mac_disp_2H,
                                       design = test_matrix_2H)

alv_mac_disp_2H <- estimateGLMTagwiseDisp(y   = alv_mac_disp_2H,
                                       design = test_matrix_2H)

names(alv_mac_disp_2H)

# Plot the dispersion
png(filename = "BCV_alv_mac_test_matrix_2H.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotBCV(alv_mac_disp_2H)

dev.off()

# Show the calculated dispersion
alv_mac_disp_2H$common.dispersion

# And show its square root, the coefficient of biological variation
sqrt(alv_mac_disp_2H$common.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information #http://seqanswers.com/forums/showthread.php?t=5591
Tagwisedisp_2H <- cbind(alv_mac_disp_2H$genes, alv_mac_disp_2H$tagwise.dispersion)
#View(Tagwisedisp_2H)
dim(Tagwisedisp_2H)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp_2H,
             file = "alv_mac_Tagwise_dispersion_test_matrix_2H.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs # ----
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
alv_mac_fit_2H <- glmFit(y = alv_mac_disp_2H, design = test_matrix_2H)
names(alv_mac_fit_2H)
colnames(alv_mac_fit_2H$design)
dim(alv_mac_fit_2H)
# Test for differential expression between the different time points/treatments
MB_2H_lrt <- glmLRT(alv_mac_fit_2H, coef = "Cond.TimeMB_2H")
DE_MB_2H <- topTags(object        = MB_2H_lrt,
                  n             = "inf",
                  adjust.method = "BH")
#View(DE_MB_2H$table)

FDR_0.05_DE_MB_2H <- subset(DE_MB_2H$table, FDR < 0.05)
#View(FDR_0.05_DE_MB_2H)
write.table(x         = FDR_0.05_DE_MB_2H,
            file      = "FDR_0.05_DE_MB_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_MB_2H  <- subset(DE_MB_2H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_MB_2H)
write.table(x         = FDR_0.05_logFC_DE_MB_2H,
            file      = "FDR_0.05_logFC_DE_MB_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

TB_2H_lrt <- glmLRT(alv_mac_fit_2H, coef = "Cond.TimeTB_2H")
DE_TB_2H <- topTags(object        = TB_2H_lrt,
                    n             = "inf",
                    adjust.method = "BH")
#View(DE_TB_2H$table)

FDR_0.05_DE_TB_2H <- subset(DE_TB_2H$table, FDR < 0.05)
#View(FDR_0.05_DE_TB_2H)
write.table(x         = FDR_0.05_DE_TB_2H,
            file      = "FDR_0.05_DE_TB_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_TB_2H  <- subset(DE_TB_2H$table, abs(logFC) > 1 & FDR < 0.05)
View(FDR_0.05_logFC_DE_TB_2H)
write.table(x         = FDR_0.05_logFC_DE_TB_2H,
            file      = "FDR_0.05_logFC_DE_TB_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

colnames(alv_mac_fit_2H$coefficients)
#MB-TB
delta_2H_lrt <- glmLRT(alv_mac_fit_2H, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0))
DE_delta_2H <- topTags(object        = delta_2H_lrt,
                    n             = "inf",
                    adjust.method = "BH")
#View(DE_delta_2H$table)

FDR_0.05_delta_2H <- subset(DE_delta_2H$table, FDR < 0.05)
#View(FDR_0.05_delta_2H)
write.table(x         = FDR_0.05_delta_2H,
            file      = "FDR_0.05_delta_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_delta_2H  <- subset(DE_delta_2H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_delta_2H)
write.table(x         = FDR_0.05_logFC_delta_2H,
            file      = "FDR_0.05_logFC_delta_2H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

##########
# 6H
#########
Animal <-factor(alv_mac_norm$samples$Animal)
head(Animal)

Condition<-factor(alv_mac_norm$samples$Condition, levels=c("CN","MB","TB"))
head(Condition)
Condition<-relevel(Condition,ref="CN")

Time<-factor(alv_mac_norm$samples$Time)
head(Time)
Time<-relevel(Time,ref="6H")

Cond.Time.levels<-rownames(alv_mac_norm$samples)
head(Cond.Time.levels)

Cond.Time.levels<-c(rownames(alv_mac_norm$samples))
head(Cond.Time.levels)
as.data.frame(Cond.Time.levels)
Cond.Time.levels<-gsub("N.*_C","C",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_T","T",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_M","M",Cond.Time.levels) 
head(Cond.Time.levels)
Cond.Time.levels<-unique(Cond.Time.levels)
Cond.Time.levels
Cond.Time <- factor(paste(alv_mac_norm$samples$group,
                          alv_mac_norm$samples$Time,
                          sep="_"),
                    levels = c(Cond.Time.levels))
head(Cond.Time)
Cond.Time<-relevel(Cond.Time,ref="CN_6H")

test_matrix_6H <- model.matrix(~Animal + Cond.Time,
                               data = alv_mac_norm$samples)

colnames(test_matrix_6H)
dim(test_matrix_6H)

write.table(x         = test_matrix_6H,
            file      = "alv_mac_test_matrix_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

############################################################################
# Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
#                      (for multi-factor data)                             #
############################################################################

alv_mac_disp_6H <- estimateGLMCommonDisp(y    = alv_mac_norm,
                                         design  = test_matrix_6H,
                                         verbose = TRUE)

alv_mac_disp_6H <- estimateGLMTrendedDisp(y   = alv_mac_disp_6H,
                                          design = test_matrix_6H)

alv_mac_disp_6H <- estimateGLMTagwiseDisp(y   = alv_mac_disp_6H,
                                          design = test_matrix_6H)

names(alv_mac_disp_6H)

# Plot the dispersion
png(filename = "BCV_alv_mac_test_matrix_6H.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotBCV(alv_mac_disp_6H)

dev.off()

# Show the calculated dispersion
alv_mac_disp_6H$common.dispersion

# And show its square root, the coefficient of biological variation
sqrt(alv_mac_disp_6H$common.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information #http://seqanswers.com/forums/showthread.php?t=5591
Tagwisedisp_6H <- cbind(alv_mac_disp_6H$genes, alv_mac_disp_6H$tagwise.dispersion)
dim(Tagwisedisp_6H)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp_6H,
             file = "alv_mac_Tagwise_dispersion_test_matrix_6H.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs # ----
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
alv_mac_fit_6H <- glmFit(y = alv_mac_disp_6H, design = test_matrix_6H)
names(alv_mac_fit_6H)
colnames(alv_mac_fit_6H$design)
dim(alv_mac_fit_6H)
# Test for differential expression between the different time points/treatments
MB_6H_lrt <- glmLRT(alv_mac_fit_6H, coef = "Cond.TimeMB_6H")
DE_MB_6H <- topTags(object        = MB_6H_lrt,
                    n             = "inf",
                    adjust.method = "BH")
#View(DE_MB_6H$table)

FDR_0.05_DE_MB_6H <- subset(DE_MB_6H$table, FDR < 0.05)
#View(FDR_0.05_DE_MB_6H)
write.table(x         = FDR_0.05_DE_MB_6H,
            file      = "FDR_0.05_DE_MB_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


TB_6H_lrt <- glmLRT(alv_mac_fit_6H, coef = "Cond.TimeTB_6H")
DE_TB_6H <- topTags(object        = TB_6H_lrt,
                    n             = "inf",
                    adjust.method = "BH")
#View(DE_TB_6H$table)

FDR_0.05_DE_TB_6H <- subset(DE_TB_6H$table, FDR < 0.05)
#View(FDR_0.05_DE_TB_6H)
write.table(x         = FDR_0.05_DE_TB_6H,
            file      = "FDR_0.05_DE_TB_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_MB_6H  <- subset(DE_MB_6H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_MB_6H)
write.table(x         = FDR_0.05_logFC_DE_MB_6H,
            file      = "FDR_0.05_logFC_DE_MB_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_TB_6H  <- subset(DE_TB_6H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_TB_6H)
write.table(x         = FDR_0.05_logFC_DE_TB_6H,
            file      = "FDR_0.05_logFC_DE_TB_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

colnames(alv_mac_fit_6H$coefficients)
#MB-TB
delta_6H_lrt <- glmLRT(alv_mac_fit_6H, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1))
DE_delta_6H <- topTags(object        = delta_6H_lrt,
                       n             = "inf",
                       adjust.method = "BH")
#View(DE_delta_6H$table)

FDR_0.05_delta_6H <- subset(DE_delta_6H$table, FDR < 0.05)
#View(FDR_0.05_delta_6H)
write.table(x         = FDR_0.05_delta_6H,
            file      = "FDR_0.05_delta_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


FDR_0.05_logFC_delta_6H  <- subset(DE_delta_6H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_delta_6H)
write.table(x         = FDR_0.05_logFC_delta_6H,
            file      = "FDR_0.05_logFC_delta_6H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)
##########
# 24H
#########
Animal <-factor(alv_mac_norm$samples$Animal)
head(Animal)

Condition<-factor(alv_mac_norm$samples$Condition, levels=c("CN","MB","TB"))
head(Condition)
Condition<-relevel(Condition,ref="CN")

Time<-factor(alv_mac_norm$samples$Time)
head(Time)
Time<-relevel(Time,ref="24H")

Cond.Time.levels<-rownames(alv_mac_norm$samples)
head(Cond.Time.levels)

Cond.Time.levels<-c(rownames(alv_mac_norm$samples))
head(Cond.Time.levels)
as.data.frame(Cond.Time.levels)
Cond.Time.levels<-gsub("N.*_C","C",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_T","T",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_M","M",Cond.Time.levels) 
head(Cond.Time.levels)
Cond.Time.levels<-unique(Cond.Time.levels)
Cond.Time.levels
Cond.Time <- factor(paste(alv_mac_norm$samples$group,
                          alv_mac_norm$samples$Time,
                          sep="_"),
                    levels = c(Cond.Time.levels))
head(Cond.Time)
Cond.Time<-relevel(Cond.Time,ref="CN_24H")

test_matrix_24H <- model.matrix(~Animal + Cond.Time,
                               data = alv_mac_norm$samples)

colnames(test_matrix_24H)
dim(test_matrix_24H)

write.table(x         = test_matrix_24H,
            file      = "alv_mac_test_matrix_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

############################################################################
# Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
#                      (for multi-factor data)                             #
############################################################################

alv_mac_disp_24H <- estimateGLMCommonDisp(y    = alv_mac_norm,
                                         design  = test_matrix_24H,
                                         verbose = TRUE)

alv_mac_disp_24H <- estimateGLMTrendedDisp(y   = alv_mac_disp_24H,
                                          design = test_matrix_24H)

alv_mac_disp_24H <- estimateGLMTagwiseDisp(y   = alv_mac_disp_24H,
                                          design = test_matrix_24H)

names(alv_mac_disp_24H)

# Plot the dispersion
png(filename = "BCV_alv_mac_test_matrix_24H.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotBCV(alv_mac_disp_24H)

dev.off()

# Show the calculated dispersion
alv_mac_disp_24H$common.dispersion

# And show its square root, the coefficient of biological variation
sqrt(alv_mac_disp_24H$common.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information #http://seqanswers.com/forums/showthread.php?t=5591
Tagwisedisp_24H <- cbind(alv_mac_disp_24H$genes, alv_mac_disp_24H$tagwise.dispersion)
dim(Tagwisedisp_24H)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp_24H,
             file = "alv_mac_Tagwise_dispersion_test_matrix_24H.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs # ----
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
alv_mac_fit_24H <- glmFit(y = alv_mac_disp_24H, design = test_matrix_24H)
names(alv_mac_fit_24H)
colnames(alv_mac_fit_24H$design)
dim(alv_mac_fit_24H)
# Test for differential expression between the different time points/treatments
MB_24H_lrt <- glmLRT(alv_mac_fit_24H, coef = "Cond.TimeMB_24H")
DE_MB_24H <- topTags(object        = MB_24H_lrt,
                    n             = "inf",
                    adjust.method = "BH")
#View(DE_MB_24H$table)

FDR_0.05_DE_MB_24H <- subset(DE_MB_24H$table, FDR < 0.05)
#View(FDR_0.05_DE_MB_24H)
write.table(x         = FDR_0.05_DE_MB_24H,
            file      = "FDR_0.05_DE_MB_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_MB_24H  <- subset(DE_MB_24H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_MB_24H)
write.table(x         = FDR_0.05_logFC_DE_MB_24H,
            file      = "FDR_0.05_logFC_DE_MB_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


TB_24H_lrt <- glmLRT(alv_mac_fit_24H, coef = "Cond.TimeTB_24H")
DE_TB_24H <- topTags(object        = TB_24H_lrt,
                    n             = "inf",
                    adjust.method = "BH")
#View(DE_TB_24H$table)

FDR_0.05_DE_TB_24H <- subset(DE_TB_24H$table, FDR < 0.05)
#View(FDR_0.05_DE_TB_24H)
write.table(x         = FDR_0.05_DE_TB_24H,
            file      = "FDR_0.05_DE_TB_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_TB_24H  <- subset(DE_TB_24H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_TB_24H)
write.table(x         = FDR_0.05_logFC_DE_TB_24H,
            file      = "FDR_0.05_logFC_DE_TB_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

colnames(alv_mac_fit_24H$coefficients)
#MB-TB
delta_24H_lrt <- glmLRT(alv_mac_fit_24H, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0))
DE_delta_24H <- topTags(object        = delta_24H_lrt,
                       n             = "inf",
                       adjust.method = "BH")
#View(DE_delta_24H$table)

FDR_0.05_delta_24H <- subset(DE_delta_24H$table, FDR < 0.05)
#View(FDR_0.05_delta_24H)
write.table(x         = FDR_0.05_delta_24H,
            file      = "FDR_0.05_delta_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_delta_24H  <- subset(DE_delta_24H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_delta_24H)
write.table(x         = FDR_0.05_logFC_delta_24H,
            file      = "FDR_0.05_logFC_delta_24H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

##########
# 48H
#########
Animal <-factor(alv_mac_norm$samples$Animal)
head(Animal)

Condition<-factor(alv_mac_norm$samples$Condition, levels=c("CN","MB","TB"))
head(Condition)
Condition<-relevel(Condition,ref="CN")

Time<-factor(alv_mac_norm$samples$Time)
head(Time)
Time<-relevel(Time,ref="48H")

Cond.Time.levels<-rownames(alv_mac_norm$samples)
head(Cond.Time.levels)

Cond.Time.levels<-c(rownames(alv_mac_norm$samples))
head(Cond.Time.levels)
as.data.frame(Cond.Time.levels)
Cond.Time.levels<-gsub("N.*_C","C",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_T","T",Cond.Time.levels) 
Cond.Time.levels<-gsub("N.*_M","M",Cond.Time.levels) 
head(Cond.Time.levels)
Cond.Time.levels<-unique(Cond.Time.levels)
Cond.Time.levels
Cond.Time <- factor(paste(alv_mac_norm$samples$group,
                          alv_mac_norm$samples$Time,
                          sep="_"),
                    levels = c(Cond.Time.levels))
head(Cond.Time)
Cond.Time<-relevel(Cond.Time,ref="CN_48H")

test_matrix_48H <- model.matrix(~Animal + Cond.Time,
                                data = alv_mac_norm$samples)

colnames(test_matrix_48H)
dim(test_matrix_48H)

write.table(x         = test_matrix_48H,
            file      = "alv_mac_test_matrix_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

############################################################################
# Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
#                      (for multi-factor data)                             #
############################################################################

alv_mac_disp_48H <- estimateGLMCommonDisp(y    = alv_mac_norm,
                                          design  = test_matrix_48H,
                                          verbose = TRUE)

alv_mac_disp_48H <- estimateGLMTrendedDisp(y   = alv_mac_disp_48H,
                                           design = test_matrix_48H)

alv_mac_disp_48H <- estimateGLMTagwiseDisp(y   = alv_mac_disp_48H,
                                           design = test_matrix_48H)

names(alv_mac_disp_48H)

# Plot the dispersion
png(filename = "BCV_alv_mac_test_matrix_48H.png",
    width    = 1366,
    height   = 768,
    units    = "px")

plotBCV(alv_mac_disp_48H)

dev.off()

# Show the calculated dispersion
alv_mac_disp_48H$common.dispersion

# And show its square root, the coefficient of biological variation
sqrt(alv_mac_disp_48H$common.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information #http://seqanswers.com/forums/showthread.php?t=5591
Tagwisedisp_48H <- cbind(alv_mac_disp_48H$genes, alv_mac_disp_48H$tagwise.dispersion)
dim(Tagwisedisp_48H)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp_48H,
             file = "alv_mac_Tagwise_dispersion_test_matrix_48H.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs # ----
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
alv_mac_fit_48H <- glmFit(y = alv_mac_disp_48H, design = test_matrix_48H)
names(alv_mac_fit_48H)
colnames(alv_mac_fit_48H$design)
dim(alv_mac_fit_48H)
# Test for differential expression between the different time points/treatments
MB_48H_lrt <- glmLRT(alv_mac_fit_48H, coef = "Cond.TimeMB_48H")
DE_MB_48H <- topTags(object        = MB_48H_lrt,
                     n             = "inf",
                     adjust.method = "BH")
#View(DE_MB_48H$table)

FDR_0.05_DE_MB_48H <- subset(DE_MB_48H$table, FDR < 0.05)
#View(FDR_0.05_DE_MB_48H)
write.table(x         = FDR_0.05_DE_MB_48H,
            file      = "FDR_0.05_DE_MB_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_MB_48H  <- subset(DE_MB_48H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_MB_48H)
write.table(x         = FDR_0.05_logFC_DE_MB_48H,
            file      = "FDR_0.05_logFC_DE_MB_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


TB_48H_lrt <- glmLRT(alv_mac_fit_48H, coef = "Cond.TimeTB_48H")
DE_TB_48H <- topTags(object        = TB_48H_lrt,
                     n             = "inf",
                     adjust.method = "BH")
#View(DE_TB_48H$table)

FDR_0.05_DE_TB_48H <- subset(DE_TB_48H$table, FDR < 0.05)
#View(FDR_0.05_DE_TB_48H)
write.table(x         = FDR_0.05_DE_TB_48H,
            file      = "FDR_0.05_DE_TB_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_logFC_DE_TB_48H  <- subset(DE_TB_48H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_DE_TB_48H)
write.table(x         = FDR_0.05_logFC_DE_TB_48H,
            file      = "FDR_0.05_logFC_DE_TB_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

colnames(alv_mac_fit_48H$coefficients)
#MB-TB
delta_48H_lrt <- glmLRT(alv_mac_fit_48H, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0))
DE_delta_48H <- topTags(object        = delta_48H_lrt,
                        n             = "inf",
                        adjust.method = "BH")
#View(DE_delta_48H$table)

FDR_0.05_delta_48H <- subset(DE_delta_48H$table, FDR < 0.05)
#View(FDR_0.05_delta_48H)
write.table(x         = FDR_0.05_delta_48H,
            file      = "FDR_0.05_delta_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)



FDR_0.05_logFC_delta_48H  <- subset(DE_delta_48H$table, abs(logFC) > 1 & FDR < 0.05)
#View(FDR_0.05_logFC_delta_48H)
write.table(x         = FDR_0.05_logFC_delta_48H,
            file      = "FDR_0.05_logFC_delta_48H.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)
##########################
# whole DE table results #
##########################
Full_DE_alv_mac <- merge(x  = DE_MB_2H$table,
                         y  = DE_TB_2H$table
                         [,(ncol(DE_TB_2H$table) - 4) :
                             ncol(DE_TB_2H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1] # Correct row names before
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]          # merging the other data frames.

colnames(Full_DE_alv_mac) <- gsub(pattern     = "MB",
                                  replacement = "_MB",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
colnames(Full_DE_alv_mac) <- gsub(pattern     = ".y$",
                                  replacement = "_TB_2hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
head(Full_DE_alv_mac)

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_delta_2H$table
                         [, (ncol(DE_delta_2H$table) - 4) :
                             ncol(DE_delta_2H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_MB_6H$table
                         [, (ncol(DE_MB_6H$table) - 4) :
                             ncol(DE_MB_6H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

colnames(Full_DE_alv_mac) <- gsub(pattern     = ".x$",
                                  replacement = "_delta_2hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
colnames(Full_DE_alv_mac) <- gsub(pattern     = ".y$",
                                  replacement = "_MB_6hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
head(Full_DE_alv_mac)



Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_TB_6H$table
                         [, (ncol(DE_TB_6H$table) - 4) :
                             ncol(DE_TB_6H$table)],
                         by = "row.names")
colnames(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_delta_6H$table
                         [, (ncol(DE_delta_6H$table) - 4) :
                             ncol(DE_delta_6H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

colnames(Full_DE_alv_mac) <- gsub(pattern     = ".x$",
                                  replacement = "_TB_6hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
colnames(Full_DE_alv_mac) <- gsub(pattern     = ".y$",
                                  replacement = "_delta_6hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
head(Full_DE_alv_mac)



Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_MB_24H$table
                         [, (ncol(DE_MB_24H$table) - 4) :
                             ncol(DE_MB_24H$table)],
                         by = "row.names")
colnames(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_TB_24H$table
                         [, (ncol(DE_TB_24H$table) - 4) :
                             ncol(DE_TB_24H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

colnames(Full_DE_alv_mac) <- gsub(pattern     = ".x$",
                                  replacement = "_MB_24hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
colnames(Full_DE_alv_mac) <- gsub(pattern     = ".y$",
                                  replacement = "_TB_24hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
head(Full_DE_alv_mac)


Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_delta_24H$table
                         [, (ncol(DE_delta_24H$table) - 4) :
                             ncol(DE_delta_24H$table)],
                         by = "row.names")
colnames(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_MB_48H$table
                         [, (ncol(DE_MB_48H$table) - 4) :
                             ncol(DE_MB_48H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

colnames(Full_DE_alv_mac) <- gsub(pattern     = ".x$",
                                  replacement = "_delta_24hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
colnames(Full_DE_alv_mac) <- gsub(pattern     = ".y$",
                                  replacement = "_MB_48hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
head(Full_DE_alv_mac)

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_TB_48H$table
                         [, (ncol(DE_TB_48H$table) - 4) :
                             ncol(DE_TB_48H$table)],
                         by = "row.names")
colnames(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

Full_DE_alv_mac <- merge(x  = Full_DE_alv_mac,
                         y  = DE_delta_48H$table
                         [, (ncol(DE_delta_48H$table) - 4) :
                             ncol(DE_delta_48H$table)],
                         by = "row.names")
head(Full_DE_alv_mac)

rownames(Full_DE_alv_mac) <- Full_DE_alv_mac[, 1]
Full_DE_alv_mac <- Full_DE_alv_mac[, -1]

colnames(Full_DE_alv_mac) <- gsub(pattern     = ".x$",
                                  replacement = "_TB_48hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
colnames(Full_DE_alv_mac) <- gsub(pattern     = ".y$",
                                  replacement = "_delta_48hr",
                                  x           = colnames(Full_DE_alv_mac),
                                  perl        = TRUE)
head(Full_DE_alv_mac)


write.table(x         = Full_DE_alv_mac,
            file      = "Full_DE_alv_mac.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

#######################
#Write single results #
#######################
write.table(x         = DE_MB_2H,
            file      = "bovis_2hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_TB_2H,
            file      = "TB_2hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_delta_2H,
            file      = "delta_2hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_MB_6H,
            file      = "bovis_6hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_TB_6H,
            file      = "TB_6hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_delta_6H,
            file      = "delta_6hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_MB_24H,
            file      = "bovis_24hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_TB_24H,
            file      = "TB_24hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_delta_24H,
            file      = "delta_24hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_MB_48H,
            file      = "bovis_48hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_TB_48H,
            file      = "TB_48hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(x         = DE_delta_48H,
            file      = "delta_48hr.txt",
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)


########################################
# Venn diagram of DE genes: FDR < 0.05 # ----
########################################

# Turn gene IDs of significant DE genes (FDR < 0.05)
# per time point into vectors
DE_TB_2H.vector <- c(rownames(FDR_0.05_DE_TB_2H))
DE_MB_2H.vector <- c(rownames(FDR_0.05_DE_MB_2H))
delta_2H.vector <- c(rownames(FDR_0.05_delta_2H))

DE_TB_6H.vector <- c(rownames(FDR_0.05_DE_TB_6H))
DE_MB_6H.vector <- c(rownames(FDR_0.05_DE_MB_6H))
delta_6H.vector <- c(rownames(FDR_0.05_delta_6H))

DE_TB_24H.vector <- c(rownames(FDR_0.05_DE_TB_24H))
DE_MB_24H.vector <- c(rownames(FDR_0.05_DE_MB_24H))
delta_24H.vector <- c(rownames(FDR_0.05_delta_24H))

DE_TB_48H.vector <- c(rownames(FDR_0.05_DE_TB_48H))
DE_MB_48H.vector <- c(rownames(FDR_0.05_DE_MB_48H))
delta_48H.vector <- c(rownames(FDR_0.05_delta_48H))

DE_TB_2H.vector <- c(rownames(FDR_0.05_logFC_DE_TB_2H))
DE_MB_2H.vector <- c(rownames(FDR_0.05_logFC_DE_MB_2H))
delta_2H.vector <- c(rownames(FDR_0.05_logFC_delta_2H))

DE_TB_6H.vector <- c(rownames(FDR_0.05_logFC_DE_TB_6H))
DE_MB_6H.vector <- c(rownames(FDR_0.05_logFC_DE_MB_6H))
delta_6H.vector <- c(rownames(FDR_0.05_logFC_delta_6H))

DE_TB_24H.vector <- c(rownames(FDR_0.05_logFC_DE_TB_24H))
DE_MB_24H.vector <- c(rownames(FDR_0.05_logFC_DE_MB_24H))
delta_24H.vector <- c(rownames(FDR_0.05_logFC_delta_24H))

DE_TB_48H.vector <- c(rownames(FDR_0.05_logFC_DE_TB_48H))
DE_MB_48H.vector <- c(rownames(FDR_0.05_logFC_DE_MB_48H))
delta_48H.vector <- c(rownames(FDR_0.05_logFC_delta_48H))


# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.plot <- venn.diagram(list(A = DE_TB_2H.vector,
                               B = DE_MB_2H.vector,
                               C = delta_2H.vector),
                          filename        = "Venn_DE_FDR_0-05.png",
                          imagetype       = "png",
                          col             = "transparent",
                          fill            = c("#ffffcc",
                                              "#225ea8",
                                              "#a1dab4"),
                          alpha           = 0.50,
                          label.col       = "#003333",
                          cex             = 10,
                          fontfamily      = "Raavi",
                          category.names  = c("MB",
                                              "TB",
                                              "delta"),
                          cat.col         = "black",
                          cat.cex         = 10,
                          cat.pos         = c(-11, 11, 0),
                          cat.dist        = c(0.21, 0.21, 0.1),
                          #cat.fontfamily  = "Raavi",
                          rotation.degree = 360,
                          margin          = 0,
                          height          = 85,
                          width           = 85,
                          units           = 'cm',
                          compression     = 'lzw',
                          resolution      = 300)


library("gplots")
venn(list(A=DE_TB_2H.vector,B=DE_MB_2H.vector,C=delta_2H.vector))
venn(list(A=DE_TB_6H.vector,B=DE_MB_6H.vector,C=delta_6H.vector))
venn(list(A=DE_TB_24H.vector,B=DE_MB_24H.vector,C=delta_24H.vector))
venn(list(A=DE_TB_48H.vector,B=DE_MB_48H.vector,C=delta_48H.vector))


#install.packages("sigora")
library("sigora")
# Set working directory and load any previously saved data
setwd("/Users/Kerri/Google Drive/Postdoc UCD /Alv mac work/EdgeR/")
sig.data<-read.csv("FDR_0.05_logFC_delta_48H.txt",sep="\t",header=TRUE)
sig.data.MB<-subset(sig.data, logFC >1)
names(sig.data.MB)
sig.data.MB.genes<-(sig.data.MB[,3])
sig.data.MB.genes<-na.omit(sig.data.MB.genes)
data(kegH)
data(nciTable)
head(nciTable)
sigRes<-sigora(kegH,level=4, queryList =sig.data.MB.genes, saveFile="myResultsMBtest.csv" )
sig.data.TB<-subset(sig.data, logFC < 1)
names(sig.data.TB)
sig.data.TB.genes<-(sig.data.TB[,3])
sig.data.TB.genes<-na.omit(sig.data.TB.genes)
data(kegH)
data(nciTable)
head(nciTable)
sigRes<-sigora(kegH,level=4, queryList =sig.data.TB.genes, saveFile="myResultsTBtest2.csv" )
########
#SIGORA
#######
install.packages("sigora")
library(sigora)
data(kegH)
a1<-genesFromRandomPathways(seed=12345,kegH,3,50)
## originally selected pathways:
a1[["selectedPathways"]]
a1[["genes"]]
## Traditional ora identifies dozens of statistically significant pathways!
ora(a1[["genes"]],kegH)
## Now let us try sigora with the same input:
sigoraRes <- sigora(GPSrepo =kegH, queryList = a1[["genes"]],level = 4)
## Again, the three originally selected pathways were:
a1[["selectedPathways"]]

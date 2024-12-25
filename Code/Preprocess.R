#!/usr/bin/env Rscript

library(ChAMP)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(limma)
library(doParallel)

#Set/create directory
wd <- "/home/ycwang27910" #My working directory
idatdir <- "/home/ycwang27910/IDAT_files_133" #Directory where IDAT files are in
resdir <- "/home/ycwang27910/FIN" #Directory for output files
dt <- format(Sys.time(), "%Y%m%d") #System date
createlog <- TRUE #Generate log files or not?

setwd(wd)
if(!dir.exists(resdir)){
      dir.create(resdir, showWarnings = FALSE, recursive = TRUE)  
}

if(createlog){
        logfile <- paste0(resdir, sprintf("/preprocess_log_%s.txt", format(Sys.time(), "%Y%m%d-%H%M%S")))
        con <- file(logfile, open = "wt")
        sink(con, type="output")
        sink(con, type="message")        
}

#Age related probes from Horvath (Genome Biol, 2013) and Hannum et al. (Mol Cell, 2013) 
rmprobe <- as.vector(unlist(read.delim(file.path(wd, "ProbesToRemove_pub.txt"), header = F)))

##################################
#Start data preprocessing on NCHC#
##################################                               

#Import idat files
myImport <- champ.import(directory = idatdir, arraytype = "EPIC")
#Load idat files with minfi (For Noob normalizaton)
myLoad <- champ.load(directory = idatdir, arraytype = "EPIC", method = "minfi", force = T)
#Normalized with Noob method and get beta value
noob <- preprocessNoob(myLoad$rgSet)
Mset_raw.betavalue <- getBeta(noob)

#Intersect CpGs with imported beta
rowname_overlap <- intersect(rownames(Mset_raw.betavalue),
                             rownames(myImport$beta)) #To ensure the probes imported by ChAMP method are identical to those imported by minfi method  
#Subset the overlapped rows from imported data
myImport$beta_subset <- subset(myImport$beta,
                               rownames(myImport$beta) %in% rowname_overlap)
myImport$detP_subset <- subset(myImport$detP,
                               rownames(myImport$beta) %in% rowname_overlap)
myImport$beadcount_subset <- subset(myImport$beadcount,
                                    rownames(myImport$beta) %in% rowname_overlap)
myImport$M_subset <- subset(myImport$M,
                            rownames(myImport$beta) %in% rowname_overlap)
norm_subset <- subset(Mset_raw.betavalue,
                      rownames(Mset_raw.betavalue) %in% rowname_overlap)
#Sort to match beta table and other files
norm_subset_sort <- norm_subset[order(row.names(norm_subset)), ]
myImport$beta_subset_sort <- myImport$beta_subset[order(row.names(myImport$beta_subset)), ]
myImport$detP_subset_sort <- myImport$detP_subset[order(row.names(myImport$detP_subset)), ]
myImport$beadcount_subset_sort <- myImport$beadcount_subset[order(row.names(myImport$beadcount_subset)), ]
myImport$M_subset_sort <- myImport$M_subset[order(row.names(myImport$M_subset)), ]
colnames(norm_subset_sort) <- colnames(myImport$beta_subset_sort)


#Filter the intersected CpGs and manually remove age-related probes
myFilter <- champ.filter(beta=norm_subset_sort,M=myImport$M_subset_sort,pd=myImport$pd,detP=myImport$detP_subset_sort,beadcount=myImport$beadcount_subset_sort,autoimpute=TRUE, filterDetP=TRUE, ProbeCutoff=0, SampleCutoff=0.1, detPcut=0.01, filterBeads=TRUE, beadCutoff=0.05, filterNoCG = TRUE,filterSNPs = TRUE, population = "EAS", filterMultiHit = TRUE, filterXY = TRUE, fixOutlier = TRUE, arraytype = "EPIC")
myFilter <- myFilter$beta[!rownames(myFilter$beta) %in% rmprobe,]
#Normalized with BMIQ method 
NBBM <- champ.norm(beta = myFilter, arraytype = "EPIC")
#Correct batch effects from slides and arrays 
myCombat <- champ.runCombat(beta = NBBM, pd = myImport$pd,
                            batchname=c("Slide","Array"))
#Save preprocessed beta values
saveRDS(myCombat, file = file.path(resdir, sprintf("NBBM_CombatSA_rmAge_%s.rds", dt)))

#Quality Control plots
champ.QC(beta = myCombat, resultsDir = file.path(resdir, sprintf("CHAMP_QCimages_AfterCombatSA_%s", dt)))
champ.SVD(beta = as.data.frame(myCombat), pd = myLoad$pd,
          resultsDir = file.path(resdir, sprintf("CHAMP_QCimages_AfterCombatSA_%s", dt)))

######################
#Exploratory analysis#
######################

#get PCA1 and PCA2 of preprocessed beta matrix
b <- as.matrix(myCombat)
pca <- prcomp(t(b), scale = TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)[1:10]
pca_df <- pca$x[,1:10]
mypca <- list(var.explained = pca.var.per, pca = pca_df)

saveRDS(mypca, file = file.path(resdir, sprintf("PCA_NBBM_%s.rds", dt)))

#Find DMPs corrected by age  

Compare <- c("HV","LC")
#FIRST THING AFTER ~ IS PHENOTYPE TO BE COMPARED AND EVERYTHING AFTER ARE ADJUSTED COVARIATES
pheno <- myLoad$pd[,c("Sample_Name","Age","Sample_Group")]
pheno$Age <- as.numeric(pheno$Age)
design <- model.matrix(~ Sample_Group + Age, data=pheno) ##Design matrix for DMP AND DMR
fit1 <- lmFit(myCombat,design)
fit2 <- eBayes(fit1)
IV=colnames(fit2$coefficients)[2]
#Differentially methylated probes
DMP <- topTable(fit2, coef = IV, adjust.method = "BH",sort.by = "P", num = Inf)
message("  You have found ",sum(DMP$adj.P.Val < 0.05), 
        " significant DMPs with a BH adjusted P-value below ", 
        0.05,".")
#Create beta table with other genomic information
data("probe.features.epic")
com.idx <- intersect(rownames(DMP),rownames(probe.features))
avg <-  cbind(rowMeans(myCombat[com.idx,which(pheno$Sample_Group==Compare[1])]),
              rowMeans(myCombat[com.idx,which(pheno$Sample_Group==Compare[2])]))
avg <- cbind(avg,avg[,2]-avg[,1]) # LC - HV
colnames(avg) <- c(paste("HV","AVG",sep="_"), paste("LC","AVG",sep="_"), "deltaBeta")
myDMP <- data.frame(DMP[com.idx,],avg,probe.features[com.idx,])

saveRDS(myDMP, file = file.path(resdir, sprintf("DMP_All_CombatSA_rmAge_%s.rds", dt)))

#Find DMRs corrected by age

RSobject <- RatioSet(myCombat, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b4.hg19"))
probe.features <- getAnnotation(RSobject)
message("<< Find DMR with Bumphunter Method >>")
registerDoParallel(cores = 3)
cpg.idx <- intersect(rownames(myCombat),rownames(probe.features))
Anno <- probe.features[cpg.idx,]
Anno <- Anno[order(Anno$chr,Anno$pos),]
cpg.idx <- rownames(Anno)

cl <- clusterMaker(Anno$chr,Anno$pos,maxGap=300)
names(cl) <- cpg.idx
bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl)>=2)))]

message("According to your data set, ",
        sum(table(cl)>=2)," clusters contains MORE THAN ",
        2," probes within ",300,
        " maxGap were detected. These clusters will be used to find DMR.\n")

Beta <- myCombat[bumphunter.idx,]
Beta <- replace(Beta,which(Beta <= 0.001),0.001)
Beta <- replace(Beta,which(Beta >= 0.999),0.999)
Y <- log((Beta/(1-Beta)),2)

set.seed(1030)
Bumps <- bumphunter(as.matrix(Y),
                    design=design,
                    chr=Anno[bumphunter.idx,]$chr,
                    pos=Anno[bumphunter.idx,]$pos,
                    cluster=cl[bumphunter.idx],
                    cutoff=NULL,
                    pickCutoff=TRUE,
                    smooth=TRUE,
                    smoothFunction=loessByCluster,
                    useWeights=FALSE,
                    permutations=NULL,
                    verbose=TRUE,
                    B=250,
                    nullMethod="bootstrap")
message("<< Calculate DMR success. >>")
DMR <- Bumps$table[which(Bumps$table$p.valueArea <= 0.05),]
message("Bumphunter detected ",nrow(DMR)," DMRs with P value <= ",0.05,".")

rownames(DMR) <- paste("DMR",1:nrow(DMR),sep="_")
DMR <- data.frame(DMR[,1:3],width=DMR[,3]-DMR[,2],strand="*",DMR[,4:14])
colnames(DMR)[1:3] <- c("seqnames","start","end") 

saveRDS(DMR, file = file.path(resdir, sprintf("DMR_CombatSA_rmAge_%s.rds", dt)))

#Get beta tables for visualization on my computer
#Count probes in each methylation levels
Chart_globalmeth <- data.frame(
        CpG_Methylation_Level=c(">70%", "30-70%", "<30%")
)

for(i in 1:ncol(myCombat)){
        
        x <- myCombat[,i]
        n <- colnames(myCombat)[i]
        H <- sum(x >= 0.7)
        M <- sum(0.3 <= x & x <= 0.7)
        L <- sum(x <= 0.3)
        c <- data.frame(c(H, M, L))
        names(c) <- n
        Chart_globalmeth <- cbind(Chart_globalmeth, c)
        
}

saveRDS(Chart_globalmeth, file.path(resdir, sprintf("GlobalMeth_All_%s.rds", dt)))

#Get beta table of Significant probes
GOMP <- myDMP %>% filter(LC_AVG > HV_AVG) %>% 
        filter(abs(deltaBeta) > 0.05, adj.P.Val < 0.05) %>% rownames() #Gain of methylation probes
LOMP <- myDMP %>% filter(LC_AVG < HV_AVG) %>%
        filter(abs(deltaBeta) > 0.05, adj.P.Val < 0.05) %>% rownames() #Loss of methylation probes
sigprobes <- c(GOMP, LOMP)

mybeta_sig <- myCombat[rownames(myCombat) %in% sigprobes,]

saveRDS(mybeta_sig, file.path(resdir, sprintf("SigBetaTable_%s.rds", dt)))





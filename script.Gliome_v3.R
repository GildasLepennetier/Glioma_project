
rm(list=ls())
source("https://raw.githubusercontent.com/GildasLepennetier/R-stuff/master/not_in.R")

setwd("~/LRZ Sync+Share/PROJECT-GILDAS/SCRIPTS/2021.06.25.RNAseq_Gliome_samples")

#### LOAD PACKAGES ####
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(BiocParallel)
library(matrixStats) # matrixStats::rowVars  << BiocManager::install("matrixStats")
options(max.print=150)
register(MulticoreParam(1))



#### ------- annot -------- ####
ANNOT<-read.xlsx("333TK_Human02_SampleAnnotation.xlsx")
rownames(ANNOT)<-ANNOT$Barcode
#select only barcode of interest
ANNOT<-subset(ANNOT, Loaded == 1)
ANNOT$intgroup = paste(ANNOT$Descr_Organ, ANNOT$Descr_Cell.Type)
ANNOT$Descr_SAP.ID <- factor(ANNOT$Descr_SAP.ID)
ANNOT$TEST=factor(paste(ANNOT$Descr_Organ,ANNOT$Descr_Cell.Type,sep="_"))

#### ------- Filter: remove Patient GMB 003, first line with weird cell count -------- ####
#ANNOT[ , c("Descr_GBM.ID","Batch","intgroup") ]
ANNOT<-ANNOT[ ! (ANNOT$TechnicalRep == "TR1" & ANNOT$Descr_GBM.ID == "GBM03") , ]


#### ------- gene matrix 1 -------- ####
COUNTS1 <- read.table("333TK_Human01_DGE_Matrix.txt",header=T,sep="\t");COUNT1_NAME="333TK_Human01"
rownames(COUNTS1)<-COUNTS1$GENE ; COUNTS1$GENE<-NULL
#filter from the annot
COUNTS1<-COUNTS1[,ANNOT$Barcode]
colnames(COUNTS1) <- paste0(colnames(COUNTS1),"-1")
#### ------- gene matrix 2 -------- ####
COUNTS2 <- read.table("333TK_Human02_DGE_Matrix.txt",header=T,sep="\t");COUNT2_NAME="333TK_Human02"
rownames(COUNTS2)<-COUNTS2$GENE ; COUNTS2$GENE<-NULL
COUNTS2<-COUNTS2[,ANNOT$Barcode]
colnames(COUNTS2) <- paste0(colnames(COUNTS2),"-2")







#### ------- Klemm 2020 -------- ####
if(T){
	COUNTS3 <- read.table("BrainTIME_RawCounts.csv",header=T,sep=",",row.names=1);COUNT3_NAME="BrainTIME"
	# ANNOT3<-read.csv("BrainTIME_ClinicalAnnotation.csv")
	# table(ANNOT3$Tissue.Type)
	# ANNOT3$Tissue.Type<-factor(ANNOT3$Tissue.Type)
	# levels(ANNOT3$Tissue.Type) <- c("BrM","Glioma","HealthyBlood","NonTumor") #re-write to avoit spaces and stuff not nice in formula
	
	#The annotation actually do not fit the counts header:
	ANNOT3 <- data.frame(Patient=colnames(COUNTS3))
	ANNOT3 <-  tidyr::separate(ANNOT3,col="Patient",into=c("Location","TissueType","Number","CellType"),sep="_",remove=F)
	
	ANNOT3$TEST=factor(paste(ANNOT3$TissueType, ANNOT3$CellType,sep="_"))
	#write.xlsx(ANNOT3,file="BrainTIME_FromCounts_Gildas.xlsx")
	
	table(ANNOT3$Location, ANNOT3$CellType)
	
}
# mdm:	monocyte-derived macrophages
# m:	microglia
# n:	non-tumor (brain) tissue
# brm:	brain metastase
# glioma:	brain tumor



#### ------- save each file in excel ------- ####
#write.xlsx(COUNTS1,"TUM_Human01_ExprMatrix.xlsx")
#write.xlsx(COUNTS2,"TUM_Human02_ExprMatrix.xlsx")
#write.xlsx(ANNOT,"TUM_Annotations.xlsx")



#### ------- barplot, count of genes with  >10 hits -------####
if(F){
	# count the number of genes (rows) with more than 10 as count
	# carefull: do not take the GENE column
	
	for( index in 3) {
		if(index==1) Countmat<-COUNTS1 ; FILENAME="GeneCount.Human_01.gt10Hits.svg" ; NAMES<-ANNOT$UniqueSampleID
		if(index==2) Countmat<-COUNTS2 ; FILENAME="GeneCount.Human_02.gt10Hits.svg" ; NAMES<-ANNOT$UniqueSampleID
		if(index==3) Countmat<-COUNTS3 ; FILENAME="GeneCount.BrainTIME.gt10Hits.svg" ; NAMES<-colnames(COUNTS3)
	
		NB_GENES=nrow(Countmat)
		LIMIT_ACCEPTABLE=NB_GENES/10
		X <- apply(X=Countmat,MARGIN = 2, FUN = function(x){sum( x>10 )}  )
		names(X) <- NAMES
		svg(filename = FILENAME,width = 18, height = 8)
		par(mar=c(18,4,4,4)+.1)
		COL=rep("gray",length(X))
		COL[X<LIMIT_ACCEPTABLE] <-"red"
		barplot(X,col=COL,las=2,main= paste0("Gene with >10 hits\n","nb genes total: ",NB_GENES," | arbitrary threashold at 10%"))
		abline(h=LIMIT_ACCEPTABLE, col="red")
		dev.off()
		
		rm(Countmat,LIMIT_ACCEPTABLE,COL,X,NB_GENES,FILENAME)
	}
}
#### ------- >>> [ACTUALLY DO NOT IGNORE THEM -> do not run this block ] Ignore Mito and Ribo genes (alternative pipeline) -------####
if(F){
	TAG_SELECT_GENES="noRiboMito"
	
	GENES <- rownames(COUNTS1) #48520
	GENES <- grep(pattern = "^RP[SL]",GENES,value=T,invert=T)
	GENES <- grep(pattern = "^MT-",GENES,value=T,invert=T) #47217
	COUNTS1 <- COUNTS1[ GENES,  ] ; rm(GENES)
	
	GENES <- rownames(COUNTS2) #45579
	GENES <- grep(pattern = "^RP[SL]",GENES,value=T,invert=T)
	GENES <- grep(pattern = "^MT-",GENES,value=T,invert=T) #44246
	COUNTS2 <- COUNTS2[ GENES,  ] ; rm(GENES)
}


#### ------- DESeq2 -------- ####
if(F){
	# Pick the reference in order to use apeglm with the coef option
	register(MulticoreParam(14))
	
	#ANNOT$TEST <- relevel(ANNOT$TEST,ref="NormalBrain_CD4")
	#ANNOT$TEST <- relevel(ANNOT$TEST,ref="NormalBrain_CD8")
	#ANNOT$TEST <- relevel(ANNOT$TEST,ref="PeripheralBlood_CD4")
	#ANNOT$TEST <- relevel(ANNOT$TEST,ref="PeripheralBlood_CD8")
	ANNOT$TEST <- relevel(ANNOT$TEST,ref="Tumor_CD8")
	
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = COUNTS1, colData = ANNOT,design= ~ TEST)
	dds <- DESeq2::DESeq(dds,parallel=T)
	dds1=dds ; rm(dds)
	
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = COUNTS2, colData = ANNOT,design= ~ TEST)
	dds <- DESeq2::DESeq(dds,parallel=T)
	dds2=dds ; rm(dds)
	
	
	ANNOT3$TEST <- relevel(ANNOT3$TEST,ref="glioma_cd8")
	#ANNOT3$TEST <- relevel(ANNOT3$TEST,ref="healthyBlood_cd4")
	#ANNOT3$TEST <- relevel(ANNOT3$TEST,ref="healthyBlood_cd8")
	
	#"glioma_cd8","brm_cd4","brm_cd45n","brm_cd8","brm_mdm","brm_mg","brm_neutrophils","glioma_cd4","glioma_cd45n","glioma_mdm",
	#"glioma_mg","glioma_neutrophils","healthyBlood_cd4","healthyBlood_cd8","healthyBlood_mdm","healthyBlood_neutrophils","nonTumor_cd45n","nonTumor_mg"
	
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = COUNTS3, colData = ANNOT3,design= ~ TEST)
	dds <- DESeq2::DESeq(dds,parallel=T)
	dds3=dds ; rm(dds)
	register(MulticoreParam(1))
}


# to do: TUM foldchange versus Braintime foldchange - OVERLAP TOOL: https://www.molbiotools.com/listcompare.php
if(F){
	# possible comparison:
	#TUM_vs_BrainTIME_CD4_T_vs_B_ovelap_993genes.txt
	# OLD GENES_select=read.table("03.DEG.WithRiboMito/CD4/TUM_vs_BrainTIME_CD4_T_vs_B_overlap_993genes.txt",header=F)
	GENES_select=read.table("03.DEG.WithRiboMito/CD4/TUM_vs_BrainTIME_CD4_T_vs_B_overlap_109genes.txt",header=F)
	GENES_select=GENES_select$V1
	DATA1=read.xlsx("03.DEG.WithRiboMito/CD4/DEG.333TK_Human01.TEST_Tumor_CD4_vs_PeripheralBlood_CD4.xlsx")
	DATA3=read.xlsx("03.DEG.WithRiboMito/CD4/DEG.BrainTIME.TEST_glioma_cd4_vs_healthyBlood_cd4.xlsx")
	
	#TUM_vs_BrainTIME_CD8_T_vs_B_ovelap_375genes.txt
	# OLD GENES_select=read.table("03.DEG.WithRiboMito/CD8/TUM_vs_BrainTIME_CD8_T_vs_B_overlap_375genes.txt",header=F)
	GENES_select=read.table("03.DEG.WithRiboMito/CD8/TUM_vs_BrainTIME_CD8_T_vs_B_overlap_328genes.txt",header=F)
	GENES_select=GENES_select$V1
	DATA1=read.xlsx("03.DEG.WithRiboMito/CD8/DEG.333TK_Human01.TEST_Tumor_CD8_vs_PeripheralBlood_CD8.xlsx")
	DATA3=read.xlsx("03.DEG.WithRiboMito/CD8/DEG.BrainTIME.TEST_glioma_cd8_vs_healthyBlood_cd8.xlsx")
	
	
	
	
	
	GENES=intersect(DATA1$GeneID,DATA3$GeneID)
	print(paste(length(GENES),"genes exists in both data"))
	#17626 genes exists in both data (same for CD4 and CD8)
	
	#reorganise
	DATA1<-subset(DATA1, GeneID %in% GENES)
	DATA3<-subset(DATA3, GeneID %in% GENES)
	rownames(DATA1)<-DATA1$GeneID
	rownames(DATA3)<-DATA3$GeneID
	DATA1<-DATA1[GENES,]
	DATA3<-DATA3[GENES,]
	
	
	
	DATA=data.frame(GeneID=GENES,
					FC_TUM=DATA1$log2FoldChange,
					FC_BrainTIME=DATA3$log2FoldChange)
	rownames(DATA)<-DATA$GeneID

	DATA$direction_TUM=ifelse(DATA$FC_TUM>0,"up","down")
	DATA$direction_BrainTIME=ifelse(DATA$FC_BrainTIME>0,"up","down")
	
	#DATA$overlap_993<-DATA$GeneID %in% GENES_select
	DATA$overlap<-DATA$GeneID %in% GENES_select
	
	GENE_LIST_CD4=c("BATF", "CARD16", "CCL20", "CCR5", "CCR10", "CD58", "CD63", "CD70", "CD74", "CD96", "CEACAM1", "CSF2", "CSF1R", 
					"CTLA4", "FOXP1", "IFNG", "IFNGR1", "IL10", "IL18R1", "IL6ST", "ITGA1", "ITGAE", "TCF7", "TIGIT", "TNFRSF9", 
					"TNFSF9", "TOX", "TREM2") #12 13B, 14, 
	DATA$gene_to_annotate=DATA$GeneID %in% GENE_LIST_CD4 # OLD DATA1$GeneID[1:25]
	
	
	
	GENE_LIST_CD8=c("ADAM19", "BACH2", "CCR1", "CCR5", "CCR7", "CCL4", "CD38", "CD81", "CXCR6", "CTLA4", "FOXP1", "IFNG", "ITGA6", 
					"ITGAE", "LEF1", "TCF7", "TNF", "TFNRSF1B", "TNFSF9", "TOX") #9, 25, 10, 14,
	DATA$gene_to_annotate=DATA$GeneID %in% GENE_LIST_CD8 # OLD DATA1$GeneID[1:25]
	
	#table(DATA$gene_to_annotate)
	#DATA[DATA$gene_to_annotate,]
	
	#DATA$color=ifelse( DATA$overlap_993,"overlap","no overlap")
	DATA$color=ifelse( DATA$overlap,"overlap","no overlap")
	
	DATA$color<-factor(DATA$color,levels=c("overlap","no overlap"))
	#74 genes actually have some NA as fold change
	#Removed 74 rows containing missing values (geom_point).
	
	#plot only when there is an overlap AND there ar ein the list
	DATA$gene_to_annotate <- DATA$gene_to_annotate & DATA$overlap
	DATA$GeneID[DATA$gene_to_annotate & DATA$overlap]
	#for the CD4 #  "CCR5"   "IFNG"   "CD74"   "TNFSF9" "FOXP1"
	#for the CD8 #   "CCR5"   "CTLA4"  "TCF7"   "LEF1"   "IFNG"   "CXCR6"  "TNF"    "ITGA6"  "TNFSF9" "CD38"   "CCR7"   "CCL4"   "ADAM19" "CCR1"   "ITGAE" "TOX"    "FOXP1"  "CD81"
	
	PLOT=ggplot(DATA) + 
		geom_point(aes(x=FC_TUM,y=FC_BrainTIME,color=color),data=subset(DATA,!overlap),size=0.3,alpha=0.3) + 
		geom_point(aes(x=FC_TUM,y=FC_BrainTIME,color=color),data=subset(DATA,overlap),size=1,alpha=1) + 
		scale_fill_manual(aesthetics="color",values=c("no overlap"="grey","overlap"="darkred")) +
		xlab("TUM") + ylab("BrainTIME") + ggtitle("Fold change of gene expression") +
		theme_minimal(26) + theme(legend.title=element_blank(),legend.position="none") #coord_fixed(ratio = 0.25) + 
	pdf("CDx_relation_plot_all.pdf",10,10);print(PLOT);dev.off()
	write.xlsx(DATA,"CDx_relation_plot_all.xlsx",overwrite = T)
	
	PLOT=ggplot(DATA) + 
		geom_point(aes(x=FC_TUM,y=FC_BrainTIME,color=color),data=subset(DATA,!overlap),size=0.3,alpha=0.3) + 
		geom_point(aes(x=FC_TUM,y=FC_BrainTIME,color=color),data=subset(DATA,overlap),size=1,alpha=1) + 
		scale_fill_manual(aesthetics="color",values=c("no overlap"="grey","overlap"="darkred")) +
		xlab("TUM") + ylab("BrainTIME") + ggtitle("Fold change of gene expression") +
		theme_minimal(26) + theme(legend.title=element_blank(),legend.position="none") + #coord_fixed(ratio = 0.25) + 
		ggrepel::geom_label_repel(aes(x=FC_TUM,y=FC_BrainTIME,label=GeneID),max.overlaps=100,
								  data=subset(DATA, gene_to_annotate),nudge_x=0,nudge_y=1)
	pdf("CDx_relation_plot_labels.pdf",10,10);print(PLOT);dev.off()

	
	
	
	
	# SUBSET_over_signif=subset(DATA, overlap)
	# 
	# TEST=cor.test(DATA$FC_TUM, DATA$FC_BrainTIME, method = "kendall")
	# TEST
	# TEST$p.value
	# TABLE=table(TUM=DATA$direction_TUM,BrainTIME=DATA$direction_BrainTIME);TABLE
	# round(TABLE/sum(TABLE)*100,2)
	# 
	# TEST=cor.test(SUBSET_over_signif$FC_TUM, SUBSET_over_signif$FC_BrainTIME, method = "kendall")
	# TEST
	# TEST$p.value
	# TABLE=table(TUM=SUBSET_over_signif$direction_TUM,BrainTIME=SUBSET_over_signif$direction_BrainTIME);TABLE
	# round(TABLE/sum(TABLE)*100,2)
	
}
#### Prepare count matrices for heatmap ####
if(F){
	# copy data COUNT
	C1 <- counts(dds1, normalized=T) #as.matrix(COUNTS1)#C1 <- C1[ base::rowSums(C1) > 10 , ]
	colnames(C1) <- paste0(colnames(C1),"-1")
	
	C2 <- counts(dds2, normalized=T) #as.matrix(COUNTS2)#C2 <- C2[ base::rowSums(C2) > 10 , ]
	colnames(C2) <- paste0(colnames(C2),"-2")
	
	#graphics::boxplot( as.vector(C1), as.vector(C2), outline=F)
	
	# Z-score = x - mean( X ) / sd( X )
	TAG_Norm="_Zscore"
	NC1 <- ( C1 - mean(C1,na.rm=T) ) / sd(C1) # center: x - mean of the column / sd
	NC2 <- ( C2 - mean(C2,na.rm=T) ) / sd(C2) # center: x - mean of the column / sd
	
	COM_GENES = intersect( rownames(C1), rownames(C2) ) ; print(paste(length(COM_GENES), "common genes between runs"))
	#29962 common genes between runs

	C <- cbind(C1[COM_GENES,],C2[COM_GENES,])
	NC <- cbind(NC1[COM_GENES,],NC2[COM_GENES,])
	rm(C1,C2)
}
#### annotation data for heatmaps: since we take both run at the same time, we have to duplicate the data ####
if(F){
	AnnotDf = data.frame( Batch=c( rep("1",nrow(ANNOT)), rep("2",nrow(ANNOT)) ),
						  Barcode = c( paste0( ANNOT$Barcode, "-1") , paste0( ANNOT$Barcode, "-2")),
						  Descr_Organ = c(ANNOT$Descr_Organ, ANNOT$Descr_Organ),
						  Descr_Cell.Type = c(ANNOT$Descr_Cell.Type, ANNOT$Descr_Cell.Type),
						  UniqueSampleID= c(ANNOT$Descr_GBM.ID, ANNOT$Descr_GBM.ID ))
	rownames(AnnotDf)<-AnnotDf$Barcode ; AnnotDf$Barcode=NULL
}
#### ------- PCA plot-------- ####
if(F){
	# Manually select
	for( index in 1:3 ) {
		#INTGROUPS colData(dds)
		if(index==1) {dds=dds1 ; COUNT_NAME=COUNT1_NAME ; INTGROUPS= c("intgroup","Descr_Organ","Descr_Cell.Type","Descr_SAP.ID")}
		if(index==2) {dds=dds2 ; COUNT_NAME=COUNT2_NAME ; INTGROUPS= c("intgroup","Descr_Organ","Descr_Cell.Type","Descr_SAP.ID")}
		if(index==3) {dds=dds3 ; COUNT_NAME=COUNT3_NAME ; INTGROUPS= c("TEST","TissueType","CellType","Location")}
		
		
		vsd <- DESeq2::vst(dds, blind=T) #rld <- rlog(dds, blind=T)
		NTOP=500
		
		for ( INTGROUP in INTGROUPS ) {
			
			# mod 1
			if(F){
				#ntop =  how many of the most variable genes should be used in calculating the PCA 
				DATA<-BiocGenerics::plotPCA(vsd,intgroup=INTGROUP,ntop=NTOP,returnData=T)
				PLOT<-ggplot(DATA) + 
					geom_point(aes(PC1,PC2,col=group),size=8,alpha = 0.7) +
					theme_minimal(16) + 
					theme(legend.title = element_blank()) +
					coord_fixed()
				OUTNAME=paste("PCA",COUNT_NAME,INTGROUP,paste0("Top",NTOP),"pdf",sep=".")
				pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
			}
			
			
			# Extra PCA with PC3
			if(T){
				MATRIX=as.matrix(assay(vsd))
				variance = MatrixGenerics::rowVars(MATRIX) #Calculates the variance for each row (column) of a matrix-like object
				MostVar=order(variance,decreasing = T)[1:NTOP]
				#variance[MostVar]
				#Select the top most variable also
				pca <- stats::prcomp(MATRIX[MostVar,], scale.=TRUE)
				
				summary(pca)
				
				#plot( pca$sdev[1:10], pch=19, cex=2, col="darkblue", xlab="PC", ylab="Standard deviation")
				plot( summary(pca)$importance[2,1:10], pch=19, cex=2, col="darkblue", xlab="PC", ylab="Proportion of Variance")
				
				DATA=as.data.frame(pca$rotation)
				
				if ( any( rownames(DATA$rotation) != colnames( dds ) ) ) { 
					print ( paste("Error in the data: samples not in the same order", INTGROUP, COUNT_NAME))
				}
				else
					{
					DATA$group <- dds[[INTGROUP]]
				
					PLOT<-ggplot(DATA) + 
						geom_point(aes(PC1,PC2,col=group),size=8,alpha = 0.7) +
						theme_minimal(16) + 
						theme(legend.title = element_blank()) 
					OUTNAME=paste("PCA-pc1-2",COUNT_NAME,INTGROUP,paste0("Top",NTOP),"pdf",sep=".")
					pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
				
					PLOT<-ggplot(DATA) + 
						geom_point(aes(PC1,PC3,col=group),size=8,alpha = 0.7) +
						theme_minimal(16) + 
						theme(legend.title = element_blank()) 
					OUTNAME=paste("PCA-pc1-3",COUNT_NAME,INTGROUP,paste0("Top",NTOP),"pdf",sep=".")
					pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
					
					PLOT<-ggplot(DATA) + 
						geom_point(aes(PC2,PC3,col=group),size=8,alpha = 0.7) +
						theme_minimal(16) + 
						theme(legend.title = element_blank()) 
					OUTNAME=paste("PCA-pc2-3",COUNT_NAME,INTGROUP,paste0("Top",NTOP),"pdf",sep=".")
					pdf(OUTNAME,width=10,height=8,useDingbats=F);print(PLOT);dev.off()
					
					#https://tem11010.github.io/Plotting-PCAs/
					
					#library(pca3d)
					#pca3d(pca, group=dds[[INTGROUP]],)
					#OUTNAME=paste("PCA-pc1-2-3",COUNT_NAME,INTGROUP,"pdf",sep=".")
					#snapshotPCA3d(OUTNAME)
					# not working This build of rgl does not include OpenGL functions.
					
					library(plotly)
					plot_ly(x=DATA$PC1, y=DATA$PC2, z=DATA$PC3, type="scatter3d", mode="markers", color=dds[[INTGROUP]])
					# use show in new windows option since WebGl is not supported by Rstudio
					
					#library(rgl)
					#rglwidget() #This build of rgl does not include OpenGL functions.
				}
			}
		}
	}
	rm(index,dds,COUNT_NAME,PLOT,OUTNAME,vsd,NTOP)
}
#### ------- DEG-------- ####
# Note: to use apeglm, we have to set the reference level of the factor before the DESeq2 function, since it doe snot work with contrasts but only with coef
if(F){
	## Note: the reference should be the first level (because) to run apeglm coef
	
	register(MulticoreParam(14))
	
	# For the TUM runs
	#COEF="TEST_Tumor_CD4_vs_NormalBrain_CD4"
	#COEF="TEST_Tumor_CD8_vs_NormalBrain_CD8"
	#COEF="TEST_NormalBrain_CD4_vs_PeripheralBlood_CD4"
	#COEF="TEST_Tumor_CD4_vs_PeripheralBlood_CD4"
	#COEF="TEST_NormalBrain_CD8_vs_PeripheralBlood_CD8"
	#COEF="TEST_Tumor_CD8_vs_PeripheralBlood_CD8"
	COEF="TEST_Tumor_CD4_vs_Tumor_CD8"
	
	# For the BrainTIME
	#COEF="TEST_glioma_cd4_vs_healthyBlood_cd4"
	#COEF="TEST_glioma_cd8_vs_healthyBlood_cd8"
	#### COEF="TEST_glioma_cd4_vs_glioma_cd8"
	
	
	for( index in 1:2 ) {
		
		if(index==1) { dds=dds1 ; COUNT_NAME<-COUNT1_NAME} #"333TK_Human01"
		if(index==2) { dds=dds2 ; COUNT_NAME<-COUNT2_NAME} #"333TK_Human02"
		if(index==3) { dds=dds3 ; COUNT_NAME<-COUNT3_NAME} #"BrainTIME"
		
		#resultsNames(dds) # lists the coefficients 
		if(COEF %ni% resultsNames(dds) ) {print("Change the level of the factors to run apeglm") }
		#res <- results(dds,contrast=c("TEST","Tumor_CD4","NormalBrain_CD4"),alpha=0.05) #var, target, ref
		#res <- lfcShrink(dds,contrast=c("TEST","Tumor_CD4","NormalBrain_CD4"),type="normal",parallel=T) # shrink log fold changes
		res <- lfcShrink(dds,coef=COEF,type="apeglm",parallel=T) # shrink log fold changes
		res <- as.data.frame(res) %>% tibble::rownames_to_column("GeneID") %>% arrange(padj)
		res$signif <- res$padj<0.05
		OUTNAME=paste("DEG",COUNT_NAME,COEF,"xlsx",sep=".")
		write.xlsx(res,OUTNAME,overwrite=T) ; print(OUTNAME) ; print(getwd())
		rm(dds,res,OUTNAME)
	} ; rm(index,COEF)
	#https://www.molbiotools.com/listcompare.html
}
# using 'apeglm' for LFC shrinkage. If used in published research, please cite:
# Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
# sequence count data: removing the noise and preserving large differences.
# Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895



#### ------- Heatmaps (run1 + run2) - careful with batch effect-------- ####
if(F){
	#gene list - choice
	RUN=1 # the DEG will be from the run 1 tests, but the heatmap will have both run 1 + 2
	TOP=30
	
	# CD4
	if(T){
		FILE=paste0("03.DEG.NoRiboMito/CD4/T_vs_Blood/DEG.333TK_Human0",RUN,".TEST_Tumor_CD4_vs_PeripheralBlood_CD4.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD4" #& AnnotDf$Descr_Organ %in% c("Tumor","PeripheralBlood") 
	} # Tumor vs blood
	if(T){
		FILE=paste0("03.DEG.NoRiboMito/CD4/Norm_vs_Blood/DEG.333TK_Human0",RUN,".TEST_NormalBrain_CD4_vs_PeripheralBlood_CD4.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD4" #& AnnotDf$Descr_Organ %in% c("NormalBrain","PeripheralBlood") 
	} # Norm vs blood
	if(T){
		FILE=paste0("03.DEG.NoRiboMito/CD4/T_vs_Norm/DEG.333TK_Human0",RUN,".TEST_Tumor_CD4_vs_NormalBrain_CD4.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD4" #& AnnotDf$Descr_Organ %in% c("Tumor","NormalBrain") 
	} # Tumor vs normal brain
	
	#CD8
	if(T){
		FILE=paste0("03.DEG.NoRiboMito/CD8/T_vs_Blood/DEG.333TK_Human0",RUN,".TEST_Tumor_CD8_vs_PeripheralBlood_CD8.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD8" #& AnnotDf$Descr_Organ %in% c("Tumor","PeripheralBlood") 
	} # Tumor vs blood
	if(T){
		FILE=paste0("03.DEG.NoRiboMito/CD8/Norm_vs_Blood/DEG.333TK_Human0",RUN,".TEST_NormalBrain_CD8_vs_PeripheralBlood_CD8.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD8" #& AnnotDf$Descr_Organ %in% c("NormalBrain","PeripheralBlood") 
	} # Norm vs blood
	if(T){
		FILE=paste0("03.DEG.NoRiboMito/CD8/T_vs_Norm/DEG.333TK_Human0",RUN,".TEST_Tumor_CD8_vs_NormalBrain_CD8.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD8" #& AnnotDf$Descr_Organ %in% c("Tumor","NormalBrain") 
	} # Tumor vs normal brain
	
	# PLOT
	if(T){
		DATA<-read.xlsx( FILE )
		GENES<-DATA %>% filter(padj < 0.05) %>% 
			mutate( logFCtag = ifelse( log2FoldChange >0, "up", "down" ) ) %>%
			group_by(logFCtag) %>%	top_n(n = TOP,wt=dplyr::desc(padj)) %>% pull(GeneID)
		
		
		GENES_NOT_INTERSECT <- GENES [ GENES %ni% rownames(NC) ]
		print( paste( "genes not insersecting the run1 and run2:", length(GENES_NOT_INTERSECT)))
		print(paste0(GENES_NOT_INTERSECT,collapse=", "))
		GENES <- GENES[ GENES %in% rownames(NC) ]
		
		# if(length(GENES) > 100){
		# 	GENES_byVar <- GENES [ order( matrixStats::rowVars( as.matrix(NC[GENES,SELECTION]) ) , decreasing = T ) ][1:100]
		# 	GENES <- GENES_byVar ; rm(GENES_byVar)
		# 	TAG_genecount="_top100"
		# }
		
		OUTNAME=gsub(".xlsx$",paste0(TAG_Norm,"_top",TOP,"up+down","run1+2.DEG_FromRun1",".clustR",".pdf"),FILE) ; rm(FILE,DATA) #".clustR+C",
		pdf(OUTNAME,width=20,height=40)
		pheatmap::pheatmap( NC[ GENES, SELECTION], #cutree_cols=2,cutree_rows=2,
							cellwidth=10,cellheight=10,scale="row", #scale="none" "row",
							annotation_col = AnnotDf,cluster_rows=T,cluster_cols=F)
		dev.off() ; print(OUTNAME)
	}
}
#### ------- Heatmaps (run1 OR run2) (now +blood) -------- ####
if(T){
	#gene list - choice for Run 1 OR 2
	RUN=1
	TOP=30
	if(RUN==1){NC=NC1}
	if(RUN==2){NC=NC2}

	# CD4
	if(T){
		FILE=paste0("03.DEG.WithRiboMito/CD4/DEG.333TK_Human0",RUN,".TEST_Tumor_CD4_vs_PeripheralBlood_CD4.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD4" #& AnnotDf$Descr_Organ %in% c("Tumor","PeripheralBlood")

	} # Tumor vs blood
	if(T){
		FILE=paste0("03.DEG.WithRiboMito/CD4/DEG.333TK_Human0",RUN,".TEST_NormalBrain_CD4_vs_PeripheralBlood_CD4.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD4" #& AnnotDf$Descr_Organ %in% c("NormalBrain","PeripheralBlood") 
	} # Norm vs blood
	if(T){
		FILE=paste0("03.DEG.WithRiboMito/CD4/DEG.333TK_Human0",RUN,".TEST_Tumor_CD4_vs_NormalBrain_CD4.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD4" #& AnnotDf$Descr_Organ %in% c("Tumor","NormalBrain") 
	} # Tumor vs normal brain
	
	#CD8
	if(T){
		FILE=paste0("03.DEG.WithRiboMito/CD8/DEG.333TK_Human0",RUN,".TEST_Tumor_CD8_vs_PeripheralBlood_CD8.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD8" #& AnnotDf$Descr_Organ %in% c("Tumor","PeripheralBlood") 
	} # Tumor vs blood
	if(T){
		FILE=paste0("03.DEG.WithRiboMito/CD8/DEG.333TK_Human0",RUN,".TEST_NormalBrain_CD8_vs_PeripheralBlood_CD8.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD8" #& AnnotDf$Descr_Organ %in% c("NormalBrain","PeripheralBlood") 
	} # Norm vs blood
	if(T){
		FILE=paste0("03.DEG.WithRiboMito/CD8/DEG.333TK_Human0",RUN,".TEST_Tumor_CD8_vs_NormalBrain_CD8.xlsx")
		SELECTION = AnnotDf$Descr_Cell.Type == "CD8" #& AnnotDf$Descr_Organ %in% c("Tumor","NormalBrain") 
	} # Tumor vs normal brain
	
	#CD4_vs_CD8
	if(T){
		FILE="03.DEG.WithRiboMito/CD4_vs_CD8/DEG.333TK_Human01.TEST_Tumor_CD4_vs_Tumor_CD8.xlsx"
		SELECTION = AnnotDf$Descr_Cell.Type %in% c("CD4","CD8")
	}

	
	
	# PLOT - using NC the normalized matrix
	if(T){
		DATA<-read.xlsx( FILE )
		
		OUTNAME=gsub(".xlsx$",paste0(TAG_Norm,"_top",TOP,"up+down",".clustR",".pdf"),FILE)
		# Extra step of removing the ribosomal proteins: change the name of the pdf too
		if(T){
			DATA <- subset( DATA, GeneID %ni% grep("^RP[LS]",DATA$GeneID,value=T) )
			DATA <- subset( DATA, GeneID %ni% grep("^MT-",DATA$GeneID,value=T) )
			OUTNAME=gsub(".xlsx$",paste0(TAG_Norm,"_top",TOP,"up+down",".clustR","_NoRiboMito.pdf"),FILE)# ; rm(FILE,DATA) #"+blood", ".clustR+C"
		}
		GENES<-DATA %>% filter(padj < 0.05) %>% 
			mutate( logFCtag = ifelse( log2FoldChange >0, "up", "down" ) ) %>%
			group_by(logFCtag) %>%	top_n(n = TOP,wt = dplyr::desc(padj)) %>% pull(GeneID)
		SELECTION <- SELECTION [ AnnotDf$Batch == RUN ]
		
		pdf(OUTNAME,width=15,height=40)
		pheatmap::pheatmap( NC[ GENES, SELECTION],cutree_cols=3,cutree_rows=2,
							cellwidth=10,cellheight=10,scale="row", #scale="none",
							annotation_col = AnnotDf[,c("Descr_Cell.Type","Descr_Organ")], labels_col=AnnotDf[which(SELECTION),"UniqueSampleID"],
							cluster_rows=T,cluster_cols=F)
		dev.off() ; print(OUTNAME)
		
		rm(FILE,DATA)
	}
	#pdf_crop_auto -r
	
}
#### ------- GSEA : GCT expression + CLS pheno ####
# GSEA
if(F){
	#https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Using_RNA-seq_Datasets_with_GSEA
	
	# genes should be ranked for each condition?
	dds=dds1 ; COUNT_NAME=COUNT1_NAME
	
	dds=dds2 ; COUNT_NAME=COUNT2_NAME
	
	# now: -NoTechRepli
	
	DF = cbind(GeneID=rownames(dds), Description=NA, as.data.frame(counts(dds,normalized=TRUE)))
	GENES=rownames(dds)
	if(T){ 
		TAG_sel="WithRiboMito" 
	}else{
		GENES <- grep(pattern = "^RP[SL]",GENES,value=T,invert=T)
		GENES <- grep(pattern = "^MT-",GENES,value=T,invert=T)
		TAG_sel="NoRiboMito"
	}
	DF<-DF[GENES,]

	nbGenes=length(GENES)
	nbSamples=length(colnames(dds))
	nbClasses=length(unique(ANNOT$TEST))
	
	#1.2
	#nbGenes nbSamples
	FILENAME = paste0( c(COUNT_NAME,TAG_sel,"NormExpr","gct") ,collapse=".")
	write("#1.2",file=FILENAME,append=F,sep="\t")
	write(c( nbGenes , nbSamples),file=FILENAME,sep="\t",append=T)
	write.table(x=DF,file=FILENAME,sep="\t",quote=F,row.names=F,append=T)
	
	#nbSamples nbClasses 1
	FILENAME = paste0( c(COUNT_NAME,TAG_sel,"Pheno","cls") ,collapse=".")
	VALUES=c(nbSamples,nbClasses,"1")
	write(VALUES,file=FILENAME,append=F,sep="\t",ncolumns=length(VALUES))
	VALUES=c("#", as.character(unique(ANNOT$TEST) ) )
	write(VALUES,file=FILENAME,append=T,sep="\t",ncolumns=length(VALUES))
	VALUES=as.character(ANNOT$TEST)
	write(VALUES,file=FILENAME,append=T,sep="\t",ncolumns=length(VALUES))
	
	
	
}
#### session_info ####
writeLines(capture.output(devtools::session_info()), "sessionInfo.txt")



#### GSEA ####

# the following code select the leading edge genes with 
# FDR.q.val < 0.05
# NormalBrain_CD4		0
# Tumor_CD4				73
# NormalBrain_CD8		0
# Tumor_CD8				139

if(F){
	#Molecular Signatures Database (MSigDB) -> https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
	#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C7
	#nano /home/ga94rac/bin/gsea
	# setting -Xmx16g			#16Go RAM instead of 4 to avoid Java heap space error
	
	#### Making the list of the dataset to extract ####
	# find jul22/ -name "gsea_report_for*.tsv"
	# jul22/CD8_T_vs_N_Run1.Gsea.1626977111817/gsea_report_for_NormalBrain_CD8_1626977111817.tsv
	# jul22/CD8_T_vs_N_Run1.Gsea.1626977111817/gsea_report_for_Tumor_CD8_1626977111817.tsv
	# jul22/CD4_T_vs_N_Run1.Gsea.1626959579200/gsea_report_for_NormalBrain_CD4_1626959579200.tsv
	# jul22/CD4_T_vs_N_Run1.Gsea.1626959579200/gsea_report_for_Tumor_CD4_1626959579200.tsv
	
	LIST<-list.files("/home/ga94rac/gsea_home/output/jul28",pattern="gsea_report_for_.*tsv",full.names=T,recursive=T)
	
	FILES<-lapply(LIST, function(x){
		df<-read.table(x,header=T,sep="\t")
		df$source=x
		return(df)
		})
	names(FILES)<-gsub("_[^_]+.tsv","",gsub("gsea_report_for_","",basename(LIST)))
	FILES_Filter<-lapply(FILES, function(df){
		subset(df, df$FDR.q.val < 0.01)
	})
	lapply(FILES_Filter,nrow)
	
	# NormalBrain_CD4		0
	# Tumor_CD4				73
	# NormalBrain_CD8		0
	# Tumor_CD8				139
	FILES_NAMES<-lapply(FILES_Filter,function(x){ if(nrow(x)>0){ return( paste(dirname(x$source), paste0(x$NAME,".tsv"),sep="/") ) }else( return(NA)) } ) 
	# list of list
	# load all the signif reports
	SELECTED_FILES<-lapply(FILES_NAMES, function( group_name ){
		LIST<-lapply( group_name, function( file_name ){
			#print(paste0("x file= ", file_name))
			if( is.na( file_name ) ){print(paste0("skipped: ", file_name));return(NA)}
			#if( file_name == "/.tsv" ){print(paste0("skipped: ",file_name));return(NA)}
			#print(paste0("currently loading: ",file_name))
			df<-read.table(file_name,header=T,sep="\t")
			#df$source=file_name
			df$name= gsub(".tsv","",basename(file_name))
			return(df)
		})
		if(length(LIST)>1){
			names(LIST)<- lapply(LIST, function(x){ x$name[1]} )
		}
		return(LIST)
	})
	
	names(SELECTED_FILES)
	# the following are NA, no significance
	SELECTED_FILES$NormalBrain_CD4
	SELECTED_FILES$NormalBrain_CD8
	# the following are significant
	SELECTED_FILES$Tumor_CD4
	SELECTED_FILES$Tumor_CD8
	
	saveRDS(SELECTED_FILES,"SELECTED_FILES.Rds")

}



#### LEM ####
if(F){
	#A ) Selection of leading edge genes
	#Those are in the files with ".tsv" extension, as CORE_ENRICHMENT == "Yes"
	
	SELECTED_FILES<-readRDS("SELECTED_FILES.Rds")
	
	#we have to ignore the empty data to avoid error (element 1 and 3 of the list: NormalBrain CD4 and CD8)
	#select for the core enriched
	SELECTED_GENES<-lapply(SELECTED_FILES[c(2,4)], function(dfs){
		lapply(dfs, function(df){
			return( df[ df$CORE.ENRICHMENT == "Yes" , ] )
		})
	})
	
	
	
	names(SELECTED_FILES)
	FINAL_CD4 <- do.call(rbind, SELECTED_FILES[[2]] )
	FINAL_CD8 <- do.call(rbind, SELECTED_FILES[[4]] )
	
	#save this data to avoid recomputing
	#saveRDS(SELECTED_GENES,"SELECTED_GENES.Rds")
	#SELECTED_GENES<-readRDS("SELECTED_GENES.Rds")
	
	#save this to check manually, eventually
	write.xlsx( FINAL_CD4 , "FINAL.Tumor_CD4.xlsx")
	write.xlsx( FINAL_CD8 , "FINAL.Tumor_CD8.xlsx")
	#write.xlsx( data.frame(GeneSymbol=unique(FINAL_CD4[ FINAL_CD4$CORE.ENRICHMENT == "Yes" , "SYMBOL" ])) , "FINAL_UniqLeadingEdge.Tumor_CD4.xlsx")
	#write.xlsx( data.frame(GeneSymbol=unique(FINAL_CD8[ FINAL_CD8$CORE.ENRICHMENT == "Yes" , "SYMBOL" ])) , "FINAL_UniqLeadingEdge.Tumor_CD8.xlsx")
	
	#clean
	rm(SELECTED_FILES,SELECTED_GENES)
	
	# load the genes for the Matrix M (unique genes from leading edge of significants gene sets.)
	# We need the shrunken values if the gene is in the leading edge, and zero if not in leading edge
	
	# CD4
	LeadingEdge_CD4 = unique(FINAL_CD4[ FINAL_CD4$CORE.ENRICHMENT == "Yes" , "SYMBOL" ])
	AllGenesInSets_CD4=FINAL_CD4[ , "SYMBOL" ] #should not be uniq since we will find groups?
	#get shrunken data for the condition
	DEG_CD4 <- read.xlsx("03.DEG.WithRiboMito/CD4/DEG.333TK_Human01.TEST_Tumor_CD4_vs_NormalBrain_CD4.xlsx")
	rownames(DEG_CD4)<-DEG_CD4$GeneID
	# get the values for the matrix
	CD4 <- DEG_CD4[ AllGenesInSets_CD4, "log2FoldChange" ]
	names(CD4)<-AllGenesInSets_CD4
	CD4[ AllGenesInSets_CD4 %ni% LeadingEdge_CD4 ] <- 0
	#make matrix
	Mat_CD4 <- matrix(data=rep( CD4,length(LeadingEdge_CD4) ),byrow=T,
					  nrow=length(LeadingEdge_CD4),
					  ncol=length(AllGenesInSets_CD4) )
	rownames(Mat_CD4)<-paste0("LeadEdge_",LeadingEdge_CD4)
	colnames(Mat_CD4)<-AllGenesInSets_CD4
	#Mat_CD4[1:5,1:15]
	#AllGenesInSets_CD4 [ AllGenesInSets_CD4 %ni% LeadingEdge_CD4 ] #"PNKP" not leading edge ; "MELK" LE
	#Mat_CD4[ 1:5, c("MELK","PNKP")]
	rm(LeadingEdge_CD4,AllGenesInSets_CD4,DEG_CD4,CD4)
	
	
	# CD8
	LeadingEdge_CD8 = unique(FINAL_CD8[ FINAL_CD8$CORE.ENRICHMENT == "Yes" , "SYMBOL" ])
	AllGenesInSets_CD8=FINAL_CD8[ , "SYMBOL" ] #should not be uniq since we will find groups?
	#get shrunken data for the condition
	DEG_CD8 <- read.xlsx("03.DEG.WithRiboMito/CD8/DEG.333TK_Human01.TEST_Tumor_CD8_vs_NormalBrain_CD8.xlsx")
	rownames(DEG_CD8)<-DEG_CD8$GeneID
	# get the values for the matrix
	CD8 <- DEG_CD8[ AllGenesInSets_CD8, "log2FoldChange" ]
	names(CD8)<-AllGenesInSets_CD8
	CD8[ AllGenesInSets_CD8 %ni% LeadingEdge_CD8 ] <- 0
	#make matrix
	Mat_CD8 <- matrix(data=rep( CD8,length(LeadingEdge_CD8) ),byrow=T,
					  nrow=length(LeadingEdge_CD8),
					  ncol=length(AllGenesInSets_CD8) )
	rownames(Mat_CD8)<-paste0("LeadEdge_",LeadingEdge_CD8)
	colnames(Mat_CD8)<-AllGenesInSets_CD8
	Mat_CD8[1:5,1:5]
	rm(LeadingEdge_CD8,AllGenesInSets_CD8,DEG_CD8,CD8)
	
	#B) "non-negative matrix factorization" using fold change (actually a grouping of genes by membership)
	
	#https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
	#install.packages("NMF")
	library(NMF)
	nmfAlgorithm(version='R')
	meth <- "nsNMF"
	
	data(esGolub)
	esGolub
	esGolub <- esGolub[1:200,]
	
	#nmf.getOption("default.seed")
	res <- nmf(x = Mat_CD4, rank = 1, method = meth, seed=123456) # ?? specification of the factorization rank. 
	
	#get the GENES vs GENES matrix, using shrunken log Fold Change
	
	object <- nmfObject(object)
	
}


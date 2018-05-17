library(ctc)
library(heatmap.plus)
library(GSVA)

tumor_subtyping <- function(indication=indication,
							inFile=inFile,
							datalib_dir=datalib_dir,
							prefix=NA){
	
	centroid_clustering_indications<-c('LIHC','LUSC','LUAD','SKCM','PRAD')

	
	if (indication=="BRCA"){
		##----------------------------------------------------------------------------------------------
		print (paste0("Begin ",indication," samples subtyping...\n"))
		
		## Perform PAM50 classification on BRCA samples
		paramDir<- paste(datalib_dir,"datalib/breast_cancer/PAM50_classification/bioclassifier_R",sep="/") # the location of unchanging files such as the function library and main program
		inFileElements<- unlist(strsplit(inFile,"/"))
		inputDir<- paste(inFileElements[1:length(inFileElements)-1],collapse="/")  # the location of the data matrix, and where output will be located

		inputFile<- inFileElements[length(inFileElements)] # the input data matrix as a tab delimited text file
		
		if (is.na(prefix)){
			short<-indication # short name that will be used for output files
		}else{
			short<-paste0(prefix,"_",indication)
		}
		
		calibrationParameters<- NA 	#the column of the "mediansPerDataset.txt" file to use for calibration; 
									#NA will force centering within the test set & -1 will not do any 
									#adjustment (when adjustment performed by used)

		hasClinical<-FALSE 	#may include tumor size as second row, with 'T' as the gene name, 
							#and encoded as binary (0 for size <= 2cm or 1 for size > 2cm)
							#set this variable to FALSE if tumor size is not available

		collapseMethod<-"mean" # can be mean or iqr (probe with max iqr is selected)
							   # typically, mean is preferred for long oligo and
							   # iqr is preferred for short oligo platforms


		####
		# run the assignment algorithm
		####
		
		source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
		source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))
		
		PAM50_subtypePrediction(paramDir=paramDir,
									short=short,
									calibrationParameters=NA,
									stdArray=T,
									hasClinical=FALSE,
									collapseMethod="mean",
									inputDir=inputDir,
									inputFile=inputFile)
									
		print (paste0("Finish ",indication," samples subtyping.\n"))
		##----------------------------------------------------------------------------------------------
	}
	
	else if (indication=="COAD"||indication=="READ"||indication=="COAD&READ"){
		##----------------------------------------------------------------------------------------------
		print (paste0("Begin ",indication," samples subtyping...\n"))
		
		## Perform CMS classification on COAD&READ samples
		library(Biobase)
		library(CMSclassifier)

		## extract CMS gene list
		cms_genes<-read.table(paste(datalib_dir,"datalib/colorectal_cancer/CMSclassify_EntrezID_to_Genesymbol.txt",sep="/"),header=T,sep="\t",check.names=F)

		exp_tumor<-read.table(inFile,
							sep="\t",header = TRUE,row.names = 1,check.names=FALSE)

		## Orgnize the data for the input for CMS classifier
		overlap_genes<-intersect(cms_genes$Gene_symbol,rownames(exp_tumor))
		print (length(overlap_genes))

		id<-match(overlap_genes,rownames(exp_tumor))
		cms_exp_tumor<-exp_tumor[id,]
		id<-match(overlap_genes,cms_genes$Gene_symbol)
		rownames(cms_exp_tumor)<-cms_genes$EntrezID[id]

		cms_fn<-gsub(".txt","_CMS.txt",inFile)
		cms_fn<-gsub(".tsv","_CMS.tsv",cms_fn)
		write.table(cbind(feature=rownames(cms_exp_tumor),cms_exp_tumor),cms_fn,row.names=F,quote=F,sep="\t")

		crc_dist<-function(vec){
		  x<-as.numeric(vec[1:4])
		  i<-which(x==max(x))
		  y<-x[-i]
		  j<-which(x==max(y))
		  max_x<-x[i]
		  second_max_x<-x[j]
		  diff<-max_x-second_max_x
		  if (diff>=0.1){
			paste(c("CMS1","CMS2","CMS3","CMS4")[i],collapse=",")
		  }else{
			paste(vec[6])
		  }
		}

		crc_exp<-read.table(cms_fn,sep="\t",header = TRUE,row.names = 1,check.names=FALSE)
		
		## Perform random forest model based classification
		Rfcms <- CMSclassifier::classifyCMS(crc_exp,method="RF")[[3]]
		Rfcms$`curated.RF.predictedCMS`<-apply(Rfcms,1,crc_dist)
		
		inFileElements<- unlist(strsplit(inFile,"/"))
		inputDir<- paste(inFileElements[1:length(inFileElements)-1],collapse="/")
		
		if (is.na(prefix)){
			predFile<- paste0(inputDir,"/",indication,"_CMS_RF_prediction.txt")
		}else{
			predFile<- paste0(inputDir,"/",prefix,"_",indication,"_CMS_RF_prediction.txt")
		}
		write.table(cbind(sample_ID=rownames(Rfcms),Rfcms),predFile,row.names=F,quote=F,sep="\t")
		
		print (paste0("Finish ",indication," samples subtyping.\n"))
		##----------------------------------------------------------------------------------------------
	}
	
	else if (indication %in% centroid_clustering_indications){
		##----------------------------------------------------------------------------------------------
		print (paste0("Begin ",indication," samples subtyping...\n"))
		
		## Perform centroid clustering to infer subtypes on indications that have centroid matrix available
		
		inFileElements<- unlist(strsplit(inFile,"/"))
		inputDir<- paste(inFileElements[1:length(inFileElements)-1],collapse="/")
		
		if (is.na(prefix)){
			short<-paste0(inputDir,"/",indication,"_centroid") # short name that will be used for output files
		}else{
			short<-paste0(inputDir,"/",prefix,"_",indication,"_centroid")
		}

		if (indication=="LUAD"){
	      trainCentroids=paste0(datalib_dir,"/datalib/centroids/LAD_predictor_centroids_from_wilkerson.2012.txt")
		}else{
		  trainCentroids=paste0(datalib_dir,"/datalib/centroids/",indication,"_subtype_centroids_from_TCGA-8-11-2017.txt")
		}
		centroid_cluster(predFile=inFile,trainCentroids=trainCentroids,prefix=short,plot=F)
		
		print (paste0("Finish ",indication," samples subtyping.\n"))
		##----------------------------------------------------------------------------------------------
	}
	
	else if (indication=="STAD"){
		##----------------------------------------------------------------------------------------------
		print (paste0("Begin ",indication," samples subtyping...\n"))
		
		## Perform GSVA based subtype estimation on STAD samples
		## extract signature gene list
		sig_genes <- readLines(paste(datalib_dir,"datalib/gastric_cancer/gastric_cancer_gene_signature.txt",sep="/"))
		sig_cluster<-list()
		for (i in 1:length(sig_genes)){
		  data<-unlist(strsplit(sig_genes[i],"\t"))
		  cluster<-data[1]
		  genes<-data[2:length(data)]
		  sig_cluster[[cluster]]<-genes
		}

		stad_exp<-read.table(inFile,
							 sep="\t",header = TRUE,check.names=FALSE)
		rownames(stad_exp)<-stad_exp[,1]
		exp<-as.matrix(stad_exp[,-1])

		## gsva with bootstrap
		gsva.ori<-gsva(exp,sig_cluster,method="gsva",rnaseq=FALSE)
		obs<-t(gsva.ori$es.obs)

		## perform classification of STAD samples
		class<-vector()
		emt_label = which(colnames(obs)=="MSS_EMT")
		msi_label = which(colnames(obs)=="MSI")
		p53_label = which(colnames(obs)=="MSS_TP53")
		pro_label = which(colnames(obs)=="Proliferation")
		mlh1_label = which(colnames(obs)=="MLH1-mRNA")
		cdh1_label = which(colnames(obs)=="CDH1-mRNA")

		for (i in 1:nrow(obs)){
		  if (obs[i,msi_label]>0 && obs[i,emt_label]<=0){
			class<-append(class,"MSI")
		  }
		  else if (obs[i,msi_label]<=0 && obs[i,emt_label]>0){
			class<-append(class,"MSS_EMT")
		  }
		  else{
			if (obs[i,p53_label]>0){
			  class<-append(class,"MSS_TP53+")
			}
			else{
			  class<-append(class,"MSS_TP53-")
			}
		  }
		}
		
		inFileElements<- unlist(strsplit(inFile,"/"))
		inputDir<- paste(inFileElements[1:length(inFileElements)-1],collapse="/")
		
		if (is.na(prefix)){
			predFile<- paste0(inputDir,"/",indication,"_gsva_classification.txt")
		}else{
			predFile<- paste0(inputDir,"/",prefix,"_",indication,"_gsva_classification.txt")
		}
		write.table(cbind(sampleID=rownames(obs),obs,Subtype=class),
					predFile,
					sep="\t",row.names=F,quote=F)
					
		print (paste0("Finish ",indication," samples subtyping.\n"))
		##----------------------------------------------------------------------------------------------
	}
	
	else if (indication=="OV"){
		##----------------------------------------------------------------------------------------------
		print (paste0("Begin ",indication," samples subtyping...\n"))
		
		## Perform GSVA based subtype estimation on OV samples
		## extract ovarian cancer subtype gene list
		subtype_genes <- readLines(paste(datalib_dir,"datalib/ovarian_cancer/OV-TP.mRNAarray.subclassmarkers.p0.01.fc0.5.txt",sep="/"))
		subtype_cluster<-list()
		for (i in 2:length(subtype_genes)){
		  data<-unlist(strsplit(subtype_genes[i],"\t"))
		  cluster<-data[1]
		  gene<-data[2]
		  if (!cluster %in% names(subtype_cluster)){
			subtype_cluster[[cluster]]<-vector()
			subtype_cluster[[cluster]]<-append(subtype_cluster[[cluster]],gene)
		  }else{
			subtype_cluster[[cluster]]<-append(subtype_cluster[[cluster]],gene)
		  }
		}
		##----------------------------------------------------------------------------------------------

		subtypes <- list('1'='Proliferative',
						 '2'='Mesenchymal',
						 '3'='Differentiated',
						 '4'='ImmuneReactive')

		ov_exp<-read.table(inFile,
						   sep="\t",header = TRUE,check.names=FALSE)
		rownames(ov_exp)<-ov_exp[,1]
		exp<-as.matrix(ov_exp[,-1])

		## gsva with bootstrap
		gsva.ori<-gsva(exp,subtype_cluster,method="gsva",rnaseq=FALSE)
		obs<-t(gsva.ori$es.obs)
		class<-apply(obs,1,function(x){i=(which(x==max(x)));paste(subtypes[[names(subtype_cluster)[i]]])})
		
		inFileElements<- unlist(strsplit(inFile,"/"))
		inputDir<- paste(inFileElements[1:length(inFileElements)-1],collapse="/")
		
		if (is.na(prefix)){
			predFile<- paste0(inputDir,"/",indication,"_gsva_classification.txt")
		}else{
			predFile<- paste0(inputDir,"/",prefix,"_",indication,"_gsva_classification.txt")
		}
		write.table(cbind(sampleID=rownames(obs),obs,Subtype=class),
					predFile,
					sep="\t",row.names=F,quote=F)
					
		print (paste0("Finish ",indication," samples subtyping.\n"))
		##----------------------------------------------------------------------------------------------
	}
	
	else{
		quit(paste0("Sorry, the subtyping of ",indication," is not supported at this moment!\n"))
		
	}

}


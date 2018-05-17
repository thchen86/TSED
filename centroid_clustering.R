
##---------------------------------------------------------------
## This is used to create universal/subtype-specific centroid matrix

## INPUT:
##     1) data_mat: the quantile-normalized TPM value matrix, the row is genes and the column is sample names
##     2) gene_vec: a vector of interesting genes to be used for the centroid clustering 
##     3) subtype_mat: a matrix with sample names and associated subtype information, the first col has to be the sample name, 
##        the second col has to be the subtype labels, afterward cols do not matter to this analysis
## OUTPUT:
##     1) subtype_centroids: a centroid matrix contains the gene median centroid for each subtype, the row is genes
##        and the column is subtypes
##     2) universal_centroids: a centroid matrix contains the gene median centroid for the whole training dataset, the row
##        is genes and the column is specific data types

create_centroid<-function(data_mat,
                          gene_vec,
                          subtype_mat){
  
  res<-list()
  
  overlap_genes<-intersect(gene_vec,rownames(data_mat))
  gene_vec<-gene_vec[match(overlap_genes,gene_vec)]
  data_mat<-data_mat[match(overlap_genes,rownames(data_mat)),]
  
  ## create a centroid vector based on the data_mat and serve as a reference centroid for future analysis
  overall_centroid_matrix<-matrix(nrow=length(gene_vec),ncol=0)
  rownames(overall_centroid_matrix)<-gene_vec
  
  overall_centroid_matrix<-cbind(overall_centroid_matrix,apply(data_mat,1,median))
  colnames(overall_centroid_matrix)<-c("log2.quantNorm.TPM.centroid")
  
  ## create the centroid matrix that contains the centroid matrix for each subtype
  ctr_data_mat<-medianCtr(data_mat)
  uniq_labels<-sort(unique(subtype_mat[,2]))
  
  subtype_centroid_matrix<-matrix(nrow=length(gene_vec),ncol=length(uniq_labels))
  rownames(subtype_centroid_matrix)<-gene_vec
  colnames(subtype_centroid_matrix)<-uniq_labels
  
  for (label in uniq_labels){
    label_index<-which(subtype_mat[,2]==label)
    samples<-subtype_mat[label_index,1]
    
    overlap_samples<-intersect(samples,colnames(ctr_data_mat))
    sample_index<-match(overlap_samples,colnames(ctr_data_mat))
    
    gene_centroids<-apply(ctr_data_mat[,sample_index],1,mean)
    
    subtype_centroid_matrix[,which(colnames(subtype_centroid_matrix)==label)]<-gene_centroids
  }
  
  res$subtype_centroids<-subtype_centroid_matrix
  res$universal_centroids<-overall_centroid_matrix
  
  return (res)
}
##---------------------------------------------------------------

##---------------------------------------------------------------
#1. Download normalized expression data
#2. Log transform expression estimates
#3. Optionally median center or appropriately adjust each probeset, which carries the population assumption
#4. Map probesets to Entrez gene names
#5. Format as tab delimited text, with first row of sample names and first column of gene names
#---- This software then provides the following steps
#5. For probesets that map to identical Entrez gene names, select the one with highest IQR (for Affy, select mean for Agilent)
#6. Extract the genes of interest
#7. Calculate Spearman's rank correlation between each sample and each subtype centroid (in pam50_centroids.txt)
#8. Assign the class of the most highly correlated centroid to each sample

## INPUT:
#     1) trainCentroid: the centroid matrix file, the rows are genes and the cols are the subtypes
#     2) predFile: the data matrix that will be used to predict subtypes
#     3) stdArray: just for visualization, and only set to F if many missing genes
#     4) prefix: the prefix-name of the output files
#     5) calibrationParameters<- NA 	#the column of the calibrationFile to use for calibration; 
#        NA will force centering within the test set & -1 will not do any adjustment (when adjustment performed by used)
#     6) calibrationFile: supplied if and only if you set calibrationParameters!=NA|-1 
#     7) collapseMethod<-"mean" # can be mean or iqr (probe with max iqr is selected)
#        typically, mean is preferred for long oligo and iqr is preferred for short oligo platforms
#     8) plot: used if want to generate the plots

centroid_cluster<-function(trainCentroids,
                           predFile,
                           stdArray=T,
                           prefix,
                           calibrationParameters=NA,
                           calibrationFile=NA,
                           collapseMethod="mean",
                           plot=F){
  
  # load the centroids for classifcation
  pamout.centroids<-read.table(trainCentroids,sep="\t",header=T,row.names=1,check.names=F)
  
  pdfname1<-paste(prefix,"_predictionScores_RankCorrelation_1.pdf",sep="")
  clustername<-paste(prefix,"_normalized_heatmap",sep="")
  outFile<- paste(prefix,"_classification_scores.txt",sep="")
  
  y<-readarray(predFile,hr=1,method=collapseMethod,impute=F)
  
  # normalization
  if(is.na(calibrationParameters)){
    y$xd<-medianCtr(y$xd)
  }else{
    if(calibrationParameters != -1){
      medians<-readarray(calibrationFile,hr=1)
      print(paste("calibration to:",dimnames(medians$xd)[[2]][calibrationParameters]))
      tm<-overlapSets(medians$xd,y$xd)
      y$xd<-(tm$y-tm$x[,calibrationParameters])
      #y$xd<-(tm$y-tm$x[,calibrationParameters])/tm$x[,15]
    }
  }
  
  num.missing<- NA
  
  if(stdArray){
    y$xd<-standardize(y$xd)
  }
  
  # assign the subtype scores
  out<-sspPredict(pamout.centroids,classes="",y$xd,std=F,distm="spearman",centroids=T)
  out$distances<- -1*out$distances
  
  call.conf<-c()
  for(j in 1:length(out$predictions)){
    call.conf[j]<- 1-cor.test(out$testData[,j],out$centroids[,which(colnames(pamout.centroids)==out$predictions[j])],method="spearman")$p.value
  }
  call.conf<-round(call.conf,2)
  
  # write output files
  outtable<-cbind(out$distances, out$predictions, call.conf)
  dimnames(outtable)[[2]]<-c(colnames(pamout.centroids),"Call","Confidence")
  
  write.table(outtable,outFile,sep="\t",col.names=NA,quote=F)
  
  if (plot){
    # make some plots for evaluation
    colors<-c("red","hotpink","darkblue","cyan","green","yellow","brown","magenta")
    subtypeColors<-out$predictions
    for (i in 1:length(colnames(pamout.centroids))){
      subtype<-colnames(pamout.centroids)[i]
      subtypeColors[subtypeColors==subtype]<-colors[i]
    }
    conf.colors<-call.conf
    conf.colors[call.conf>=0.95]<-"black"
    conf.colors[call.conf<0.95]<-"red"
    
    pdf(paste(clustername,".pdf",sep=""))
    myHeatmap(out$testData,cbind(subtypeColors,conf.colors),file=paste(clustername,".cdt",sep=""),rowNames=rownames(out$testData))
    dev.off()
    
    pdf(pdfname1,height=10,width=12)
    pars<-par(no.readonly=T)
    myplot(out,prefix,colnames(pamout.centroids))
    dev.off()
  }
}
##---------------------------------------------------------------

##---------------------------------------------------------------
## Below are centroid clustering related functions

# function for boxplots of correlation by subtype
myplot<-function(y,prefix,subtypes){
  par(mfrow=c(round(length(subtypes)/2+0.1),2),mar=c(5,3,2,2),las=3)
  y$prediction<-factor(y$prediction,levels=subtypes)
  for (i in 1:length(subtypes)){
    subtype<-subtypes[i]
    boxplot(y$distances[,i]~y$prediction,border=8,ylab=paste(subtype," Correlation",sep=""),main=subtype)
    stripchart(y$distances[,i]~y$prediction,vertical=T,method="jitter",pch=3,add=T)
  }
}

medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}

pcaEA<-function(x,classes,size=1,showLegend=T,legendloc="topright",mainStr="",startPC=1,stopPC=2,showNames=T,showClasses=F,axisExpansion=0,groupColors=NA){
  
  features<- dim(x)[1]
  samples<- dim(x)[2]
  sampleNames<- dimnames(x)[[2]]
  featureNames<-dimnames(x)[[1]]
  x<-apply(x,2,as.numeric)
  
  #principal components plots
  data.pca<-prcomp(as.matrix(x))
  
  # Proportion of total variance distributed over 10 first components:
  tmp<-data.pca$sdev[1:10]^2/sum(data.pca$sdev^2)
  
  gr.labels<-as.vector(t(classes))
  gr.labels.fac<-factor(as.vector(t(classes)),exclude="")
  nlabels<-nlevels(gr.labels.fac)
  legendLabels<-vector()
  legendColors<-vector()
  for(k in 1:nlabels){
    group<-levels(gr.labels.fac)[k]
    legendLabels[k]<-group
    if(length(groupColors)>1){
      gr.labels[gr.labels.fac==group]<-groupColors[k]
      legendColors[k]<-groupColors[k]
    }else{
      gr.labels[gr.labels.fac==group]<-k
      legendColors[k]<-k
    }
  }
  
  if(length(groupColors)==1){
    gr.labels<-as.numeric(gr.labels)
  }
  
  #plot 2pcs by each other
  i<-startPC
  j<-stopPC
  
  #graphing parameters
  par(lab=c(3,4,3))
  par(mgp=c(.3,.5,.0))
  par(mai=c(.5,.5,.5,.5))
  par(xaxt="n",yaxt="n")
  
  strM<-mainStr
  strX<-paste("PC",i,paste("(",round(tmp[i],4)*100,"%)",sep=""),sep=" ")
  strY<-paste("PC",j,paste("(",round(tmp[j],4)*100,"%)",sep=""),sep=" ")
  xmin<-min(data.pca$rotation[,i])-abs(axisExpansion*min(data.pca$rotation[,i]))
  xmax<-max(data.pca$rotation[,i])+abs(axisExpansion*max(data.pca$rotation[,i]))
  ymin<-min(data.pca$rotation[,j])-abs(axisExpansion*min(data.pca$rotation[,j]))
  ymax<-max(data.pca$rotation[,j])+abs(axisExpansion*max(data.pca$rotation[,j]))
  plot(data.pca$rotation[,i],data.pca$rotation[,j], xlab=strX, ylab=strY, 
       main=strM, col=gr.labels,cex=size,pch="",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
  if(showNames){
    text(data.pca$rotation[,i],data.pca$rotation[,j],labels=names(data.pca$rotation[,i]),cex=size*.6)
  }else{
    if(showClasses){
      text(data.pca$rotation[,i],data.pca$rotation[,j],labels=gr.labels.fac,cex=size*.6)
    }else{
      points(data.pca$rotation[,i],data.pca$rotation[,j],col=gr.labels,cex=size*1.5,pch=19)
    }
    if(showLegend){
      legend(legendloc,legend=legendLabels,col=legendColors,pch=19,x.intersp=.3,yjust=.5,bty="n",cex=size)
    }
  }
}

cols <- function(lowi = "green", highi = "red", ncolors = 20) {
  low <- col2rgb(lowi)/255
  high <- col2rgb("black")/255
  col1 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                       high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  low <- col2rgb("black")/255
  high <- col2rgb(highi)/255
  col2 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                       high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  col<-c(col1[1:(ncolors-1)],col2)
  return(col)
}

myHeatmap<-function(x,t.colors=NA,fileName="cluster.cdt",linkage="complete",distance="pearson",contrast=2,returnSampleClust=F,rowNames=NA,rightMar=7,bottomMar=1,colNames=NA){
  
  temp<-hclust2treeview(x,method=distance,file=fileName,link=linkage,keep.hclust=T)
  gTree<-temp[[1]]
  sTree<-temp[[2]]
  
  imageVals<-x
  imageVals[x > contrast] <- contrast
  imageVals[x < -1 * contrast] <- -1 * contrast
  
  if(sum(is.na(t.colors))>0){
    heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
            col=cols(),labCol=colNames, scale="none",
            margins=c(bottomMar,rightMar),labRow=rowNames)
  }else{
    if(length(t.colors)>dim(imageVals)[2]){
      heatmap.plus(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
                   col=cols(),labCol=colNames,labRow=rowNames,scale="none",
                   ColSideColors=t.colors, margins=c(bottomMar,rightMar))
    }else{
      heatmap(imageVals,Rowv=as.dendrogram(gTree),Colv=as.dendrogram(sTree),
              col=cols(),labCol=colNames,labRow=rowNames,scale="none",
              ColSideColors=as.vector(t(t.colors)), margins=c(bottomMar,rightMar))
    }
  }
  if(returnSampleClust){
    return(sTree)
  }
}


standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}


overlapSets<-function(x,y){
  
  # subset the two lists to have a commonly ordered gene list
  x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
  y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]
  
  #and sort such that thing are in the correct order
  x<-x[sort.list(row.names(x)),]
  y<-y[sort.list(row.names(y)),]
  
  return(list(x=x,y=y))
}


collapseIDs<-function(x,allids=row.names(x),method="mean"){
  
  allids<-as.vector(allids)
  ids<- levels(as.factor(allids))
  x.col<- NULL
  
  if(length(ids)==dim(x)[1]){ 
    dimnames(x)[[1]]<-allids
    return(x) 
  }
  
  for(i in 1:length(ids)){
    if(sum(allids==ids[i])>1){
      indices <- allids==ids[i] 
      if(method=="mean"){
        vals<-apply(x[indices,],2,mean,na.rm=T)
      }
      if(method=="median"){
        vals<-apply(x[indices,],2,median,na.rm=T)
      }
      if(method=="stdev"){   
        temp<- x[indices,]
        stdevs<- apply(temp,1,sd,na.rm=T)
        vals<- temp[match(max(stdevs),stdevs),]
      }
      if(method=="iqr"){   
        temp<- x[indices,]
        iqrs<- apply(temp,1,function(x){quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T)})
        vals<- temp[match(max(iqrs),iqrs),]
      }
      x.col <- rbind(x.col,vals)
    }else{
      x.col <- rbind(x.col,x[allids==ids[i],])
    }
  }
  
  dimnames(x.col)<- list(ids,dimnames(x)[[2]])
  return(x.col)
  
}


readarray<-function(dataFile,designFile=NA,hr=1,impute=T,method="mean"){
  
  headerRows <- hr
  
  x<-read.table(dataFile,sep="\t",header=F,fill=T,stringsAsFactors=FALSE)
  
  if(headerRows==1){
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]
    classes<-NULL
    ids<-x[,1]
    xd<-x[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)	
  }else{
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]
    
    classes<-x[1:(headerRows-1),]
    dimnames(classes)[[1]]<-classes[,1]
    classes<-classes[,-1]
    classes[classes==""]<-NA
    classes<-t(classes)
    rownames(classes)<-sampleNames
    classes<-as.data.frame(classes)
    
    xd<-x[(-1:-(headerRows-1)),]
    ids<-as.vector(t(xd[,1]))
    xd<-xd[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)
  }
  
  features<- dim(xd)[1]
  samples<- dim(xd)[2]
  geneNames<-rownames(xd)
  xd<-apply(xd,2,as.numeric)
  rownames(xd)<-geneNames
  colnames(xd)<-sampleNames
  
  if(!is.na(designFile)){
    x<-read.table(designFile,sep="\t",header=T,row.names=1,fill=T,,stringsAsFactors=FALSE)
    xd<-xd[,sort.list(colnames(xd))]
    xd<-xd[,colnames(xd) %in% rownames(x)]
    x<-x[rownames(x) %in% colnames(xd),]
    x<-x[sort.list(rownames(x)),]
    classes<-as.data.frame(x)
  }
  
  if(sum(apply(xd,2,is.na))>0 & impute){
    library(impute)
    allAnn<-dimnames(xd)
    data.imputed<-impute.knn(as.matrix(xd))$data
    xd<-data.imputed[1:features,]
    dimnames(xd)<-allAnn
  }
  
  return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}


sspPredict<-function(x,classes="",y,nGenes="",priors="equal",std=F,distm="euclidean",centroids=F){
  
  dataMatrix<-x
  features<- dim(x)[1]
  samples<- dim(x)[2]
  sampleNames<- dimnames(x)[[2]]
  featureNames<- dimnames(x)[[1]]
  
  #parse the test file - same as train file but no rows of classes
  tdataMatrix<-y
  tfeatures<- dim(y)[1]
  tsamples<- dim(y)[2]
  tsampleNames<- dimnames(y)[[2]]
  tfeatureNames<- dimnames(y)[[1]]
  
  #dimnames(tdataMatrix)[[2]]<-paste("x",seq(1,471))
  temp <- overlapSets(dataMatrix,tdataMatrix)
  dataMatrix <- temp$x
  tdataMatrix <- temp$y
  sfeatureNames<-row.names(dataMatrix)
  
  # standardize both sets
  if(std){
    dataMatrix<-standardize(dataMatrix)
    tdataMatrix<-standardize(tdataMatrix)
  }
  
  if(!centroids){
    thisClass <- as.vector(classes[,1])
    nClasses<-nlevels(as.factor(thisClass))
    classLevels<-levels(as.factor(thisClass))
    for(j in 1:nClasses){
      thisClass[thisClass==classLevels[j]] <- j
    }
    thisClass<-as.numeric(thisClass)
    dataMatrix <- dataMatrix[,!(is.na(thisClass))]
    thisClass <- thisClass[!(is.na(thisClass))]
    
    scores<-apply(dataMatrix,1,bwss,thisClass)
    trainscores<-vector()	
    for(j in 1:dim(dataMatrix)[1]){			
      trainscores[j]<-scores[[row.names(dataMatrix)[j]]]$bss / scores[[row.names(dataMatrix)[j]]]$wss
    }
    
    dataMatrix<-dataMatrix[sort.list(trainscores,decreasing=T),]
    tdataMatrix<-tdataMatrix[sort.list(trainscores,decreasing=T),]	
    
    if(nGenes==""){
      nGenes<-dim(dataMatrix)[1]
    }
    print(paste("Number of genes used:",nGenes))
    
    dataMatrix<-dataMatrix[1:nGenes,]
    tdataMatrix<-tdataMatrix[1:nGenes,]
    
    centroids<-matrix(,nrow=nGenes,ncol=nClasses)
    for(j in 1:nClasses){
      centroids[,j]<-apply(dataMatrix[,thisClass==j],1,mean)
    }
    dimnames(centroids)<-list(row.names(dataMatrix),NULL)
    
  }else{
    nGenes<-dim(dataMatrix)[1]
    print(paste("Number of genes used:",nGenes))
    centroids<-dataMatrix
    nClasses<-dim(centroids)[2]
    classLevels<-dimnames(centroids)[[2]]
  }
  
  distances<-matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
  for(j in 1:nClasses){
    if(distm=="euclidean"){
      distances[,j]<-dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
    }
    if(distm=="correlation" | distm=="pearson"){
      distances[,j]<-t(-1*cor(cbind(centroids[,j],tdataMatrix),use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
    }
    if(distm=="spearman"){
      distances[,j]<-t(-1*cor(cbind(centroids[,j],tdataMatrix),method="spearman",use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
    }
  }
  
  scores<-apply(distances,1,min)
  prediction<-vector(length=tsamples)
  for(i in 1:tsamples){
    prediction[i]<-classLevels[match(scores[i],distances[i,])]
  }
  names(prediction)<-tsampleNames
  
  return(list(predictions=prediction,testData=tdataMatrix,distances=distances,centroids=centroids))
}

##---------------------------------------------------------------
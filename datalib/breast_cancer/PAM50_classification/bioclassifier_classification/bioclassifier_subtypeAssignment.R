###
# input variables for the subtype prediction script
###
library(ctc)
library(heatmap.plus)

paramDir<- "D:/Projects/Public_data_harmony/TCGA-8-11-2017_expression_subtyping/datalib/breast_cancer/PAM50_classification/bioclassifier_R" # the location of unchanging files such as the function library and main program
inputDir<- "D:/Projects/Public_data_harmony/TCGA-8-11-2017_expression_subtyping/subtype_results/BRCA"  # the location of the data matrix, and where output will be located

inputFile<- "TCGA-8-11-2017_BRCA_log2_TPM_quanNorm_for_classification.txt" # the input data matrix as a tab delimited text file
short<-"TCGA-8-11-2017_BRCA" # short name that will be used for output files

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

Package: TSED (tumor subtyping based on gene expression data)
Current version: v1.1
Author: Tenghui Chen (thchenb@gmail.com)

General utility:
	source("expression_data_tumor_subtyping.R")
	source("centroid_clustering.R")
	tumor_subtyping(indication=indication,inFile=inFile,prefix=NA,datalib_dir=datalib_dir)

Important message about the input and output files:
	1. At current version, we only support the cancer subtyping using the following types of gene expression data:
		1) BRCA: log2 (Microarray, TPM, FPKM, RPKM);
		2) COAD: log2 (Microarray, TPM, FPKM, RPKM);
		3) READ: log2 (Microarray, TPM, FPKM, RPKM);
		4) STAD: log2 (TPM);
		5) OV: log2 (TPM);
		6) LIHC: log2 (TPM);
		7) LUSC: log2 (TPM);
		8) LUAD: log2 (TPM);
		9) SKCM: log2 (TPM);
		10) PRAD: : log2 (TPM).		
	2. Regarding the input file format, please refer to the example input file in the example_data folder;
	3. Regarding the output file format, please refer to the example output file in the example_data folder.
	
Parameter setting:
	1. indication: the cancer indication name;
		At this version, 10 cancer indications were supported: 'LIHC','LUSC','LUAD','SKCM','PRAD','BRCA','STAD','OV','COAD','READ'. 
		The names were inherited from TCGA acronyms.
	2. inFile: the full path of input file;
	3. prefix: the prefix of the cancer subtype output file, default is NA;
	4. datalib_dir: the directory where the datalib folder locates in user's system.
	
	Note that the output file will reside in the same folder where the input file locates.

	
Brief message about the subtyping methods:
	1. For subtypes such as LIHC, LUSC, LUAD, SKCM and PRAD, the centroids for each subtype were provided. 
		The centroids were generated based on TCGA publications;
	2. For subtypes such as STAD and OV, the subtype was inferred based on GSVA results.
		The subtype specific genes were obtained from TCGA publications;
	3. For subtype such as BRCA, the subtype was inferred based on PAM50 methods. 
		In our package, the core algorithm of original PAM50 was retained, but the usage manner was modified. 
		You should expect the identical results from using our package and the original PAM50 package;
	4. For subtype such as COAD and READ, the subtype was inferred based on CMSclassifier. 
		In our package, we directly called the CMSclassifer for the classification.

Please contact Tenghui Chen (thchenb@gmail.com) if you may have any questions, any feedbacks would be appreciated.




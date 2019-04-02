## 
## Creator: K2
## Script Version: 1 
## Script status: Beta
## Description: Used for filtering results from Oncotator
## Support: Oncotator1.9.8.0+
## Steps: 1. Select variants with 'PASS' 2. Filter with Impact (High, Moderate) 3. Filter with LR or SVM ("D" or "-") 
## 4.  Filter with CADD score (exculde values which is "< 15") 5. Filter with Allele Frequency (Using General AF to exclude variant that greater than  10% (>0.1) )

# ----------------------------- Global Function ------------------------------------ #

## Please don't change it until you know what to do. 

rbind.match.columns <- function(input1=0, input2=0) {
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  if( is.null(n.input1) ){
    n.input1 = 0 
  }
  
  if( is.null(n.input2) ){
    n.input2 = 0 
  }
  
  if ( n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}

rbind.all.columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  if(! is.null(x.diff)){
    x[, c(as.character(y.diff))] <- NA        
  }
  
  if(! is.null(y.diff)){
    y[, c(as.character(x.diff))] <- NA
  }
  return(rbind(x, y))
}


# Step 0:  Configuration  

# Name this analysis
analysis_name <- "Glioblastoma_MoreFilter"

# Oncotator Version
OncotatorVersion <- "1.9.9.0" # 1.9.8.0,1.9.9.0

# select filtering steps
filterStatus <- TRUE
filterImapct <- TRUE
filterPredctionTools <- TRUE
filterCADDScores <- TRUE
filterAF <- TRUE



# All population databases 
col<- c("ExAC_ESP_AF_GLOBAL","ExAC_KG_AF_GLOBAL","1000gp3_EAS_AF","1000gp3_AF",
        "1000gp3_AFR_AF","1000gp3_AMR_AF","1000gp3_EUR_AF","dbNSFP_1000Gp1_AF",
        "dbNSFP_1000Gp1_AFR_AF", "dbNSFP_1000Gp1_AMR_AF","dbNSFP_1000Gp1_ASN_AF", 
        "dbNSFP_1000Gp1_EUR_AF", "dbNSFP_ESP6500_AA_AF","dbNSFP_ESP6500_EA_AF","ExAC_AF")

col<- c("1000gp3_AF","ExAC_AF")

# Thredshold scores
AF_freq <- 0.1
CADD_Score <- 15

# Input: Read whole directory  or file path 

fileNames <- Sys.glob("/home/note/MAFPlayground/Glioblastoma_all_merged_ready.maf.txt")

# Name the output ( if MergeResult TRUE )
file_result_name <-paste("/home/note/MAFPlayground/",analysis_name,".44WES.ready.filtered.maf",sep="")  

MergeResult <- TRUE ### <----- Not Finish , use merge by default
ToDiakoku <- FALSE

# -------------------------- local variable ----------------#
variantCount_Begining = 0
variantCount_Status = 0
variantCount_Imapct = 0
variantCount_PredctionTools =0
variantCount_CADDScores =0
variantCount_AF = 0
variantCount_Remainin = 0


# -------------------------- Let the hunt begin -----------------------------------
if (!exists("maf_file")){
  maf_file  <- vector()
}else{
  rm(maf_file)
}

for (fileName in fileNames) {
  
  print(paste0("Current working fileName: ", fileName))
  
    # if the merged dataset does exist, append to it
    temp_input_file <- read.csv(fileName, sep='\t', header=TRUE, comment.char = "#",na.strings=c("","NA","NaN", " ",".")  ,stringsAsFactors=FALSE ,check.names=FALSE)
    maf_begin <-temp_input_file
    if(filterStatus){
      ## 21 filters
      ##FILTER=<ID=artifact_in_normal,Description="artifact_in_normal">
      ##FILTER=<ID=bad_haplotype,Description="Variant near filtered variant on same haplotype.">
      ##FILTER=<ID=base_quality,Description="alt median base quality">
      ##FILTER=<ID=chimeric_original_alignment,Description="NuMT variant with too many ALT reads originally from autosome">
      ##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">
      ##FILTER=<ID=contamination,Description="contamination">
      ##FILTER=<ID=duplicate_evidence,Description="evidence for alt allele is overrepresented by apparent duplicates">
      ##FILTER=<ID=fragment_length,Description="abs(ref - alt) median fragment length">
      ##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">
      ##FILTER=<ID=low_avg_alt_quality,Description="Low average alt quality">
      ##FILTER=<ID=mapping_quality,Description="ref - alt median mapping quality">
      ##FILTER=<ID=multiallelic,Description="Site filtered because too many alt alleles pass tumor LOD">
      ##FILTER=<ID=n_ratio,Description="Ratio of N to alt exceeds specified ratio">
      ##FILTER=<ID=orientation_bias,Description="Orientation bias (in one of the specified artifact mode(s) or complement) seen in one or more samples.">
      ##FILTER=<ID=panel_of_normals,Description="Blacklisted site in panel of normals">
      ##FILTER=<ID=read_orientation_artifact,Description="orientation bias detected by the orientation bias mixture model">
      ##FILTER=<ID=read_position,Description="median distance of alt variants from end of reads">
      ##FILTER=<ID=str_contraction,Description="Site filtered due to contraction of short tandem repeat region">
      ##FILTER=<ID=strand_artifact,Description="Evidence for alt allele comes from one read direction only">
      ##FILTER=<ID=strict_strand_bias,Description="Evidence for alt allele is not represented in both directions">
      ##FILTER=<ID=t_lod,Description="Mutation does not meet likelihood threshold">
      
      
      if(OncotatorVersion == "1.9.9.0"){
        ## 21 filters
      maf_status  <- temp_input_file[which(
                                            temp_input_file$artifact_in_normal=='PASS'
                                           & temp_input_file$bad_haplotype=='PASS'
                                           & temp_input_file$base_quality=='PASS' 
                                           & temp_input_file$chimeric_original_alignment=='PASS' 
                                           & temp_input_file$clustered_events=='PASS' 
                                           & temp_input_file$contamination=='PASS'
                                           & temp_input_file$duplicate_evidence=='PASS'
                                           & temp_input_file$fragment_length=='PASS' 
                                           & temp_input_file$germline_risk=='PASS' 
                                           & temp_input_file$low_avg_alt_quality=='PASS' 
                                           & temp_input_file$mapping_quality=='PASS'
                                           & temp_input_file$multiallelic=='PASS' 
                                           & temp_input_file$orientation_bias=='PASS' 
                                           & temp_input_file$panel_of_normals=="PASS"
                                           & temp_input_file$read_orientation_artifact=='PASS'
                                           & temp_input_file$orientation_bias=='PASS' 
                                           & temp_input_file$read_position=='PASS'
                                           & temp_input_file$str_contraction=='PASS'
                                           & temp_input_file$strand_artifact=='PASS'
                                           & temp_input_file$t_lod=='PASS'
                                           & temp_input_file$strict_strand_bias=='PASS' 
                                          ),]
      }else if(OncotatorVersion == "1.9.8.0"){
        maf_status <- temp_input_file[ (temp_input_file$filter=='PASS'),]
      }else{
        
      }
      
      temp_input_file<-maf_status
      print(paste("Filter by Status completed " ))
    }
    
    if(filterImapct){
      # High:transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification
      # Moderate: inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant
 
        maf_imapct <- temp_input_file[which(temp_input_file$Ensembl_so_term=='transcript_ablation'
                                            | temp_input_file$Ensembl_so_term=='frameshift_variant' 
                                            |  temp_input_file$Ensembl_so_term=='stop_gained' 
                                            | temp_input_file$Ensembl_so_term=='splice_donor_variant'
                                            | temp_input_file$Ensembl_so_term=='splice_acceptor_variant'
                                            | temp_input_file$Ensembl_so_term=='stop_lost' 
                                            |  temp_input_file$Ensembl_so_term=='start_lost' 
                                            | temp_input_file$Ensembl_so_term=='transcript_amplification'
                                            | temp_input_file$Ensembl_so_term=='inframe_insertion' 
                                            |  temp_input_file$Ensembl_so_term=="missense" 
                                            |  temp_input_file$Ensembl_so_term=='inframe_deletion'
                                            |  temp_input_file$Ensembl_so_term=='protein_altering_variant'),]   
        

      
      temp_input_file<-maf_imapct
      print(paste("Filter by Impact completed " ))
    }
    
    ## ++++ Filter by PredctionTools +++++
    if(filterPredctionTools){
      # dbNSFP annotation, Our logistic regression (LR) and support vector machine (SVM) based ensemble prediction score, which incorporated 10 scores (SIFT, PolyPhen-2 HDIV, PolyPhen-2 HVAR, GERP++, MutationTaster, Mutation Assessor, FATHMM, LRT, SiPhy, PhyloP) and the maximum frequency observed in the 1000 genomes populations.
      # Larger value means the SNV is more likely to be damaging.
      
      maf_SVM <- temp_input_file[which(grepl("D", temp_input_file$dbNSFP_RadialSVM_pred) | is.na(temp_input_file$dbNSFP_RadialSVM_pred)),]
      
      maf_LR <- maf_SVM[which(grepl("D", maf_SVM$dbNSFP_LR_pred) | is.na(maf_SVM$dbNSFP_LR_pred)),]

      temp_input_file<-maf_LR
      print(paste("Filter by Predction Tools completed " ))
    }
    
    ## ++++ Filter by CADDScores +++++
    if(filterCADDScores){
      
      # Preprocess CAD , because some rows contain two score in one cell (e.g 3.59|3.59), so we must replace it with only one value (e.g. 3.59).
      temp_input_file$dbNSFP_CADD_phred <- gsub("\\|.*", "", temp_input_file$dbNSFP_CADD_phred)
     # FixCADD<-temp_input_file
      
      # Covert CADD_phred to numeric
      temp_input_file$dbNSFP_CADD_phred <- as.numeric(as.character(temp_input_file$dbNSFP_CADD_phred))
     # FixCADD_1 <-temp_input_file
      
      maf_CADD <- temp_input_file[ which(temp_input_file$dbNSFP_CADD_phred >=15 | is.na(temp_input_file$dbNSFP_CADD_phred) ) ,]
      
      temp_input_file<-maf_CADD
      print(paste("Filter by CADD score completed " ))
    }
    
    ## ++++ Filter by Freqency +++++
    if(filterAF){
      
      
      afdata <- temp_input_file
      afdata <- afdata[,col]
      afdata[col] <- sapply(afdata[col],as.numeric)

      afdata <- afdata <= AF_freq
      afdata[is.na(afdata)] <- TRUE # keep NA or NULL value 
      
      temp_input_file <- temp_input_file[which(apply(as.data.frame(afdata), 1,all)),]    
      
      print(paste("Filter by AF completed " ))
    }
  
      maf_file <- rbind.all.columns(maf_file, temp_input_file)    
      rm(temp_input_file)
      print("********* Binding Complete ************* ")
      

      
}



if(ToDiakoku){
 # maf_status[( !duplicated(maf_status$Genome_Change) &  !duplicated(maf_status$Annotation_Transcript) &  !duplicated(maf_status$Codon_Change) &  !duplicated(maf_status$cDNA_Change)  &  !duplicated(maf_status$Transcript_Strand)  )	, ]
 
   maf_file<-maf_file[( !duplicated(maf_file$Genome_Change) )	, ]
}

write.table(maf_file, file = file_result_name, sep = "\t",row.names = FALSE, col.names = TRUE , quote = FALSE, append = FALSE)

print(paste0("********* Done ************* " ,file_result_name))
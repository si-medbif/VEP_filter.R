# VEP_filter.R
Filtration SUDS variant from VEP by R

###How to run
## !! limitation this script using only the annotation from VEP
#1) Import a dataset from VEP and name it "vepdata".
#2) Import a dataset from snpeff and name it "snpeffdata".
#3) If candidate gene dataset from VEP skip step 3 and step 4, then step 5 create dataname to data3
#4) set working directory 
setwd("D:/")
#4)Press the "Source" button on the top right corner of this panel.


##Method for filtration variant from VEP
#Step1 Recode impact variant to numeric and then sort by descending
vepdata$sort_impact <- as.integer(ifelse(vepdata$IMPACT=="LOW","1",
                                         ifelse(vepdata$IMPACT=="MODIFIER","2",
                                                ifelse(vepdata$IMPACT=="MODERATE","3",
                                                       ifelse(vepdata$IMPACT=="HIGH","4","")))))
vepdata <- vepdata[order(-vepdata$sort_impact),] 


#Step 2 Remove duplicate by location
data <- vepdata[!duplicated(vepdata$Location), ]


#Step 3 Split position by location (if you want to add genotype from snpeff) 
#library(tidyr)
#data2 <- separate(data, col = "Location",into = c("Chromosome", "POS","POS2"),
#                  sep = "[:-]",fill = "right")


# Step 4 Merge data from VEP and SnpEff 
#data3 <- merge(data2, snpeffdata, by = "POS", all= TRUE)


#Step 5 Filter by MetaSVM (D = Deleterious, T = Tolerable)
## If MetaSVM reported T, exclude from data   
### !!! If candidate gene create the data name "data3" 
data3 <- data
svmdata <- data3[which(data3$MetaSVM_pred == "D" |data3$MetaSVM_pred == "-"),]


#Step 6  Filter Clinical significance("CLIN_SIG") by Clin Var
## If CLIN_SIG reported benign, likely benign, exclude from dataset 
uncertain <- grepl("uncer", svmdata$CLIN_SIG) #find "uncertain_significance"
notpro <- grepl("pro", svmdata$CLIN_SIG)  #find "not_provided" 
benign <- grepl("ben", svmdata$CLIN_SIG)  #find "benign or likely_benign"
patho <- grepl("patho", svmdata$CLIN_SIG) #find "pathogenic or likely_pathogenic"
na <- grepl("-", svmdata$CLIN_SIG) #find "-"
select_clinsig = patho|((uncertain|na)&!notpro&!benign) #logical for filter out benign and likely benign 
clinsig<- cbind.data.frame(svmdata, select_clinsig)
clinsig <-clinsig[which(clinsig$select_clinsig == "TRUE"),]


#step 7 select column AF 
afdata <- clinsig
afdata <- afdata[,grep("AF",names(afdata))]
drop <- c("1000Gp1_AFR_AC", "ExAC_AFR_AC")
afdata <- afdata[,!(names(afdata) %in% drop)]
col<- c("AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF", 
        "gnomAD_AMR_AF",  "gnomAD_ASJ_AF",  "gnomAD_EAS_AF",  "gnomAD_FIN_AF",  "gnomAD_NFE_AF","gnomAD_OTH_AF",
        "gnomAD_SAS_AF","1000Gp1_AF","1000Gp1_AFR_AF", "1000Gp1_AMR_AF","1000Gp1_ASN_AF", "1000Gp1_EUR_AF",
        "ARIC5606_AA_AF", "ARIC5606_EA_AF", "ESP6500_AA_AF","ESP6500_EA_AF","ExAC_AF","ExAC_AFR_AF",
        "ExAC_AMR_AF","ExAC_Adj_AF","ExAC_EAS_AF","ExAC_FIN_AF","ExAC_NFE_AF","ExAC_SAS_AF")
afdata[col] <- sapply(afdata[col],as.numeric)
sapply(afdata, class)


# Step 8 Filter out MAF > 10% (depend on criteria of your study) in MAF in many studies
afdata <- afdata <= 0.1
afdata[is.na(afdata)] <- TRUE
mafdata <- clinsig[which(apply(as.data.frame(afdata), 1,all)),]


# step 9 filter out modifier and low impact for vep (if you want to filter)
filterimpact <- mafdata[which(mafdata$IMPACT == "HIGH"| mafdata$IMPACT == "MODERATE"),]


#Step 10 save output to working dicterory 
write.csv(filterimpact,"test_candidate_SUDS_variant.csv",na = "")



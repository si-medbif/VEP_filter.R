##
## Creator: K2
## Script Version: 1
## Script status: Production (inherit form vep_filter.chunky.template)
## Description: Used for filtering breast cancer germline which produced from VEP
## Support: VEP95.2, dbSNP 3.5a, Whole Genome Sequnceing, Hg19, GATK-V.4.1.0
##

library(doParallel)

# ------------ Global Funtion ------------------ #
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

# -------------------------------------------- #

# ------------ Global Variable ------------------ #

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

filename <- "b212wgs_annotation_output.vep.chr17.hg19.gz"
index = 0
counter = 0
chunks = 100000
con = file(filename, "r")
skip_headerLine = 343
data_header <- names(read.csv(filename, nrows=1, sep="\t", skip = skip_headerLine ,header = TRUE , check.names = FALSE))
ready_to_write <- data.frame(matrix(ncol =  length(data_header), nrow = 0))
names(ready_to_write) <- data_header

# ----------------------------------------------- #

# ------------ Script Start ------------------ #
system.time(
repeat{

  tryCatch({
    vepdata <-read.table(con, nrows=chunks, header=FALSE, fill = FALSE, sep="\t",comment.char = "#" , check.names = FALSE ,quote = "",as.is = TRUE )
    # do processing on dataChunk (i.e adding header, converting data type)

    colnames(vepdata) <- data_header
    finalMatrix <- foreach(i = nrow(vepdata), .combine=cbind) %dopar% {

      #Step1 Recode impact variant to numeric and then sort by descending
      # vepdata$sort_impact <- as.integer(ifelse(vepdata$IMPACT=="LOW","1",
      #                                          ifelse(vepdata$IMPACT=="MODIFIER","2",
      #                                                 ifelse(vepdata$IMPACT=="MODERATE","3",
      #                                                        ifelse(vepdata$IMPACT=="HIGH","4","")))))
      # vepdata <- vepdata[order(-vepdata$sort_impact),]

      filterimpact <- vepdata[which(vepdata$IMPACT %in% c("HIGH","MODERATE")),]


      #Step 2  Filter Clinical significance("CLIN_SIG")
      ## If CLIN_SIG reported benign and likely benign,then exclude them from dataset
#      uncertain <- grepl("uncer", filterimpact$CLIN_SIG) #find "uncertain_significance"
#      notpro <- grepl("pro", filterimpact$CLIN_SIG)  #find "not_provided"
#      benign <- grepl("ben", filterimpact$CLIN_SIG)  #find "benign or likely_benign"
#      patho <- grepl("patho", filterimpact$CLIN_SIG) #find "pathogenic or likely_pathogenic"
#      na <- grepl("-", filterimpact$CLIN_SIG) #find "-"
#      select_clinsig = patho|((uncertain|na)&!notpro&!benign) #logical statement to filter out benign and likely benign
       clinsig <- filterimpact[!grepl("2|3",filterimpact$clinvar_clnsig),]
       ## known pathogenic variants from clinvar
       clinvar_pathogenic <- vepdata[grepl("4|5",vepdata$clinvar_clnsig) & !grepl("2|3",vepdata$clinvar_clnsig),]

      #Step 3 Filter by MetaSVM (D = Deleterious, T = Tolerable) on MODERATE Impact records
      ## If MetaSVM reported T, exclude it from data

      # Covert CADD_phred to numeric
      clinsig$CADD_phred <- as.numeric(as.character(clinsig$CADD_phred))
      # Keep all (MetaSVM == "D" or "-") AND (CADD Phred Scale >= 15)
      temp_data <- clinsig[ !(clinsig$IMPACT == "MODERATE" & (clinsig$CADD_phred < 15 | clinsig$MetaSVM_pred == "T")  )  , ]
      #temp_data <- clinsig[ !(clinsig$IMPACT == "MODERATE" & (clinsig$CADD_phred < 15)  )  , ]
      #temp_data <- temp_data[ !(temp_data$IMPACT == "MODERATE" & (temp_data$MetaSVM_pred == "T")  )  , ]

      #step 4 Filter out common variants
      #select column AF
      afdata <- temp_data[ ,grep("^AF$|*_AF",names(temp_data))]
      col<- c("AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF",
              "gnomAD_AMR_AF",  "gnomAD_ASJ_AF",  "gnomAD_EAS_AF",  "gnomAD_FIN_AF",  "gnomAD_NFE_AF","gnomAD_OTH_AF",
              "gnomAD_SAS_AF","1000Gp3_AF","1000Gp3_AFR_AF", "1000Gp3_AMR_AF", "1000Gp3_EUR_AF","1000Gp3_SAS_AF",
              "ESP6500_AA_AF","ESP6500_EA_AF","ExAC_AF","ExAC_AFR_AF",
              "ExAC_AMR_AF","ExAC_Adj_AF","ExAC_EAS_AF","ExAC_FIN_AF","ExAC_NFE_AF","ExAC_SAS_AF")
      afdata[col] <- sapply(afdata[col],as.numeric)

      #Filter out a common variants with allele frequency in any population  that is greater than 0.1
      afdata <- afdata <= 0.1
      afdata[is.na(afdata)] <- TRUE
      mafdata <- temp_data[which(apply(as.data.frame(afdata), 1,all)),]


      temp_data <- rbind(clinvar_pathogenic,mafdata) # combined known_pathogenic + prioritized variants.

    }

    ready_to_write <- rbind.match.columns(ready_to_write, finalMatrix)

    #check if file end has been reached and break from repeat
    if(nrow(vepdata) < chunks){
      break
    }

    #increment the index to read next chunk
    index <- index + 1
    print(paste('Processing rows:', index * chunks))


  }, error=function(err) {
    ## matching condition message only works when message is not translated
    if (identical(conditionMessage(err), "no lines available in input"))
      break
    else stop(err)
  })
}
)
## Cleans all rows with NAs in every columns (Unexpected event)
ready_to_write<-ready_to_write[rowSums(is.na(ready_to_write)) != ncol(ready_to_write), ]

write.csv(ready_to_write,"VEPfilter.1.chr17.csv",na = "",row.names = F)
close(con)

#stop cluster
stopCluster(cl)


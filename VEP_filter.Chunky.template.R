## 
## Creator: K2
## Script Version: 1
## Script status: Template
## Read file -  chunk by chunk
## Support VEP95.2 , dbSNP 3.5a
## 


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

filename <- "b212wgs_annotation_output.vep.chr17.hg19.gz"
index <- 0
counter <- 0
chunks <- 50000
con = file(filename, "r")
skip_headerLine = 343

data_header <- names(read.csv(filename, nrows=1, sep="\t", skip = skip_headerLine ,header = TRUE , check.names = FALSE))
ready_to_write <- data.frame(matrix(ncol =  length(data_header), nrow = 0))
names(ready_to_write) <- data_header

repeat{
  
  tryCatch({
    vepdata <-read.table(con, nrows=chunks, header=FALSE, fill = FALSE, sep="\t",comment.char = "#" , check.names = FALSE)
    # do processing on dataChunk (i.e adding header, converting data type) 
    
    colnames(vepdata) <- data_header
    
    #Step1 Recode impact variant to numeric and then sort by descending
    vepdata$sort_impact <- as.integer(ifelse(vepdata$IMPACT=="LOW","1",
                                             ifelse(vepdata$IMPACT=="MODIFIER","2",
                                                    ifelse(vepdata$IMPACT=="MODERATE","3",
                                                           ifelse(vepdata$IMPACT=="HIGH","4","")))))
    vepdata <- vepdata[order(-vepdata$sort_impact),] 
    
    #Step 2 Remove duplicate by location
    #data <- vepdata[!duplicated(vepdata$Location), ]
    
    
    #Step 5 Filter by MetaSVM (D = Deleterious, T = Tolerable)
    ## If MetaSVM reported T, exclude from data   
    ### !!! If candidate gene create the data name "data3" 
    data3 <- vepdata
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
    
    
    # step 9 filter out modifier and low impact for vep (if you want to filter)
    filterimpact <- clinsig[which(clinsig$IMPACT == "HIGH"| clinsig$IMPACT == "MODERATE"),]
    
    ready_to_write <- rbind.match.columns(ready_to_write, filterimpact)
    
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
write.csv(ready_to_write,"filter.csv",na = "",row.names = F)
close(con)
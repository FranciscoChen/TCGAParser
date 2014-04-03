#!/usr/local/env Rscript
#
# Automatization of TCGA file downloading and barcode pairing.
# Data frame filtering and storage in PostgreSQL.
#
# Francisco Chen, GPL v2
# 18th Feb 2014

#' Get matching links in a URL
#' 
#' Given a URL and a regular expression, returns all links in the URL matching the regular expression along its last modified date.
#' 
#' @usage getLinks(URL, REGEXPR)
#' @param URL Link to a website to scrape.
#' @param REGEXPR Regular expression to match to the scraped contents.
#' @details Requires RCurl package to getURL from https. 
#' @return Returns a list[[i]][[1|2]] with i matches and two sublevels for every match: 1 is the link and 2 is the last modified date (if any). 
#' If the given URL does not exist or there are no matches, returns nothing.
getLinks <- function (URL, REGEXPR){
  if(exists("getURL") == FALSE){
    require (RCurl)
  }
  if (url.exists(URL)){
    scrapedData <- getURL(URL)
    links <- strsplit(scrapedData, "href=\"")[[1]]
    linkIndices <- grep(REGEXPR, links)
    if (length(linkIndices > 0)){
      list <- vector ("list", length(linkIndices))
      for (i in 1:length(linkIndices)){
        list[[i]][[1]] <- regmatches(links[linkIndices[i]], regexpr(REGEXPR, links[linkIndices[i]]))
        list[[i]][[2]] <- regmatches(links[linkIndices[i]], regexpr("[[:alnum:]]{4}-[[:alnum:]]{2}-[[:alnum:]]{2}[[:blank:]][[:alnum:]]{2}[[:punct:]][[:alnum:]]{2}", links[linkIndices[i]]))
      }
      return (list)
    }
  }
}

#' Get a data frame from a URL
#' 
#' Download a file as a temporary .txt file in the working directory and read it into R as a data frame.
#' 
#' @usage getDataframe(URL, header = TRUE, sep = "\t", quote = "\"'", dec = ".", 
#' row.names, col.names, as.is = !stringsAsFactors, na.strings = "NA", 
#' colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, 
#' fill = !blank.lines.skip, strip.white = FALSE, blank.lines.skip = TRUE, 
#' comment.char = "#", allowEscapes = FALSE, flush = FALSE, 
#' stringsAsFactors = default.stringsAsFactors(), fileEncoding = "", 
#' encoding = "unknown", text)
#' @param URL Link to a text file.
#' @details The file is read with tab as separators and expecting a header. All other arguments are as in read.table.
#' @return Returns a data frame or nothing.
#' @seealso read.table
getDataframe <- function (URL, header = TRUE, sep = "\t", quote = "\"'", dec = ".", 
                          row.names, col.names, as.is = !stringsAsFactors, na.strings = "NA", 
                          colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, 
                          fill = !blank.lines.skip, strip.white = FALSE, blank.lines.skip = TRUE, 
                          comment.char = "#", allowEscapes = FALSE, flush = FALSE, 
                          stringsAsFactors = default.stringsAsFactors(), fileEncoding = "", 
                          encoding = "unknown", text) {
  if(exists("getBinaryURL") == FALSE){
    require (RCurl)
  }
  if (url.exists(URL)){
    # Downloads the data as .txt
    binary <- getBinaryURL(URL, ssl.verifypeer=FALSE)
    connection <- file(paste(writePath, "/", "tempData.txt", sep = ""), open = "wb")
    writeBin(binary, connection)
    close(connection)
    # Reads the .txt file.
    if (file.exists(paste(writePath, "/", "tempData.txt", sep = "")) == TRUE){
      data <- read.table(paste(writePath, "/", "tempData.txt", sep = ""), header, sep, quote, dec, 
                         row.names, col.names, as.is, na.strings, 
                         colClasses, nrows, skip, check.names, 
                         fill, strip.white, blank.lines.skip, 
                         comment.char, allowEscapes, flush, 
                         stringsAsFactors, fileEncoding, 
                         encoding, text)
      file.remove(paste(writePath, "/", "tempData.txt", sep = ""))
    }
    if (exists("data") == TRUE){
      return (data)
    }
  }
}

#' Fix a data frame
#' 
#' Format a data frame to a common standard. Varies from array to array.
#' 
#' @usage fixDataframe(dataframe, array)
#' @param dataframe Data frame to fix.
#' @param array Array from which the data was processed. 
#' Currently supported arrays: humanmethylation450, illuminahiseq_rnaseqv2.
#' @return Returns the fixed data frame or nothing.
fixDataframe <- function(dataframe, array) {
  switch(array,
         humanmethylation450 = {
           # Looks for "Beta_value" in the file head and returns its position.
           matrixPosition <- grep("Beta_value", as.matrix(head(dataframe)), ignore.case = TRUE)
           # If "Beta_value" was found, it will calculate its row (grepRow) and column (grepCol).
           if (length(matrixPosition) > 0) {
             if (grepl("Beta_value", as.matrix(head(dataframe))[matrixPosition], ignore.case = TRUE)){ 
               
               matrixRows <- nrow(as.matrix(head(dataframe))) # How many total rows are in the header?
               grepCol <- 1
               
               while (matrixPosition > matrixRows){
                 matrixPosition <- matrixPosition - matrixRows
                 grepCol <- grepCol +1
               }
               
               grepRow <- matrixPosition
               
               names(dataframe) <- as.matrix(dataframe[grepRow,]) # Row names are substituded for the row that contains "Beta_value".
               data <- dataframe[-grepRow,] # Removes the row that contained "Beta_value".
               data <- data[1:485512,] # Removes all the rs probes.
               return (data)
             }
           }
           stop ("'Beta_value' not found in the dataframe head: ",as.matrix(head(dataframe)))
         },
         illuminahiseq_rnaseqv2 = {
           return (dataframe)
         })
}

#' Initialise variables
#' 
#' Declares and assigns variables like URLs and regular expressions.
#' 
#' @usage init(cancer, array)
#' @param cancer The cancer type we intend to work with, abbreviated as in the TCGA.
#' @param array Array from which the data was processed.  
#' @seealso https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm
#' Currently fully supported arrays: humanmethylation450, illuminahiseq_rnaseqv2.
init <- function (cancer, array){
  accessRoot <<- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
  cancerName <<- cancer
  arrayName <<- array
  switch(array,
         humanmethylation450 = {
           arrayPath <<- "/cgcc/jhu-usc.edu/humanmethylation450/methylation/"
           archiveRegexpr <<- "jhu-usc.edu_[[:alpha:]]*.HumanMethylation450.Level_3.[0-9]+.6.[0-9]/"
           fileRegexpr <<-  "jhu-usc.edu_[[:alpha:]]*.HumanMethylation450.[[:alnum:]]*.lvl-3.TCGA-[[:alnum:]]*-[[:alnum:]]*-0[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*.txt"
           barcodeRegexpr <<- "TCGA-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*"
           value <<- "Beta_value"
           magetabArchiveRegexpr <<- "jhu-usc.edu_[[:alpha:]]{2,4}.HumanMethylation450.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "jhu-usc.edu_[[:alpha:]]{2,4}.HumanMethylation450.[[:digit:]]+.[[:digit:]]+.0.sdrf.txt"
         },
         agilentg4502a_07_3 = {
           arrayPath <<- "/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/"
           magetabArchiveRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.AgilentG4502A_07_3.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.AgilentG4502A_07_3.sdrf.txt"
           },
         illuminaga_rnaseq = {
           arrayPath <<- "/cgcc/unc.edu/illuminaga_rnaseq/rnaseq/"
           magetabArchiveRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaGA_RNASeq.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaGA_RNASeq.[[:digit:]]+.sdrf.txt"
         },
         illuminaga_rnaseqv2 = {
           arrayPath <<- "/cgcc/unc.edu/illuminaga_rnaseqv2/rnaseqv2/"
           magetabArchiveRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaGA_RNASeqV2.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaGA_RNASeqV2.[[:digit:]]+.[[:digit:]]+.0.sdrf.txt"
         },
         illuminahiseq_rnaseq = {
           arrayPath <<- "/cgcc/unc.edu/illuminahiseq_rnaseq/rnaseq/"
           magetabArchiveRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaHiSeq_RNASeq.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaHiSeq_RNASeq.[[:digit:]]+.sdrf.txt"
         },
         illuminahiseq_rnaseqv2 = {
           arrayPath <<- "/cgcc/unc.edu/illuminahiseq_rnaseqv2/rnaseqv2/"
           archiveRegexpr <<- "unc.edu_[[:alpha:]]*.IlluminaHiSeq_RNASeqV2.Level_3.[0-9]+.7.[0-9]/"
           fileRegexpr <<-  "unc.edu.[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*-[[:alnum:]]*.[[:alnum:]]*.bt.exon_quantification.txt"
           barcodeRegexpr <<- "[[:alnum:]]{8}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}"
           value <<- "RPKM"
           magetabArchiveRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaHiSeq_RNASeqV2.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaHiSeq_RNASeqV2.[[:digit:]]+.[[:digit:]]+.0.sdrf.txt"
         },
         illuminahiseq_totalrnaseqv2 = {
           arrayPath <<- "/cgcc/unc.edu/illuminahiseq_totalrnaseqv2/totalrnaseqv2/"
           magetabArchiveRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaHiSeq_TotalRNASeqV2.mage-tab.[[:digit:]]+.[[:digit:]]+.0/"
           magetabFileRegexpr <<- "unc.edu_[[:alpha:]]{2,4}.IlluminaHiSeq_TotalRNASeqV2.[[:digit:]]+.[[:digit:]]+.0.sdrf.txt"
         }
         )
}

#' Download data from the TCGA
#' 
#' Downloads all the data related to a cancer type and an array as specified in init.
#' 
#' @include init
#' @include getLinks
#' @include getDataframe
#' @include fixDataframe
#' @usage download(cancer, array)
#' @details Takes considerable time and memory to complete. It will generate a .txt file for each TCGA folder.
download <- function(cancer, array){
  init(cancer, array)
  archiveList <- getLinks(paste(accessRoot, cancerName, arrayPath, sep = ""), archiveRegexpr)
  
  # First time: initialize the dataframe.
  if (file.exists(paste(cancerName, "_dataframe_", arrayName, ".txt", sep = "")) == FALSE){
    fileList <- getLinks(paste(accessRoot, cancerName, arrayPath, archiveList[[1]][[1]], sep = ""), fileRegexpr)
    reference <- getDataframe(paste(accessRoot, cancerName, arrayPath, archiveList[[1]][[1]], fileList[[1]][[1]], sep = ""))
    reference <- fixDataframe(reference, array)
    reference[colnames(reference) == value] <- NULL
    switch(array,
           humanmethylation450 = {
             reference[colnames(reference) == "Gene_Symbol"] <- NULL
             annotations <- getDataframe("www.ncbi.nlm.nih.gov//geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL16304&id=47833&db=GeoDb_blob89")
             reference <- cbind(reference, annotations$Closest_TSS)
             reference <- cbind(reference, annotations$Distance_closest_TSS)
             reference <- cbind(reference, annotations$Closest_TSS_gene_name)
           },
           illuminahiseq_rnaseqv2 = {
             reference[colnames(reference) == "raw_counts"] <- NULL
             reference[colnames(reference) == "median_length_normalized"] <- NULL
             exon <- "^chr([0123456789XY]+):([[:digit:]]*)-([[:digit:]]*):([[:punct:]])$"
             reference$chromosome <- apply(reference[,1, drop = FALSE], 2, function(x) as.character(gsub(exon, "\\1", x)))
             reference$start <- apply(reference[,1, drop = FALSE], 2, function(x) as.character(gsub(exon, "\\2", x)))
             reference$end <- apply(reference[,1, drop = FALSE], 2, function(x) as.character(gsub(exon, "\\3", x)))
             reference$strand <- apply(reference[,1, drop = FALSE], 2, function(x) as.character(gsub(exon, "\\4", x)))
           }
    )
    write.table(reference, file = paste(writePath, "/", cancerName, "_dataframe_", arrayName, ".txt", sep = ""), sep = "\t", quote = FALSE)
    rm(reference)
    gc()
  }
  
  for (i in 1:length(archiveList)){
    #Get all files in every folder
    fileList <- getLinks(paste(accessRoot, cancerName, arrayPath, archiveList[[i]][[1]], sep = ""), fileRegexpr)
    if (length(fileList) > 0){
      df <- read.table(paste(readPath, "/", cancerName, "_dataframe_", arrayName, ".txt", sep = ""), header = TRUE, sep = "\t")
      for (j in 1: length(fileList)){
        #Append every file to our dataframe
        URL <- paste(accessRoot, cancerName, arrayPath, archiveList[[i]][[1]], fileList[[j]][[1]], sep = "")
        newdf <- getDataframe(URL)
        newdf <- fixDataframe(newdf, array)
        colnames(newdf)[colnames(newdf) == value] <- regmatches(URL, regexpr(barcodeRegexpr, URL))
        newdf <- newdf[,regmatches(URL, regexpr(barcodeRegexpr, URL)), drop = FALSE]
        df <- cbind(df, newdf)
      }
      #Finally, save the table
      archiveName <- regmatches(archiveList[[i]][[1]], regexpr(substr(archiveRegexpr, 1, nchar(archiveRegexpr)-1),archiveList[[i]][[1]]))
      write.table(df, file = paste(writePath, "/", cancerName, "_dataframe_", archiveName, ".txt", sep = ""), append = FALSE, sep = "\t", quote = FALSE)
      rm(df)
      gc()
    }
  }
}

#' Get TCGA barcodes
#' 
#' Gets a dataframe containing all the barcodes assigned to every sample in a specific cancer type and array.
#' 
#' @include init
#' @include getLinks
#' @include getDataframe
#' @usage getBarcodes(cancer, array)
#' @return A data frame of the mage-tab file or NA.
getBarcodes <- function (cancer, array){
  init(cancer, array)
    archives <- getLinks(paste(accessRoot, cancerName, arrayPath, sep = ""), magetabArchiveRegexpr)
  if (length(archives) > 0){
    if (length(archives) > 1){
      # Get the most recent one.
      dates <- NULL
      for (j in 1:length(archives)){
        dates[[j]] <- archives[[j]][[2]]
      }
      dates <- sort(dates, decreasing = TRUE)
      archives <- archives[grep(dates[1], archives)]
    }
    files <- getLinks(paste(accessRoot, cancerName, arrayPath, archives[[1]][[1]], sep = ""), magetabFileRegexpr)
    if (length(files) > 0){
      if (length(files) > 1){
        # Get the most recent one.
        dates <- NULL
        for (j in 1:length(files)){
          dates[[j]] <- files[[j]][[2]]
        }
        dates <- sort(dates, decreasing = TRUE)
        files <- files[grep(dates[1], files)] 
      }
      return (getDataframe(paste(accessRoot, cancerName, arrayPath, archives[[1]][[1]], files[[1]][[1]], sep = "")))
    }
  }
  return (NA)
}

#' Pair TCGA barcodes
#' 
#' Pairs barcodes between two arrays for a specific cancer type.
#' 
#' @usage pairBarcodes (cancer, array1, array2)
#' @details Filters for the first 16 TCGA barcode characters, which refer to participant, sample type and vial. Duplicates are removed, along with normal and control tissue samples.
#' Finally, it will write into disk the data frames for each array.
pairBarcodes <- function (cancer, array1, array2){
  a <- getBarcodes(cancer, array1)
  a$Barcodes <- strtrim(a$Comment..TCGA.Barcode., 16)
  a <- subset(a, duplicated(Barcodes) == FALSE)
  
  b <- getBarcodes(cancer, array2)
  b$Barcodes <- strtrim(b$Comment..TCGA.Barcode., 16)
  b <- subset(b, duplicated(Barcodes) == FALSE)
  b <- subset(b, Barcodes %in% a$Barcodes)
  b <- subset(b, grepl("TCGA-[[:alnum:]]*-[[:alnum:]]*-0[[:alnum:]]*", Barcodes))
  
  a <- subset(a, Barcodes %in% b$Barcodes)
  a <- subset(a, grepl("TCGA-[[:alnum:]]*-[[:alnum:]]*-0[[:alnum:]]*", Barcodes))
  
  write.table(a, file = paste(writePath, "/", cancer, "_commonPatients_", array1, ".txt", sep = ""), append = FALSE, sep = "\t", quote = FALSE)
  write.table(b, file = paste(writePath, "/", cancer,"_commonPatients_", array2, ".txt", sep = ""), append = FALSE, sep = "\t", quote = FALSE)
}

#' Filter TCGA data by paired barcodes
#'
#' Filters cancer samples by paired barcodes and prepares the table to be copied to PostgreSQL.
#' 
#' @include init
#' @include getLinks
#' @usage filterBarcodes(cancer, array)
#' @details Uses the files generated by download and pairBarcodes. 
#' It will generate and save three data frames: one containing information for every probe, one containing the sample values for each probe, and one containing the sample barcodes (reference).
filterBarcodes <- function (cancer, array){
  init(cancer, array)
  archiveList <- getLinks(paste(accessRoot, cancer, arrayPath, sep = ""), archiveRegexpr)
  all <- read.table(paste(readPath, "/", cancerName, "_dataframe_", arrayName, ".txt", sep = ""), header = TRUE, sep = "\t")
  all <- all[,1, drop = FALSE]
  
  pairedBarcodes <- read.table(paste(readPath, "/", cancerName, "_commonPatients_", arrayName, ".txt", sep = ""), sep = "\t", header = TRUE)
  for (i in 1:length(archiveList)){
    archiveName <- regmatches(archiveList[[i]][[1]], regexpr(archiveRegexpr, archiveList[[i]][[1]]))
    archiveName <- substr(archiveName, 1, nchar(archiveName)-1)
    if (file.exists(paste(cancerName, "_dataframe_", archiveName, ".txt", sep = "")) == TRUE){
      a <- read.table(paste(readPath, "/", cancerName, "_dataframe_", archiveName, ".txt", sep = ""), header = TRUE, sep = "\t")
      colnames(a) <- gsub("[.]","-", colnames(a))
      switch(array,
             humanmethylation450 = {
               a <- subset (a, select = colnames(a) %in% pairedBarcodes$Comment..TCGA.Barcode.)
             },
             illuminahiseq_rnaseqv2 = {
               colnames(a)[grep(barcodeRegexpr, colnames(a))] <- regmatches(colnames(a), regexpr(barcodeRegexpr, colnames(a)))
               a <- subset (a, select = colnames(a) %in% pairedBarcodes$Extract.Name)
               colnames(a) <- as.character(lapply(colnames(a), function(x) x <- as.character(pairedBarcodes[grep(x, pairedBarcodes$Extract.Name),]$Comment..TCGA.Barcode.)))
             }
      )
      all <- cbind (all, a)
      rm(a)
    }
  }
  probeinfo <- read.table(paste(readPath, "/", cancerName, "_dataframe_", arrayName, ".txt", sep = ""), header = TRUE, sep = "\t")
  write.table(probeinfo, paste(writePath, "/", cancerName, "_probeinfo_", arrayName, "_to_SQL.txt", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  rm(probeinfo)
  all <- all[,-1]
  all <- all[,order(colnames(all))]
  samples <- all[0,]
  colnames(samples) <- gsub("-","_", colnames(samples))
  write.table(samples, paste(writePath, "/", cancerName, "_sampleinfo_", arrayName, "_to_SQL.txt", sep = ""), quote = FALSE, sep = "\t")
  id <- read.table(paste(readPath, "/", cancerName, "_dataframe_", arrayName, ".txt", sep = ""), header = TRUE, sep = "\t")
  id <- id[,1, drop = FALSE]
  all <- cbind (id, all)
  write.table(all, paste(writePath, "/", cancerName, "_", arrayName, "_to_SQL.txt", sep = ""), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  rm(all)
  gc()
}

#' Analyze the NA distribution
#' 
#' Analyze the NA distribution in the filtered TCGA samples, per probe and sample (uses the files from filterBarcodes).
#' Also, get the standard deviation for every probe.
#' 
#' @include init
#' @usage analyzeNAs(cancer, array)
#' @details It will write two text files to disk, one containing the number of NAs and standard deviation of every probe and the other containing the number of NAs for every sample.
analyzeNAs <- function (cancer, array){
  init(cancer, array)
  df <- read.table(paste(readPath, "/", cancerName, "_", arrayName, "_to_SQL.txt", sep = ""), sep = "\t", header = TRUE)
   
  probeTable <- as.data.frame (matrix(nrow = nrow(df), ncol = 3))
  sampleTable <- as.data.frame (matrix(nrow = 2, ncol = ncol(df[,2:ncol(df)])))
    
  probeTable[,1] <- df[,1]
  probeTable[,2] <- apply(df[,2:ncol(df)], 1, function(z) sum(is.na(z)))
  probeTable[,3] <- apply(df[,2:ncol(df)], 1, function(z) sd(z, na.rm=TRUE))
    
  write.table(probeTable, paste(writePath, "/", cancerName, "_", arrayName, "_probe_NA_sd.txt", sep = ""), append = FALSE, sep = "\t", quote = FALSE)

  samples <- read.table(paste(readPath, "/", cancerName, "_sampleinfo_", arrayName, "_to_SQL.txt", sep = ""), sep = "\t", header = TRUE)
  sampleTable[1,] <- colnames(samples)
  sampleTable[2,] <- apply(df[,2:ncol(df)], 2, function(z) sum(is.na(z)))
    
  write.table(sampleTable, paste(writePath, "/", cancerName, "_", arrayName, "_sample_NA.txt", sep = ""), append = FALSE, sep = "\t", quote = FALSE) 
}

#' Print code for PostgreSQL
#' 
#' Print the code to create the table in the PostgreSQL database (needs the file from pairBarcodes).
#' 
#' @include init
#' @usage printPostgreSQLCode(cancer, array)
#' @details The print will be displayed on the console.
printPostgreSQLCode <- function (cancer, array){
  init(cancer, array)
  samples <- read.table(paste(readPath, "/", cancerName, "_sampleinfo_", arrayName, "_to_SQL.txt", sep = ""), sep = "\t", header = TRUE)
  
  code <- paste ("CREATE TABLE ", cancerName, ".", arrayName, "( \n", sep = "")
  switch(array,
         humanmethylation450 = {
           varchar <- "10"
         },
         illuminahiseq_rnaseqv2 = {
           varchar <- "18"
         }
  )
  probe <- paste ("probe VARCHAR(", varchar, "),\n", sep = "")
  code <- paste (code, probe)
  for (i in 1:length(colnames(samples))){
    expression <- paste (colnames(samples)[i], " FLOAT(4),\n", sep = "", collapse = NULL)
    code <- paste (code, expression)
  }
  code <- paste (code, paste("CONSTRAINT pk_", cancerName, "_probe PRIMARY KEY (probe), \n", sep = ""))
  code <- paste (code, paste ("FOREIGN KEY (probe) REFERENCES ", cancerName, ".", arrayName, "_probeinfo(probe) \n );", sep = ""))
  
  return (cat(code))
}


#' Generate a TCGA barcode pairing table
#' 
#' For each tumor's humanmethylation450/ get the number of matching barcodes in agilentg4502a_07_3/, illuminaga_rnaseq/, illuminaga_rnaseqv2/, illuminahiseq_rnaseq/ and illuminahiseq_rnaseqv2/.
#' 
#' @include getLinks
#' @usage generatePairingTable()
#' @details Filters for the first 16 TCGA barcode characters, which refer to participant, sample type and vial. Duplicates are removed, along with normal and control tissue samples.
#' Finally, it will write the result in the specified path.
generatePairingTable <- function () {
  cancerNames <- getLinks("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/", "^[[:alpha:]]{2,4}/")
  
  barcodePairing <- as.data.frame (matrix(nrow = length(cancerNames), ncol = 7))
  rownames(barcodePairing) = cancerNames
  colnames(barcodePairing) <- c("humanmethylation450", "agilentg4502a_07_3", "illuminaga_rnaseq", "illuminaga_rnaseqv2", "illuminahiseq_rnaseq", "illuminahiseq_rnaseqv2", "illuminahiseq_totalrnaseqv2")
  
  for (i in 1:length(cancerNames)){
    
    tryCatch(data1 <- getBarcodes(cancerNames[[i]][[1]], "humanmethylation450"), error = function(e) {warning(e, immediate. = TRUE)})
    tryCatch(data2 <- getBarcodes(cancerNames[[i]][[1]], "agilentg4502a_07_3"), error = function(e) {warning(e, immediate. = TRUE)})
    tryCatch(data3 <- getBarcodes(cancerNames[[i]][[1]], "illuminaga_rnaseq"), error = function(e) {warning(e, immediate. = TRUE)})
    tryCatch(data4 <- getBarcodes(cancerNames[[i]][[1]], "illuminaga_rnaseqv2"), error = function(e) {warning(e, immediate. = TRUE)})
    tryCatch(data5 <- getBarcodes(cancerNames[[i]][[1]], "illuminahiseq_rnaseq"), error = function(e) {warning(e, immediate. = TRUE)})
    tryCatch(data6 <- getBarcodes(cancerNames[[i]][[1]], "illuminahiseq_rnaseqv2"), error = function(e) {warning(e, immediate. = TRUE)})
    tryCatch(data7 <- getBarcodes(cancerNames[[i]][[1]], "illuminahiseq_totalrnaseqv2"), error = function(e) {warning(e, immediate. = TRUE)})
    
    if (exists("data1", mode = "list")){
      data1$Comment..TCGA.Barcode. <- strtrim(data1$Comment..TCGA.Barcode., 16)
      data1 <- subset(data1, duplicated(Comment..TCGA.Barcode.) == FALSE)
      data1 <- subset (data1, grepl("TCGA-[[:alnum:]]*-[[:alnum:]]*-0[[:alnum:]]*", Comment..TCGA.Barcode.))
      barcodePairing[i,1] <- length(data1$Comment..TCGA.Barcode.)
      
      if (exists("data2", mode = "list")){
        data2$Comment..TCGA.Barcode. <- strtrim(data2$Comment..TCGA.Barcode., 16)
        data2 <- subset(data2, duplicated(Comment..TCGA.Barcode.) == FALSE)
        data2 <- subset (data2, Comment..TCGA.Barcode. %in% data1$Comment..TCGA.Barcode.)
        barcodePairing[i,2] <- length(data2$Comment..TCGA.Barcode.)
      }
      
      if (exists("data3", mode = "list")){
        data3$Comment..TCGA.Barcode. <- strtrim(data3$Comment..TCGA.Barcode., 16)
        data3 <- subset(data3, duplicated(Comment..TCGA.Barcode.) == FALSE)
        data3 <- subset (data3, Comment..TCGA.Barcode. %in% data1$Comment..TCGA.Barcode.)
        barcodePairing[i,3] <- length(data3$Comment..TCGA.Barcode.)
      }
      
      if (exists("data4", mode = "list")){
        data4$Comment..TCGA.Barcode. <- strtrim(data4$Comment..TCGA.Barcode., 16)
        data4 <- subset(data4, duplicated(Comment..TCGA.Barcode.) == FALSE)
        data4 <- subset (data4, Comment..TCGA.Barcode. %in% data1$Comment..TCGA.Barcode.)
        barcodePairing[i,4] <- length(data4$Comment..TCGA.Barcode.)
      }
      
      if (exists("data5", mode = "list")){
        data5$Comment..TCGA.Barcode. <- strtrim(data5$Comment..TCGA.Barcode., 16)
        data5 <- subset(data5, duplicated(Comment..TCGA.Barcode.) == FALSE)
        data5 <- subset (data5, Comment..TCGA.Barcode. %in% data1$Comment..TCGA.Barcode.)
        barcodePairing[i,5] <- length(data5$Comment..TCGA.Barcode.)
      }
      
      if (exists("data6", mode = "list")){
        data6$Comment..TCGA.Barcode. <- strtrim(data6$Comment..TCGA.Barcode., 16)
        data6 <- subset(data6, duplicated(Comment..TCGA.Barcode.) == FALSE)
        data6 <- subset (data6, Comment..TCGA.Barcode. %in% data1$Comment..TCGA.Barcode.)
        barcodePairing[i,6] <- length(data6$Comment..TCGA.Barcode.)
      }
      
      if (exists("data7", mode = "list")){
        data7$Comment..TCGA.Barcode. <- strtrim(data7$Comment..TCGA.Barcode., 16)
        data7 <- subset(data7, duplicated(Comment..TCGA.Barcode.) == FALSE)
        data7 <- subset (data7, Comment..TCGA.Barcode. %in% data1$Comment..TCGA.Barcode.)
        barcodePairing[i,7] <- length(data7$Comment..TCGA.Barcode.)
      }
    }
  }
  write.table(barcodePairing, paste(writePath, "/", "barcodePairing.txt", sep = ""), append = FALSE, sep = "\t", quote = FALSE)
}

readPath <- "~/"
writePath <- "~/"
## Collect arguments
args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
} else {  
## Evaluate arguments
  for (a in 2:length(args)){
    eval(parse(text = args[a]))
  }
}

## Help section
if("--help" %in% args) {
  cat("
      TCGA Parser help page
      
      All arguments will be evaluated as is, in the specified order:
      
      The arguments can be used to define variables and call functions.

      For further information about functions and variables see the TCGA Parser documentation.
      
      readPath and writePath are \"~/\" by default.
      
      Examples:
      \"readPath <- \\\"$PATH\\\" #Establishes the path from which files are read.\"
      \"writePath <- readPath #Establishes the path where files will be written as the readPath.\"
      \"download(\\\"coad\\\", \\\"humanmethylation450\\\") #Downloads all the colon adenocarcinoma Illumina Infinium Human DNA Methylation 450 arrays data from the TCGA.\"
      
      --help      - Display this help page
      \n")
}
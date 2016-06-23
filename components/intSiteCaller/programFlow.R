libs <- c('stringr', 'ShortRead', 'BSgenome', 'yaml')
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

# Determin the path to the python currently used by conda
python <- system('which python', intern=T)

# Load configuration file
config <<- yaml.load_file('INSPIIRED.yml')

# List to translate cluster job names to job ids.
# This will be used in future versions of the software when job waiting is reintroduced.
# This object is in the global name space because several functions will need to access it.
clusterJobIds <<- list()

writeLog <- function(...)
{
   arguments <- list(...)
  
   for(i in arguments)
   {
      if ( typeof(i) == "character" )
      {
         write(i,file=config$logFile,append=T)
      }
      else
      {
         w <- try(write.table(i,file=config$logFile,append=T,sep="\t",quote=F))
         if (class(w) == "try-error"){
           write.table("Could not write requested data item\n",file=config$logFile,append=T,sep="\t",quote=F) }
      }
   }
}

runProcess <- function(queue="normal", cpus=1, maxmem=NULL, wait=NULL, jobName=NULL, logFile=NULL, command=NULL)
{	  
   if (config$parallelize == 'bsub')
   { 
      cmd <- paste0("bsub -q ", queue, " -n ", as.character(cpus), " -M ", maxmem) 
    
      if(!is.null(wait))    cmd <- paste0(cmd, " -w \"ended(", wait, ")\"")
      if(!is.null(jobName)) cmd <- paste0(cmd, " -J \"", jobName, "\"")
      if(!is.null(logFile)) cmd <- paste0(cmd, " -o ", logFile)
    
      cmd <- paste0(cmd, " ", command)

      jobId <- system(cmd, intern=TRUE)
      writeLog(paste0('runProcess(); command: ', command))
   } else if (config$parallelize == 'qsub') {
      r <- sample(11111111:99999999, 1)

      if (is.null(jobName)) jobName <- paste0('intSite', r)
      if (is.null(logFile)) logFile <- paste0('intSite', r, '.log')
     
      write('#!/bin/bash', file=paste0(r,'.qsub'),append=F)
      write(paste0('#PBS -q ', queue),   file=paste0(r,'.qsub'),append=T)
      write(paste0('#PBS -N ', jobName), file=paste0(r,'.qsub'),append=T)

      if (!is.null(wait)) write(paste0('#PBS -W depend=afterany:', clusterJobIds[wait]), file=paste0(r,'.qsub'),append=T)

      write(paste0('#PBS -o ', logFile), file=paste0(r,'.qsub'),append=T)
      write(paste0('#PBS -e ', r, '.err'), file=paste0(r,'.qsub'),append=T)
      write(paste0('#PBS -l procs=', cpus, ',mem=', maxmem, 'mb'), file=paste0(r,'.qsub'),append=T)

      if (config$forceQsubPath) write('PATH=$PBS_O_PATH', file=paste0(r,'.qsub'),append=T)

      write('cd "$PBS_O_WORKDIR"', file=paste0(r,'.qsub'),append=T)
      write(command,file=paste0(r,'.qsub'),append=T)
      
      cmd <- paste0('qsub ', r, '.qsub')
       
      jobId <- system(cmd, intern=TRUE)
      writeLog(paste0('runProcess(); command: ', cmd))
   } else { 
      # While running serially, submit process to the system and wait for it to complete.
      writeLog(paste0('runProcess() (serial) command: ', command))
      system(command, wait=TRUE)
   }
}

#takes a textual genome identifier (ie. hg18) and turns it into the correct
#BSgenome object
get_reference_genome <- function(reference_genome) {
  pattern <- paste0("\\.", reference_genome, "$")
  match_index <- which(grepl(pattern, installed.genomes()))
  if (length(match_index) != 1) {
    write("Installed genomes are:", stderr())
    write(installed.genomes(), stderr())
    stop(paste("Cannot find unique genome for", reference_genome))
  }
  BS_genome_full_name <- installed.genomes()[match_index]
  library(BS_genome_full_name, character.only=T)
  get(BS_genome_full_name)
}

#' align sequences
#' Note: it is important not to change the blat parameters.
#' The parameters were optimized after lengthy experimentations.
#' Leave them as they are unless there is a specific reason other than
#' curoisity. Hard coded for a reason.
#' 
#' To try different blat parameters, create a file named blatOverzRide.txt
#' in the root analysis folder with the blat command template such as
#' 
#' [@node063 I1]$ cat blatOverRide.txt
#' blat %s.2bit %s %s.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead
#' [@node063 I1]$
#' 
alignSeqs <- function( dataN ){
    
    Sys.sleep(1)
   
    toAlign <- get(load("toAlign.RData"))
    alignFile <- toAlign[dataN]
    
    alias <- strsplit(alignFile, "/")[[1]][1]
    
    completeMetadata <- get(load("completeMetadata.RData"))
    genome <- completeMetadata[completeMetadata$alias==alias,"refGenome"]
    indexPath <- paste0(genome, ".2bit")
    
    blatTemplate <- "blat %s.2bit %s %s.psl -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead"
    if( file.exists("blatOverRide.txt") ) {
        blatTemplate <- readLines("blatOverRide.txt")
        writeLog("Blat parameters were overridden by file blatOverRide.txt")
    }
    cmd <-sprintf(blatTemplate, genome, alignFile, alignFile)
    unlink(paste0(alignFile, c(".psl", ".psl.gz")), force=TRUE)
    system(cmd)
    
    system(paste0("gzip ", alignFile, ".psl"))
    file.create(paste0('alignSeqs-', dataN, '.done'))
}

callIntSites <- function( dataN ){
  Sys.sleep(1)
 
  codeDir <- get(load("codeDir.RData"))
  source(file.path(codeDir, "intSiteLogic.R"))

  completeMetadata <- get(load("completeMetadata.RData"))[dataN,]
  print(t(completeMetadata), quote=FALSE)  

  status <- tryCatch(eval(as.call(append(processAlignments,
                                         unname(as.list(completeMetadata[c("alias", "minPctIdent",
                                                                           "maxAlignStart", "maxFragLength",
                                                                           "refGenome")]))))),
                     error=function(e){ writeLog(paste0('Caught error: ', e$message)) })

  save(status, file="callStatus.RData") #working directory is changed while executing getTrimmedSeqs
  file.create(paste0('../callIntSites-', dataN, '.done'))
}

demultiplex <- function(){
  I1 <- readFasta(list.files("Data", pattern="correctedI1-.", full.names=T))
  
  completeMetadata <- get(load("completeMetadata.RData"))
  
  I1 <- I1[as.vector(sread(I1)) %in% completeMetadata$bcSeq]
  samples <- completeMetadata[match(as.character(sread(I1)), completeMetadata$bcSeq), "alias"]
  
  #only necessary if using native data - can parse out description w/ python
  I1Names <-  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
  
  unlink("Data/demultiplexedReps", recursive=TRUE,  force=TRUE)
  suppressWarnings(dir.create("Data/demultiplexedReps"))

  # JKE  
  writeLog('Starting to demultiplex R1')
  R1 <- readFastq("Data/Undetermined_S0_L001_R1_001.fastq.gz")
  demultiplex_reads(R1, "R1", I1Names, samples, completeMetadata)
  writeLog('completed demultiplexing R1')  

  writeLog('Starting to demultiplex R2')
  R2 <- readFastq("Data/Undetermined_S0_L001_R2_001.fastq.gz")
  demultiplex_reads(R2, "R2", I1Names, samples, completeMetadata)
  writeLog('completed demultiplexing R2')

  file.create('demultiplex.done')
}



#' write fastq for each barcode and each sample
#' @param reads fastq reads as parsed by readFastq()
#' @param suffix either "R1" or "R2"
demultiplex_reads <- function(reads, suffix, I1Names, samples, completeMetadata) {
    RNames <- sapply(strsplit(as.character(ShortRead::id(reads)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!
    names(RNames) <- NULL

    reads <- reads[match(I1Names, RNames)]
    reads <- split(reads, samples)
    for (i in 1:length(reads)){

        writeLog(paste0('Demultiplexing ', suffix, ' read: ', i, '/', length(reads)))

        barcode.i <- completeMetadata$bcSeq[ completeMetadata$alias == names(reads)[i] ]
        stopifnot(length(barcode.i)==1)
        alias_by_barcode <- completeMetadata$alias[ completeMetadata$bcSeq == barcode.i ]
        stopifnot(length(alias_by_barcode)>=1)
        fqFiles <- paste0("Data/demultiplexedReps/", alias_by_barcode, "_", suffix, ".fastq.gz")
        cat(barcode.i, "\t", paste(fqFiles, collapse=" "), "\n" )

        null <- sapply(fqFiles, function(fq) writeFastq(reads[[i]], fq, mode="w") )
    }  
}


errorCorrectBC <- function(){
  library("ShortRead")

  codeDir <- get(load("codeDir.RData"))
  completeMetadata <- get(load("completeMetadata.RData"))
  jobID <- get(load("jobID.RData"))
  
  I1 <- readFastq("Data/Undetermined_S0_L001_I1_001.fastq.gz")
  I1 <- trimTailw(I1, 2, "0", 12)
  I1 <- I1[width(I1)==max(width(I1))]
  I1 <- split(I1, ceiling(seq_along(I1)/500000))

  for(chunk in names(I1))
  {
    writeFasta(I1[[chunk]], file=paste0("Data/trimmedI1-", chunk, ".fasta"))
  }
 
  for (i in 1:length(I1))
  { 
     runProcess(jobName=sprintf("ErrorCorrectWorker_%s-%s", jobID, i),
                maxmem=1000,
                logFile=paste0('logs/errorCorrectWorkerOutput', i, '.txt'),
                command=paste0(python, ' ', codeDir, "/errorCorrectIndices/processGolay.py ", i))
  }

  # Wait for all Golay correction jobs to be compelted.
  for (i in 1:length(I1))
  {
    repeat
    {
       if (file.exists(paste0("Data/correctedI1-", i,".done"))) break
       Sys.sleep(1)
    }
  }

  runProcess(jobName=sprintf("Demultiplex_%s", jobID),
             maxmem=50000, 
             logFile="logs/demultiplexOutput.txt",
             command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); demultiplex();\""))

  # JKE
  writeLog('Waiting for demultiplex.done ...')  
  repeat
  {
     if (file.exists('demultiplex.done')) break
     Sys.sleep(1)
  }
  writeLog('demultiplex() completed.')


  for (i in 1:nrow(completeMetadata))
  {
     runProcess(jobName=sprintf("TrimReads_%s-%s", jobID, i),
                maxmem=16000,
                logFile=paste0('logs/trimOutput', i, '.txt'),
                command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); trimReads(", i, ");\""))
  }

  for (i in 1:nrow(completeMetadata))
  {
    writeLog(paste0('Waiting for trimReads-', i, '.done ...'))
    repeat
    {
       if (file.exists(paste0('trimReads-', i, '.done'))) break
       Sys.sleep(1)
    }
  }

  writeLog('Starting postTrimReads()')
  runProcess(jobName=sprintf("PostTrimProcessing_%s", jobID),
             maxmem=8000,
             logFile="logs/postTrimOutput.txt",
             command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); postTrimReads();\""))
}


postTrimReads <- function(){

  Sys.sleep(1)
  library("BSgenome")
  library("rtracklayer") #needed for exporting genome to 2bit
  completeMetadata <- get(load("completeMetadata.RData"))

  codeDir <- get(load("codeDir.RData"))
  jobID <- get(load("jobID.RData"))
  
  numAliases <- nrow(completeMetadata)
  
  toAlign <- list.files(".", "R[12]-.*fa$", recursive=TRUE)
  toAlign <- toAlign[order(-file.info(toAlign)$size)]
  save(toAlign, file="toAlign.RData", compress=FALSE)
  numFastaFiles <- length(toAlign)

  #make temp genomes
  genomesToMake <- unique(completeMetadata$refGenome)
  
  for(genome in genomesToMake){
    export(get_reference_genome(genome), paste0(genome, ".2bit"))
  }
    
  for (i in 1:numFastaFiles)
  {
     runProcess(jobName=sprintf("AlignSeqs_%s-%s", jobID, i),
                maxmem=12000,
                logFile=paste0('logs/alignOutput', i, '.txt'),
                command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); alignSeqs(", i, ");\""))
  }
 
  for (i in 1:numFastaFiles)
  {
     writeLog(paste0('Waiting for alignSeqs-', i, '.done ...'))
     repeat
     {
       if (file.exists(paste0('alignSeqs-', i, '.done'))) break
       Sys.sleep(1)
     }
  }

  for (i in 1:nrow(completeMetadata))
  { 
      runProcess(jobName=sprintf("CallIntSites_%s-%s", jobID, i),
                 ### maxmem=120000, #multihits suck lots of memory
                 maxmem=50000,
                 logFile=paste0('logs/callSitesOutput',i,'.txt'),
                 command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); callIntSites(",i,");\""))
  }

  for (i in 1:nrow(completeMetadata))
  {
     writeLog(paste0('Waiting for callIntSites-', i, '.done ...'))
     repeat
     {
       if (file.exists(paste0('callIntSites-', i, '.done'))) break
       Sys.sleep(1)
     }
  }

  check_error()

  #writeLog('Calling check_error()')
  #runProcess(jobName=sprintf("ErrorCheck_%s", jobID),
  #           maxmem=4000,
  #           logFile="logs/errorCheck.txt",
  #           command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); check_error();\""))
}

trimReads <- function( dataN ){
    
  Sys.sleep(1)

  codeDir <- get(load("codeDir.RData"))
  source(file.path(codeDir, "intSiteLogic.R"))

  completeMetadata <- get(load("completeMetadata.RData"))[dataN,]
  
  alias <- completeMetadata$alias
  print(t(as.data.frame(completeMetadata)), quote=FALSE)
  
  suppressWarnings(dir.create(alias, recursive=TRUE))
  
  status <- tryCatch(eval(as.call(append(getTrimmedSeqs,
                                         unname(as.list(completeMetadata[c("qualityThreshold", "badQualityBases",
                                                                           "qualitySlidingWindow", "primer", "ltrBit",
                                                                           "largeLTRFrag", "linkerSequence", "linkerCommon",
                                                                           "mingDNA", "read1", "read2", "alias", "vectorSeq")]))))),
                     error=function(e){ writeLog(paste0("Caught error: ", e$message))  })
  
  writeLog("saving trimStatus.RData\n")
  save(status, file="trimStatus.RData") #working directory is changed while executing getTrimmedSeqs

  file.create(paste0('../trimReads-', dataN, '.done'))
}

processMetadata <- function(){

  jobID <- parsedArgs$jobID
  
  #expand codeDir to absolute path for saving
  codeDir <- normalizePath(parsedArgs$codeDir)

  source(file.path(codeDir, 'linker_common.R'))
  source(file.path(codeDir, 'read_sample_files.R'))

  #setting R's working dir also sets shell location for system calls, thus
  #primaryAnalysisDir is propagated without being saved

  setwd(parsedArgs$primaryAnalysisDir)

  save(jobID, file=paste0(getwd(), "/jobID.RData"))

  save(codeDir, file=paste0(getwd(), "/codeDir.RData"))

  sample_file <- 'sampleInfo.tsv'
  proc_file <- "processingParams.tsv"
  if ( ! file.exists(proc_file)) { # have to use default
      default <- "default_processingParams.tsv"
      proc_file <- file.path(codeDir, default)
  }
  completeMetadata <- read_sample_processing_files(sample_file, proc_file)
  
  completeMetadata$read1 <- paste0(getwd(), "/Data/demultiplexedReps/", completeMetadata$alias, "_R1.fastq.gz")
  completeMetadata$read2 <- paste0(getwd(), "/Data/demultiplexedReps/", completeMetadata$alias, "_R2.fastq.gz")

  stopifnot(all(c("qualityThreshold", "badQualityBases", "qualitySlidingWindow",
                  "primer", "ltrBit", "largeLTRFrag", "linkerSequence", "linkerCommon",
                  "mingDNA", "read1", "read2", "alias", "vectorSeq", "minPctIdent",
                  "maxAlignStart", "maxFragLength", "gender") %in% names(completeMetadata)))
  
  stopifnot(all( file.exists(completeMetadata$vectorSeq) ))
  
  
  ## check primer, ltrBit, largeLTRFrag consistency
  ## largeLTRFrag should start with RC(primer+ltrBit)
  rc.primer <- as.character(
      reverseComplement(DNAStringSet(completeMetadata$primer)))
  rc.ltrbit <- as.character(
      reverseComplement(DNAStringSet(completeMetadata$ltrBit)))
  
  rc.primerltrbitInlargeLTR <- mapply(function(x,y, z) grepl(y, x) | grepl(z, x),
                                      x=completeMetadata$largeLTRFrag,
                                      y=rc.primer,
                                      z=rc.ltrbit)
  
  if(!all(rc.primerltrbitInlargeLTR)) {
      print(data.frame(PLTR=names(rc.primerltrbitInlargeLTR),
                       rc.primerltrbitInlargeLTR))
      stop()
  }
  
  
  save(completeMetadata, file="completeMetadata.RData")

  suppressWarnings(dir.create("logs"))

  #error-correct barcodes - kicks off subsequent steps
  runProcess(jobName=paste0("ErrorCorrect_", jobID),
             maxmem=20000,
             logFile="logs/errorCorrectOutput.txt",
             command=paste0("Rscript -e \"source('", codeDir, "/programFlow.R'); errorCorrectBC();\""))
}


check_error <- function(errFile="error.txt") {
    writeLog(paste0("Errors if any were written to file ", errFile))
    cmd <- "grep -i \"exit\\|halt\\|huge\" logs/*.txt"
    err <- system(cmd, intern=TRUE)
    cmd <-  "grep -i max logs/*.txt | grep -i memory | awk '{print $1, $(NF-1)}' | sort -k2nr"
    mem <- system(cmd, intern=TRUE)
    if (length(err)==0) err <- "No obvious error found"
    write(c(err, "\nMemory usage", mem), errFile)

    if (!config$debug) system("rm *.qsub *.err *.done *.e *.o")
}


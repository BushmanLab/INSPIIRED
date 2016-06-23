#### check results after running intSiteCaller on simulated data ####
#' Note that reads simulated have their truth in the qnames for example
#' @M03249:1:000-SIMchr1p52699700:1:1:33:4 1:N:0:0 means
#' chr1+52699700, with length=33 down stream
#'
#' Strategy:
#' 1. read completeMetadata.RData
#' 2. load truth from Data/demultiplexedReps/GTSP0308-1_R1[2].fastq.gz
#' 3. load results from $alias/sites....
#' 4. compare alignments with different tolerance
#' 4. summarize

#' return collection of global arguments
get_args <- function() {
    suppressMessages(library(argparse))
    
    codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
    if( length(codeDir)==0 ) codeDir <- "."
    
    parser <- ArgumentParser(formatter_class='argparse.RawTextHelpFormatter')
    parser$add_argument("-c", "--codeDir", type="character", nargs=1,
                        default=codeDir,
                        help="Directory of code")
    parser$add_argument("-p", "--workDir", type="character", nargs=1,
                        default=".",
                        help="Directory of code")
    parser$add_argument("-e", "--err", type="integer", nargs=1,
                        default=5,
                        help="tolerance for alignment")
    parser$add_argument("-n", "--nproc", type="integer", nargs=1,
                        default=5,
                        help="tolerance for alignment")
    
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
    
    args$workDir <- normalizePath(args$workDir, mustWork=TRUE)
    return(args)
}
args <- get_args()
print(t(as.data.frame(args)), quote=FALSE)

libs <- c("stringr",
          "plyr",
          "dplyr",
          "RMySQL",
          "GenomicRanges",
          "ShortRead",
          "BiocParallel",
          "ggplot2")
null <- suppressMessages(sapply(libs, require, character.only=TRUE))

options(stringsAsFactors=FALSE)
options(dplyr.width = Inf)
#' increase output width to console width
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
   options(width=as.integer(howWide))
}
wideScreen()

get_metadata <- function() {
    ##df1 <- read.table(file.path(args$workDir, "sampleInfo.tsv"), header=TRUE)
    ##df2 <- read.table(file.path(args$workDir, "processingParams.tsv"), header=TRUE)
    ##df <- merge(df1, df2)
    df <-get(load("completeMetadata.RData"))
    
    ## demultiplexed R1 R2 fastq files
    df$R1fastq <- file.path(args$workDir, "Data/demultiplexedReps",
                            sprintf("%s_%s.fastq.gz", df$alias, "R1"))
    df$R2fastq <- file.path(args$workDir, "Data/demultiplexedReps",
                            sprintf("%s_%s.fastq.gz", df$alias, "R2"))
    
    ## pipeline results
    df$allSites <- file.path(args$workDir, df$alias, "allSites.RData")
    df$sites <-    file.path(args$workDir, df$alias, "sites.final.RData")
    df$rawSites <- file.path(args$workDir, df$alias, "rawSites.RData")
    df$multihit <- file.path(args$workDir, df$alias, "multihitData.RData")
    df$primerID <- file.path(args$workDir, df$alias, "primerIDData.RData")
    
    
    return(df)
}
##metadata <- get_metadata()    

#' get the original machine undemultiplexed I1, R1, R2 fastq.gz files
get_machine_file <- function(dir="Data") {
    fastqFiles <- list.files(path=dir, pattern="fastq.gz", full.names=TRUE)
    I1File <-  grep("_I1_", fastqFiles, value=TRUE)
    R1File <-  grep("_R1_", fastqFiles, value=TRUE)
    R2File <-  grep("_R2_", fastqFiles, value=TRUE)
    stopifnot(length(I1File)==1)
    stopifnot(length(R1File)==1)
    stopifnot(length(R2File)==1)
    return(data.frame(uI1=I1File, uR1=R1File, uR2=R2File) )
}
##get_machine_file(dir=file.path(args$workDir,"Data"))


#' Load truth from file
#' @param truthFile
#' @return data.frame of sample, chr, strand, position, proportion
#' @example
#' load_truth_from_file("truth.csv")
load_truth_from_file <- function(truthFile="truth.csv") {
    df <- read.csv("truth.csv")
    return(df)
}
##truth <- load_truth_from_file()




#' Load results from RData files, see intSiteUploader for help
#' @param meta dataframe of metadata
#' @return dataframe of siteID, multihitID, chr, strand, position, breakpoint, width, count qname, qid
#' @example
#' load_uniqueSites_from_RData()
#' load_multiSites_from_RData()
#' 
.load_uniqueSites_from_RData <- function(meta=metadata) {
    stopifnot(c("sites", "allSites") %in% colnames(meta))
    stopifnot(nrow(meta)==1)
    
    message("Loading\t", meta$sites, "\t", meta$allSites) 
    sites.final <- get(load(meta$sites))
    allSites <-    get(load(meta$allSites))
    
    sites <- plyr::ldply(1:length(sites.final), function(i)
        {
            i.gr <- allSites[unlist(sites.final[i]$revmap)]
            
            alias <- sites.final[i]$sampleName
            chr <- seqnames(sites.final[i])
            strand <- strand(sites.final[i])
            position <- start(flank(sites.final[i], -1, start=T))
            pcrBreakpoints <- start(flank(i.gr, -1, start=F))
            qname <- names(i.gr)
            qid <- as.integer(sub(".*:(\\d+)$", "\\1", qname))
            
            freq <- plyr::count(pcrBreakpoints)
            
            site <- data.frame(alias=as.character(alias),
                               siteID=i,
                               multihitID=0,
                               chr=as.character(chr),
                               strand=as.character(strand),
                               position=as.integer(position),
                               breakpoint=as.integer(pcrBreakpoints),
                               width=as.integer(abs(pcrBreakpoints-position)+1),
                               count=as.integer(1),
                               qname=as.character(qname),
                               qid=as.integer(qid),
                               stringsAsFactors=FALSE)
            
            return(site)
        } )       
    
    return(sites)
}
load_uniqueSites_from_RData <- function(meta=metadata) {
    sites.l <- bplapply(1:nrow(meta), function(i)
                    {
                        i.df <- try(.load_uniqueSites_from_RData(meta[i,]))
                        if( class(i.df) == "try-error" ) i.df <- data.frame()
                        return(i.df)
                    }
                        ,BPPARAM=MulticoreParam(args$nproc))
    sites <- dplyr::rbind_all(sites.l)
}   
#'
#'
#' 
.load_multiSites_from_RData <- function(meta=metadata) {
    stopifnot(c("alias", "multihit") %in% colnames(meta))
    stopifnot(nrow(meta)==1)
    message("Loading\t", meta$multihit)
    multi <- get(load(meta$multihit))
    
    multiAln <- multi$unclusteredMultihits
    
    fromCluster <- findOverlaps(multiAln,
                                multi$clusteredMultihitPositions,
                                ignore.strand=FALSE,
                                maxgap=args$err,
                                select="first")
    
    multiAln$multihitID <- fromCluster
    
    position <- start(flank(multiAln, -1, start=TRUE))
    breakpoint <- start(flank(multiAln, -1, start=FALSE))
    qname <- names(multiAln)
    qid <- as.integer(sub(".*:(\\d+)$", "\\1", qname))
    
    msite <- data.frame(
        alias=as.character(meta$alias),
        siteID=0,
        multihitID=as.integer(multiAln$multihitID),
        chr=as.character(seqnames(multiAln)),
        strand=as.character(strand(multiAln)),
        position=as.integer(position),
        breakpoint=as.integer(breakpoint),
        width=as.integer(abs(breakpoint-position)+1),
        count=as.integer(1),
        qname=as.character(qname),
        qid=as.integer(qid),
        stringsAsFactors=FALSE)
    
    return(msite)
}
load_multiSites_from_RData <- function(meta=metadata) {
    sites.l <- bplapply(1:nrow(meta), function(i)
                    {
                        i.df <- try(.load_multiSites_from_RData(meta[i,]))
                        if( class(i.df) == "try-error" ) i.df <- data.frame()
                        return(i.df)
                    }
                        ,BPPARAM=MulticoreParam(args$nproc))
    sites <- dplyr::rbind_all(sites.l)
}   


load_run_stats <- function(workDir=".") {
    
    stats.file <- list.files(workDir, pattern="^stats.RData$",
                             recursive=TRUE, full.names=TRUE)
    
    tmp.statlist <- lapply(setNames(stats.file, stats.file), function(x) {
        a <- load(x)
        get(a)
    })
    stats <- plyr:::rbind.fill(tmp.statlist)
    stats$sample <- as.character(stats$sample)
    rownames(stats) <- NULL
    
    return( c("demultiplexed"=sum(stats$barcoded),
              "LTRed"=sum(stats$LTRed),
              "linkered"=sum(stats$linkered),
              "LTRed.linkered"=sum(stats$ltredlinkered),
              "lengthTrimed"=sum(stats$lenTrim),
              "vectorTrimed"=sum(stats$vTrimed)) )
    
}


#' check total number of decoded reads
#' @param metadata, metadata
#' @return named vector of numbers of reads in R1 and R2
check_demultiplexed_reads <- function(metadata) {
    sumR1 <- summary(readFastq(metadata$R1fastq))
    sumR2 <- summary(readFastq(metadata$R2fastq))
    return(c('nR1'=as.integer(sumR1["Length"]),
             'nR2'=as.integer(sumR2["Length"])))
}


#' check run time by looking at logs
#' time start is logged in logs/errorCorrectOutput.txt
#' time end is taken as the latest from logs/callSitesOutput*.txt
#' @return time elasped in seconds
#' 
check_runtime_by_logs <- function() {
    format <- "%a %b %d %H:%M:%S %Y"
    ##x1 <- "Mon Nov  9 10:16:39 2015"
    ##x2 <- "Mon Nov  9 10:16:58 2015"
    
    cmd <- "head -15 logs/errorCorrectOutput.txt | grep Started"
    cout <- system(cmd, intern=TRUE)
    x1 <- sub("Start.*at ", "", cout)
    t1 <- strptime(x1, format=format)
    
    cmd <- "head -15 logs/callSitesOutput*.txt | grep reported | grep Results"
    cout <- system(cmd, intern=TRUE)
    x2 <- sub("Results reported on ", "", cout)
    t2 <- head(sort(strptime(x2, format=format), decreasing=TRUE), 1)
    
    tsecs <- as.integer(difftime(t2,t1, units="secs"))
    return(tsecs)
}
##check_runtime_by_logs()


#' check blat parameter from alignmentg log
#' @return time elasped in seconds
#' 
check_blat_param_by_log <- function() {
    cmd <- "grep blat logs/alignOutput1.txt | grep psl"
    cout <- system(cmd, intern=TRUE)
    
    tileSize <- as.integer(stringr::str_match(cout, "tileSize=(\\d+)")[2])
    stepSize <- as.integer(stringr::str_match(cout, "stepSize=(\\d+)")[2])
    minIdentity <- as.integer(stringr::str_match(cout, "minIdentity=(\\d+)")[2])
    maxIntron <- as.integer(stringr::str_match(cout, "maxIntron=(\\d+)")[2])
    minScore <- as.integer(stringr::str_match(cout, "minScore=(\\d+)")[2])
    
    param <- c("tileSize"=tileSize,
               "stepSize"=stepSize,
               "minIdentity"=minIdentity,
               "maxIntron"=maxIntron,
               "minScore"=minScore)
    param.name <- names(param)
    param.default <- c(11, 11, 30, 100000, 90)
    
    param <- ifelse(is.na(param), param.default, param)
    names(param) <- param.name
    
    return(param)
}
##check_blat_param_by_log()


#### load data, truth and results ####
#### site level comparison ####

metadata <- cbind(get_metadata(),
                  get_machine_file(dir=file.path(args$workDir,"Data")))

truth <- load_truth_from_file("truth.csv")
res.uniq <- load_uniqueSites_from_RData()
res.multi <- load_multiSites_from_RData()



truth.site <- subset(truth, !is.na(chr))
truth.site.gr <-  makeGRangesFromDataFrame(truth.site,
                                           start.field="position",
                                           end="position",
                                           strand.field="strand",
                                           keep.extra.columns=TRUE)


res.uniq$sample <- sub("-\\d+$", "", res.uniq$alias)
res.uniq.site <- dplyr::group_by(res.uniq, sample, alias, siteID) %>%
    dplyr::mutate(Abun=length(unique(breakpoint))) %>%
        dplyr::select(sample, alias, siteID, chr, strand, position, Abun) %>%
            dplyr::ungroup() %>% 
                dplyr::distinct()

res.uniq.site.gr <- makeGRangesFromDataFrame(res.uniq.site,
                                             start.field="position",
                                             end="position",
                                             strand.field="strand",
                                             keep.extra.columns=TRUE)

ovl <- findOverlaps(query=res.uniq.site.gr,
                    subject=truth.site.gr,
                    maxgap=args$err,
                    select="first")
res.uniq.site.gr$align <- truth.site.gr$sample[ovl]

res.uniq.site$truth <- res.uniq.site.gr$align
res.uniq.site <- plyr::arrange(res.uniq.site, sample, alias, truth, -Abun)




res.multi$sample <- sub("-\\d+$", "", res.multi$alias)
res.multi.site <- dplyr::group_by(res.multi, sample, alias, multihitID) %>%
    dplyr::mutate(Abun=length(unique(breakpoint-position))) %>%
        dplyr::select(sample, alias, multihitID, chr, strand, position, Abun) %>%
            dplyr::ungroup() %>% 
                dplyr::distinct()


res.multi.site.gr <- makeGRangesFromDataFrame(res.multi.site,
                                              start.field="position",
                                              end="position",
                                              strand.field="strand",
                                              keep.extra.columns=TRUE)

ovl <- findOverlaps(query=res.multi.site.gr,
                    subject=truth.site.gr,
                    maxgap=args$err,
                    select="first")
res.multi.site.gr$align <- truth.site.gr$sample[ovl]

res.multi.site$truth <- res.multi.site.gr$align
res.multi.site <- plyr::arrange(res.multi.site, sample, alias, truth, -Abun)

## combine multisites for a single cluster
res.multi.site.combined  <- dplyr::group_by(res.multi.site, sample, alias, multihitID) %>%
    dplyr::mutate(truth = paste(unique(truth[!is.na(truth)]), collapse=",")) %>%
        dplyr::select(sample, alias, multihitID, Abun, truth) %>%
            dplyr::ungroup() %>%
                dplyr::distinct() %>%
                    dplyr::mutate(truth=ifelse(nchar(truth)==0, NA, truth))

res.all.site <- dplyr::rbind_list(res.uniq.site, res.multi.site.combined)
res.all.site <- plyr::arrange(res.all.site, sample, alias, truth, -Abun)


cloneSum <- dplyr::group_by(res.all.site, sample) %>%
    dplyr::summarize(cloneAbun=sum(Abun[!is.na(truth)]),
                     cloneSites=length(unique(truth[!is.na(truth)])),
                     cloneAbun2Sites=cloneAbun/cloneSites,
                     fpAbun=sum(Abun[is.na(truth)]),
                     fpSites=sum(is.na(truth)),
                     fpAbun2Sites=fpAbun/fpSites)

cloneSum <- plyr::arrange(cloneSum, sample)

print(cloneSum, n=-1)

write.table(format(cloneSum, digits=2), file="cloneSummary.txt",
            sep="\t", row.names=FALSE, quote=FALSE)




#### check if uniq sites and multiHit sites overlap ####
res.multi.site.gr$siteID <- NA
res.multi.site.grl <- split(res.multi.site.gr, res.multi.site.gr$alias)

res.uniq.site.gr$multihitID <- NA
res.uniq.site.grl <- split(res.uniq.site.gr, res.uniq.site.gr$alias)

for( rep in intersect(names(res.multi.site.grl), names(res.uniq.site.grl))) {
    ##message(rep,"\t",length(res.multi.site.grl[[rep]]))
    ovl <- findOverlaps(query=res.multi.site.grl[[rep]],
                        subject=res.uniq.site.grl[[rep]],
                        maxgap=args$err,
                        select="first")
    res.multi.site.grl[[rep]]$siteID <- res.uniq.site.grl[[rep]]$siteID[ovl]
    
    ovl <- findOverlaps(query=res.uniq.site.grl[[rep]],
                        subject=res.multi.site.grl[[rep]],
                        maxgap=args$err,
                        select="first")
    res.uniq.site.grl[[rep]]$multihitID <- res.multi.site.grl[[rep]]$multihitID[ovl]
    
}




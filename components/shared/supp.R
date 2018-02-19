suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(BSgenome))
suppressMessages(library(GenomicRanges))
suppressMessages(library(colorspace))
suppressMessages(library(RSVGTipsDevice))
suppressMessages(library(foreach))
suppressMessages(library(iterators))
suppressMessages(library(rtracklayer))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomeInfoDb))

#------------- hiAnnotator -----------------------------

#' Initiate UCSC genome browser session given the freeze argument.
#'
#' @param freeze one of following: hg18, mm8, rheM, etc. Default is hg18.
#'
#' @return browser session object compatible with rtracklayer functions.
#'
#' @seealso \code{\link{getUCSCtable}}, \code{\link{makeGRanges}},
#' \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' session <- makeUCSCsession()
#' genome(session)
#' session <- makeUCSCsession("mm8")
#' genome(session)
#' }
makeUCSCsession <- function(freeze = "hg18") {
  bsession <- browserSession()
  genome(bsession) <- freeze
  bsession
}

#' Obtain a UCSC annotation table given the table & track name.
#'
#' @param tableName Name of the annotation table as it appears on UCSC browser.
#' @param trackName Name of the track annotation table as it appears in on
#' UCSC browser.
#' @param bsession UCSC session object returned by \code{\link{makeUCSCsession}}
#' or \code{\link{browserSession}}. If left NULL the function will call
#' \code{\link{makeUCSCsession}} with the provided freeze to initiate a session.
#' @param freeze one of following: hg18, mm8, rheM, etc. Default is hg18.
#' @param ... Arguments to be passed to \code{\link{ucscTableQuery}}.
#'
#' @return a dataframe containing the annotation data.
#'
#' @seealso \code{\link{makeUCSCsession}}, \code{\link{getNearestFeature}},
#' \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' refflat <- getUCSCtable("refFlat","RefSeq Genes")
#' ## same as above ##
#' refflat <- getUCSCtable("refFlat","RefSeq Genes",
#' bsession=session,freeze="hg18")
#' }
getUCSCtable <- function(tableName, trackName, bsession = NULL,
                         freeze = "hg18", ...) {
  if (is.null(bsession)) {
    bsession <- makeUCSCsession(freeze)
  }
  if (!tableName %in% tableNames(ucscTableQuery(bsession,track = trackName))) {
    stop(
      paste(
        "The provided table name:",tableName,
        "doesn't exists in track",trackName,
        "on UCSC for",freeze,"genome"
      )
    )
  }
  
  ## using getTable() instead of track() due to "No supported output types"
  ## error for certain annotation types.
  getTable(ucscTableQuery(bsession,track = trackName,table = tableName,...))
}

#' Find the column index of interest given the potential choices.
#'
#' The function finds relevant column(s) of interest from a vector of column
#' names derived from a dataframe. If no usable column is found, the function
#' spits out a relevant error or returns the index of the usable column(s).
#' This is an assistant function called by functions listed in the
#' see also section.
#'
#' @param col.names column names from a dataframe
#' @param col.options potential column names or partial names that may exist
#' in col.names
#' @param col.type type of column information the function is searching for,
#' used in construction of error messages. Default is NULL.
#' @param multiple.ok if multiple matches are found then return indices,
#' else spit an error out. Default is TRUE.
#'
#' @return the index of usable column(s) or an error if no applicable
#' column is found.
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getNearestFeature}},
#' \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' data(sites)
#' names(sites)
#' getRelevantCol(names(sites),c("chr","chromosome","tname","seqnames",
#' "chrom","contig"),"seqnames")
#' getRelevantCol(names(sites),c("ort","orientation","strand"),"strand")
getRelevantCol <- function(col.names, col.options,
                           col.type = NULL, multiple.ok = FALSE) {
  answer <- sapply(col.options,
                   function(x)
                     grep(x, col.names, ignore.case = TRUE))
  answer <- unique(as.numeric(unlist(answer)))
  if (length(answer) > 1) {
    if (!multiple.ok) {
      stop(paste(
        "More than one",col.type,"based column found:",
        paste(col.names[answer],sep = "",collapse = ", ")
      ))
    } else {
      answer
    }
  } else if (length(answer) == 0) {
    stop(paste("No",col.type,"based column found."))
  } else {
    answer
  }
}

#' Breaks two GRanges objects into chunks of N size.
#'
#' Given a query and subject GRanges objects, the function breaks query into
#' chunks of N size where each chunk has a respective subject object filtered
#' by seqnames present in the query chunk. This is a helper function used by
#' one of the annotation function in 'See Also' section where each chunk is
#' sent to a parallel node for processing.
#'
#' @param sites.rd a GRanges object.
#' @param features.rd a GRanges object.
#' @param chunkSize number of rows to use per chunk of query. Default to
#' length(sites.rd)/detectCores() or length(query)/getDoParWorkers()
#' depending on parallel backend registered.
#'
#' @return a list of GRanges objects where each element is of length 2
#' representing query & subject chunks.
#'
#' @seealso \code{\link{makeGRanges}}, \code{\link{doAnnotation}},
#' \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}},
#' \code{\link{getFeatureCounts}}.
#'
#' @export
#'
#' @examples
#' data(sites)
#' data(genes)
#' sites <- makeGRanges(sites,soloStart=TRUE)
#' genes <- makeGRanges(genes)
#' makeChunks(sites, genes)
makeChunks <- function(sites.rd, features.rd, chunkSize = NULL) {
  # do a quick check of things
  .checkArgsSetDefaults()
  rm("query","subject")
  
  # make chunks
  chunks <- breakInChunks(length(sites.rd),
                          ifelse(
                            !is.null(chunkSize),
                            length(sites.rd) / chunkSize,
                            ifelse(
                              !is.null(is.null(getDoParWorkers())),
                              length(sites.rd) / getDoParWorkers(),
                              length(sites.rd) / detectCores()
                            )
                          ))
  
  mapply(
    function(x,y) {
      new.query <- sites.rd[x:y,]
      new.query <- keepSeqlevels(new.query,
                                 value = unique(as.character(seqnames(new.query))))
      new.subject <- suppressWarnings(keepSeqlevels(features.rd,
                                                    value = seqlevels(new.query)))
      list("query" = new.query, "subject" = new.subject)
    }, start(chunks), end(chunks), SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
}

#' Clean the supplied string from punctuations and spaces.
#'
#' Function to clean the supplied string from punctuations and spaces so it
#' can be used as column headings.
#'
#' @param x string or a vector to be cleaned.
#' @param description OPTIONAL string identifying the purpose of the supplied
#' string in x to be displayed in the cleaning message. This triggers a message.
#'
#' @return cleaned string or a vector.
#'
#' @seealso \code{\link{getFeatureCounts}}, \code{\link{makeGRanges}},
#' \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' cleanColname("HIV-test")
#' cleanColname("HIV*test")
#' cleanColname("HIV-test","myAlias")
cleanColname <- function(x, description = NULL) {
  newname <- gsub("[._]+", "_", make.names(x,unique = TRUE))
  if (any(newname != x))
    if (!is.null(description))
      message("Cleaning the supplied '", description,"'")
  newname
}

#' Make a sorted GRanges object from a dataframe.
#'
#' The function converts a dataframe into a GRanges object without too much
#' hassle of renaming column names. The function finds column names that sound
#' like seqname, chromosome, start, stop, position, etc and puts them in
#' respective slots to facilitate the conversion of a dataframe to a GRanges
#' object. If more than one column that sounds like start, stop, or position
#' is present, the function will use the first match as the representative.
#' It is recommended to run this function before utilizing any other annotation functions since it will sort the object by chromosome and position for copying annotations back to their respective rows confidently.
#'
#' @param x dataframe to be converted into a GRanges object
#' @param freeze UCSC genome version of the data in x. Default is NULL.
#' This parameter is generally used to populate seqinfo slot of GRanges objects.
#' @param positionsOnly boolean flag indicating to return only position based
#' data or everything from the dataframe. Defaults to FALSE.
#' @param soloStart flag denoting whether only one position based column
#' is available. In other words, only starts are present and no stops.
#' Default=FALSE.
#' @param chromCol use the defined column name for seqname/chromosome based
#' data from the dataframe. Defaults to NULL.
#' @param strandCol use the defined column name for strand or orientation from
#' the dataframe. Defaults to NULL.
#' @param startCol use the defined column name for start coordinate from
#' the dataframe. Defaults to NULL.
#' @param stopCol use the defined column name for stop coordinate from
#' the dataframe. Defaults to NULL and not required if soloStart=TRUE.
#' @param keepFactors keep vectors/columns stored as factors? Defaults to FALSE
#'
#' @return a GRanges object converted from x.
#'
#' @seealso  \code{\link{getNearestFeature}}, \code{\link{getFeatureCounts}},
#' \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(genes)
#'
#' makeGRanges(genes, soloStart=TRUE)
#' makeGRanges(genes)
#' #makeGRanges(genes, freeze="hg18", soloStart=TRUE)
#' #makeGRanges(genes, freeze="hg18")
makeGRanges <-
  function(x, freeze = NULL, positionsOnly = FALSE, soloStart = FALSE,
           chromCol = NULL, strandCol = NULL, startCol = NULL,
           stopCol = NULL, keepFactors = FALSE) {
    ## set column names for seqnames and strand if not provided ##
    if (is.null(chromCol)) {
      colIndex <- getRelevantCol(
        names(x),
        c(
          "chr","chromosome","tname","space","chrom",
          "contig","seqnames"
        ),
        "seqnames"
      )
      chromCol <- names(x)[colIndex]
    }
    x$seqnames <- x[,chromCol]
    
    if (is.null(strandCol)) {
      colIndex <- getRelevantCol(names(x),
                                 c("ort","orientation","strand"),
                                 "strand")
      strandCol <- names(x)[colIndex]
    }
    x$strand <- x[,strandCol]
    
    if (is.null(startCol)) {
      startCol <- getRelevantCol(
        names(x),
        c("position", "intsite", "txstart",
          "start", "chromstart"),
        "start", multiple.ok = TRUE
      )
      startCol <- names(x)[startCol[1]]
    }
    x$start <- x[,startCol]
    
    ## only do stop if soloStart=F ##
    if (!as.logical(soloStart) & is.null(stopCol)) {
      stopCol <- getRelevantCol(names(x),
                                c("txend", "end", "stop", "chromend"),
                                "end", multiple.ok = TRUE)
      stopCol <- names(x)[stopCol[1]]
    }
    x$end <- x[,stopCol]
    
    ## do some testing for NAs in seqnames, start or end ##
    if (any(is.na(x$start))) {
      stop("NAs found in column containing start positions")
    }
    
    if (any(is.na(x$seqnames))) {
      stop("NAs found in column containing chromosome information")
    }
    
    ## convert any factor columns to character to avoid downstream issues ##
    if (!keepFactors) {
      factorCols <- sapply(x,is.factor)
      if (any(factorCols)) {
        for (y in names(which(factorCols))) {
          x[,y] <- as.character(x[,y])
          if (!any(is.na(suppressWarnings(as.numeric(x[,y]))))) {
            x[,y] <- as.numeric(x[,y])
          }
        }
      }
    }
    
    ## if start and end coordinates are present, sort by midpoint
    ## else if only single coordinate is present, then add  end column and sort
    if (length(startCol) > 0 & length(stopCol) > 0) {
      if (any(is.na(x$end))) {
        stop("NAs found in column containing end positions")
      }
    } else {
      x$end <- x$start
    }
    
    if (as.logical(positionsOnly)) {
      x <- x[,c("seqnames","start","end","strand")]
    }
    
    metadataCols <- setdiff(names(x),
                            c("seqnames", "start", "end", "strand", "width"))
    metadataCols <- metadataCols[!is.na(metadataCols)]
    
    sites.gr <-
      GRanges(
        seqnames = x$seqnames, IRanges(start = x$start, end = x$end),
        strand = x$strand
      )
    ## Loop through incase only one metadataCol is present which returns a
    ## vector instead of a dataframe...DataFrame(x[,metadataCols]) may not work
    for (f in metadataCols) {
      mcols(sites.gr)[[f]] <- x[,f]
    }
    
    if (!is.null(freeze)) {
      genomeLib <- grep(freeze, installed.genomes(), value = TRUE)
      if (length(genomeLib) != 0) {
        bsGenomeObject <- strsplit(genomeLib,"\\.")[[1]][2]
        chrom.info <- seqlengths(do.call(`:::`,
                                         list(genomeLib,bsGenomeObject)))
      } else {
        ## get the chromInfo file from UCSC
        z <-
          gzcon(url(
            paste0(
              "http://hgdownload.cse.ucsc.edu/goldenPath/",
              freeze, "/database/chromInfo.txt.gz"
            )
          ))
        zlines <- try(readLines(z))
        close(z)
        if (class(zlines) == "try-error")
          stop("Could not get thru to UCSC server -
               try later or drop the freeze parameter!")
        raw.data <- textConnection(zlines)
        chrom.info <- read.delim(raw.data, header = FALSE,
                                 stringsAsFactors = FALSE)[,1:2]
        chrom.info <-
          structure(chrom.info$V2, names = chrom.info$V1)
        close(raw.data)
      }
      
      # amend seqinfo slot of sites.gr #
      seqlevels(sites.gr) <- sortSeqlevels(seqlevels(sites.gr))
      seqlengths(sites.gr) <- chrom.info[seqlevels(sites.gr)]
    }
    
    sites.gr <- sort(sites.gr, ignore.strand = TRUE)
    sites.gr
  }

#' Get nearest annotation boundary for a position range.
#'
#' Given a query object, the function retrieves the nearest feature and its
#' properties from a subject and then appends them as new columns within the
#' query object. When used in genomic context, the function can be used to
#' retrieve the nearest gene 5' or 3' end relative to genomic position
#' of interest.
#'
#' @param sites.rd GRanges object to be used as the query.
#' @param features.rd GRanges object to be used as the subject or the
#' annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated
#' annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance.
#' Options are '5p','3p', 'either'(default), or 'midpoint'.
#' @param feature.colnam column name from features.rd to be used for retrieving
#' the nearest feature name. By default this is NULL assuming that features.rd
#' has a column that includes the word 'name' somewhere in it.
#' @param dists.only flag to return distances only. If this is TRUE, then
#' 'feature.colnam' is not required and only distance to the nearest feature
#' will be returned. By default this is FALSE.
#' @param parallel use parallel backend to perform calculation with
#' \code{\link[foreach]{foreach}}. Defaults to FALSE. If no parallel backend is
#' registered, then a serial version of foreach is ran using
#' \code{\link[foreach]{registerDoSEQ}}.
#' @param relativeTo calculate distance relative to query or subject.
#' Default is 'subject'. This essentially means whether to use query or subject
#' as the anchor point to get distance from!
#'
#' @return a GRanges object with new annotation columns appended at the end
#' of sites.rd.
#'
#' @note
#' \itemize{
#'   \item When side='midpoint', the distance to nearest feature is
#'   calculated by (start+stop)/2.
#'   \item If strand information doesn't exist, then everything is defaulted
#'   to '+' orientation (5' -> 3')
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered
#'   before running the function. One can use any of the following libraries
#'   compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW, doMPI,
#'   doParallel. For example: library(doMC); registerDoMC(2)
#'   \item When relativeTo="subject", the biological distance is relative to
#'   subject, meaning, the function reports the distance to query from subject
#'   (i.e. an integration site is upstream or downstream from a gene).
#'   When relativeTo="query", the distance is from the point of view of query
#'   or an integration site (i.e. gene is upstream or downstream from an
#'   integration site).
#' }
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}},
#' \code{\link{getSitesInFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene")
#' nearestGenes
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",
#' side="5p")
#' nearestGenes
#' \dontrun{
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",
#' side="3p")
#' nearestGenes
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",
#' side="midpoint")
#' ## Parallel version of getNearestFeature
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",
#' parallel=TRUE)
#' nearestGenes
#' }
getNearestFeature <- function(sites.rd, features.rd,
                              colnam = NULL, side = "either",
                              feature.colnam = NULL,
                              dists.only = FALSE, parallel = FALSE,
                              relativeTo = 'subject') {
  ## this is to avoid "no visible binding for global variable" in R CMD check
  query <- qID <- ok.chrs <- y <- freq <- NULL
  
  ## set global vars ##
  .checkArgsSetDefaults()
  
  if (!dists.only) {
    mcols(subject)$featureName <- mcols(features.rd)[,feature.colnam]
  }
  rm(features.rd)
  
  if (side %in% c('5p','3p','midpoint')) {
    ##get only 5 prime sides of features
    if (side == '5p')
      subject <- flank(subject, width = -1)
    
    ##get only 3 prime sides of features
    if (side == '3p')
      subject <- flank(subject, width = -1, start = FALSE)
    
    ##get (start+stop)/2 of features
    if (side == 'midpoint')
      ranges(subject) <- IRanges(mid(ranges(subject)), width = 1)
  }
  
  prefix <- ifelse(side == "either","",side)
  colnam <- cleanColname(colnam)
  
  ## chunksize the objects for parallel processing ##
  chunks <- if (parallel) {
    makeChunks(query, subject)
  } else {
    list(list("query" = query, "subject" = subject))
  }
  
  ## first get the nearest indices, respective tempyIDs, and distances ##
  res <- foreach(
    x = iter(chunks), .inorder = FALSE,
    .export = c("side","relativeTo"),
    .packages = c("GenomicRanges","dplyr")
  ) %dopar% {
    res.x <- as.data.frame(nearest(
      x$query, x$subject, select = "all",
      ignore.strand = TRUE
    ))
    res.x$qID <- mcols(x$query)$tempyID[res.x$queryHits]
    res.x$sID <- mcols(x$subject)$tempyID[res.x$subjectHits]
    res.x <-
      getLowestDists(x$query, x$subject, res.x, side, relativeTo)
    inner_join(res.x, count(res.x, queryHits), by = "queryHits") %>%
      rename(freq = n)
  }
  
  if (!dists.only) {
    ## for the feature of shortest indices, get the names, and strand
    ## attributes fix cases where >1 equally nearest features were returned
    ## by concatenating feature names and strand while returning one
    ## distance per query
    
    res <-
      foreach(
        x = iter(res), y = iter(sapply(chunks,"[[","subject")),
        .inorder = FALSE, .combine = rbind
      ) %dopar% {
        # make sure x & y have the respective data chunks! #
        stopifnot(all(x$sID %in% mcols(y)$tempyID))
        
        x$featureName <-
          mcols(y)[,"featureName"][x$subjectHits]
        x$strand <-
          as.character(strand(y))[x$subjectHits]
        
        # isolate non-singletons to save time & memory! #
        besties <- droplevels(subset(x,freq == 1))
        x <- droplevels(subset(x,freq > 1))
        
        x <- x %>% group_by(queryHits) %>%
          mutate(
            featureName = paste(unique(featureName),
                                collapse = ","),
            strand = paste(unique(strand), collapse =
                             ",")
          ) %>%
          ungroup %>%
          select(queryHits, qID, dist, featureName, strand) %>%
          unique
        
        # put singletons & curated non-singletons
        # back together!
        besties <- rbind(besties[,names(x)], x)
        besties <-
          arrange(besties,qID)
        
        besties
      }
    ## change column names for swift merging by .mergeAndReturn() ##
    names(res)[grepl("featureName",names(res))] <-
      paste0(prefix,colnam)
    names(res)[grepl("strand",names(res))] <-
      paste0(prefix,colnam,"Ort")
    
  } else {
    ## fix cases where >1 equally nearest features were returned by
    ## choosing 1 distance
    res <- foreach(x = iter(res), .inorder = FALSE,
                   .combine = rbind) %dopar% {
                     unique(x[,c("queryHits","qID","dist")])
                   }
  }
  
  rm(chunks)
  
  ## change distance column name for .mergeAndReturn() ##
  names(res)[grepl("dist",names(res))] <-
    paste0(prefix,colnam,"Dist")
  
  # Do a last check to make sure there is only 1 hit per qID #
  # This is useful in cases where two equally nearest distances
  # but in opposite directions are returned #
  test <- duplicated(res$qID)
  if (any(test)) {
    res <- res[!test,]
  }
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()
  
  sites.rd
}

#' Get two nearest upstream and downstream annotation boundary for a
#' position range.
#'
#' Given a query object, the function retrieves the two nearest feature
#' upstream and downstream along with their properties from a subject and
#' then appends them as new columns within the query object. When used in
#' genomic context, the function can be used to retrieve two nearest gene
#' upstream and downstream of the genomic position of interest.
#'
#' @param sites.rd GRanges object to be used as the query.
#' @param features.rd GRanges object to be used as the subject or the
#' annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated
#' annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance.
#'  Options are '5p','3p', 'either'(default), or 'midpoint'.
#' @param feature.colnam column name from features.rd to be used for retrieving
#' the nearest feature name. By default this is NULL assuming that features.rd
#' has a column that includes the word 'name' somewhere in it.
#' @param relativeTo calculate distance relative to query or subject.
#' Default is 'subject'. See documentation of  \code{\link{getNearestFeature}}
#' for more information.
#'
#' @return a GRanges object with new annotation columns appended at the end
#' of sites.rd.
#'
#' @note
#' \itemize{
#'   \item When side='midpoint', the distance to nearest feature is
#'   calculated by (start+stop)/2.
#'   \item For cases where a position is at the edge and there are no feature
#'   up/down stream since it would fall off the chromosome, the function simply
#'   returns 'NA'.
#'   \item If there are multiple locations where a query falls into,
#'   the function arbitrarily chooses one to serve as the nearest feature,
#'   then reports 2 upstream & downstream feature. That may occasionally yield
#'   features which are the same upstream and downstream, which is commonly
#'   encountered when studying spliced genes or phenomena related to it.
#'   \item If strand information doesn't exist, then everything is defaults
#'   to '+' orientation (5' -> 3')
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered
#'   before running the function. One can use any of the following libraries
#'   compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW, doMPI,
#'   doParallel. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{makeGRanges}},
#' \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene")
#' nearestGenes
#' \dontrun{
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",
#' side="5p")
#' nearestGenes
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",
#' side="3p")
#' nearestGenes
#' }
get2NearestFeature <- function(sites.rd, features.rd,
                               colnam = NULL, side = "either",
                               feature.colnam = NULL, relativeTo = "subject") {
  ## this is to avoid "no visible binding for global variable" in R CMD check
  query <- qID <- ok.chrs <- y <- freq <- NULL
  
  ## set global vars ##
  .checkArgsSetDefaults()
  
  ## make sure features.rd/subject is sorted ##
  mcols(subject)$featureName <- mcols(features.rd)[,feature.colnam]
  subject <- sort(subject)
  rm(features.rd)
  
  if (side %in% c('5p','3p','midpoint')) {
    ##get only 5 prime sides of features
    if (side == '5p')
      subject <- flank(subject, width = -1)
    
    ##get only 3 prime sides of features
    if (side == '3p')
      subject <- flank(subject, width = -1, start = FALSE)
    
    ##get (start+stop)/2 of features
    if (side == 'midpoint')
      ranges(subject) <- IRanges(mid(ranges(subject)), width = 1)
  }
  
  ## u = upstream, d = downstream
  ## thinking concept: u2.....u1.....intSite(+).....d1.....d2
  ## thinking concept: d2.....d1.....intSite(-).....u1.....u2
  ## searching concept: res.left2.....res.left1.....res....intSite....res...
  ## ..res.right1.....res.right2
  
  ## first get the nearest indices, respective tempyIDs ##
  res <- as.data.frame(nearest(query, subject, select = "all",
                               ignore.strand = TRUE))
  res$qID <- mcols(query)$tempyID[res$queryHits]
  res$qStrand <- as.character(strand(query))[res$queryHits]
  res <- getLowestDists(query, subject, res, side, relativeTo)
  
  ## perform upstream-downstream checks by testing distances
  res$u2 <- with(res, 
                 ifelse(dist < 0,
                        ifelse(qStrand == "+", subjectHits - 1, 
                               subjectHits + 1),
                        ifelse(qStrand == "+", subjectHits - 2, 
                               subjectHits + 2)
                 ))
  
  res$u1 <- with(res, 
                 ifelse(dist < 0, subjectHits,
                        ifelse(qStrand == "+", subjectHits - 1, subjectHits + 1)
                 ))
  
  res$d1 <- with(res, 
                 ifelse(dist < 0,
                        ifelse(qStrand == "+", subjectHits + 1, 
                               subjectHits - 1), 
                        subjectHits
                 ))
  
  res$d2 <- with(res, 
                 ifelse(dist < 0,
                        ifelse(qStrand == "+", subjectHits + 2, 
                               subjectHits - 2),
                        ifelse(qStrand == "+", subjectHits + 1, 
                               subjectHits - 1)
                 ))
  
  prefix <- ifelse(side == "either","Either",side)
  
  message("u = upstream, d = downstream")
  message("thinking concept: u2.....u1.....intSite(+).....d1.....d2")
  message("thinking concept: d2.....d1.....intSite(-).....u1.....u2")
  
  colnam <- cleanColname(colnam)
  
  ## add columns back to query object
  all.res <- lapply(c("u1","u2","d1","d2"), function(f) {
    message(f)
    res.nrst <- res[,c("queryHits","subjectHits","qID",f)]
    
    # make sure we haven't jumped a chromosome by shifting nearest indices
    # if we did, then record it as NA at later a stage
    
    # fix cases where chosen subjectHits are off the
    # length of subjectObject
    fixed <- with(res.nrst, ifelse(get(f) < 1 | get(f) > length(subject),
                                   subjectHits, get(f)))
    
    # do the chromosome test & tag rows which were off the subject length #
    res.nrst$qChr <- as.character(seqnames(query))[res.nrst$queryHits]
    res.nrst$sChr <- as.character(seqnames(subject))[fixed]
    rows <- res.nrst$qChr != res.nrst$sChr | res.nrst[,f] < 1 |
      res.nrst[,f] > length(subject)
    
    # overwrite subjectHits indices with that of interested motif for l
    # ater steps
    res.nrst$subjectHits <- res.nrst[,f]
    
    # remove unnecessary columns #
    res.nrst[,f] <- NULL
    res.nrst$qChr <- NULL
    res.nrst$sChr <- NULL
    
    # extract cases which fell off the chromosome but drop any queryHits
    # which found multiple nearest hits and only one of them happened to be
    # off the chromosome!
    res.nrst.bad <- droplevels(res.nrst[rows & !res.nrst$queryHits %in%
                                          res.nrst$queryHits[!rows],])
    res.nrst <- droplevels(res.nrst[!rows,])
    
    res.nrst <- getLowestDists(query, subject, res.nrst, side, relativeTo)
    res.nrst$featureName <- mcols(subject)[,"featureName"][res.nrst$subjectHits]
    res.nrst$strand <- as.character(strand(subject))[res.nrst$subjectHits]
    
    res.nrst <- res.nrst %>% group_by(queryHits, qID, dist) %>%
      summarise(featureName = paste(unique(featureName), collapse = ","),
                strand = paste(unique(strand), collapse = ",")) %>% ungroup
    
    # add back rows which fell off the edge of chromosome #
    if (any(rows)) {
      res.nrst.bad[,c("dist", "featureName", "strand")] <- NA
      res.nrst <- rbind(res.nrst, unique(res.nrst.bad[,names(res.nrst)]))
      res.nrst <- arrange(res.nrst, qID)
    }
    
    if (f == "u1") {
      coldef <- paste(prefix,colnam,"upStream1",sep = ".")
    }
    if (f == "u2") {
      coldef <- paste(prefix,colnam,"upStream2",sep = ".")
    }
    if (f == "d1") {
      coldef <- paste(prefix,colnam,"downStream1",sep = ".")
    }
    if (f == "d2") {
      coldef <- paste(prefix,colnam,"downStream2",sep = ".")
    }
    
    ## add meta columns to the result ##
    names(res.nrst)[grepl("featureName",names(res.nrst))] <- coldef
    names(res.nrst)[grepl("strand",names(res.nrst))] <-
      paste(coldef,"Ort", sep = ".")
    names(res.nrst)[grepl("dist",names(res.nrst))] <- 
      paste(coldef,"Dist", sep = ".")
    
    res.nrst
  })
  
  res <- do.call(cbind, all.res)
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()
  
  sites.rd
}

#' Get the lowest biological distance from the 5' or 3' boundaries of query
#' and subject.
#'
#' Given a query and subject with indicies from \code{\link[IRanges]{nearest}},
#' calculate the shortest biological distance to either boundaries of the query
#' and subject. This is a helper function utilized in
#' \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}
#'
#' @param query GRanges object to be used as the query which holds data for
#' 'queryHits' attribute of res.nrst.
#' @param subject GRanges object to be used as the subject which holds data for
#' 'subjectHits' attribute of res.nrst.
#' @param res.nrst a dataframe of nearest indices as returned by
#' \code{\link[IRanges]{nearest}}.
#' @param side boundary of subject/annotation to use to calculate the
#' nearest distance. Options are '5p','3p', or the default 'either'.
#' @param relativeTo calculate distance relative to query or subject.
#' Default is 'subject'. See documentation of  \code{\link{getNearestFeature}}
#' for more information.
#'
#' @return res.nrst with lowest distances appended at the end.
#'
#' @note for cases where a query has multiple nearest neighbors or overlaps
#' with >1 subjects, the function will choose the subject with the lowest
#' absolute distance.
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @export
#'
#' @examples
#' query <- GRanges("A", IRanges(c(1, 5, 12, 20), width=1),
#' strand=c("-","+","-","+"))
#' subject <- GRanges("A", IRanges(c(1,5,10,15,21), width=8:4),
#' strand=c("+", "+", "-", "-","-"))
#' res <- as.data.frame(nearest(query, subject, select="all",
#' ignore.strand=TRUE))
#' res <- getLowestDists(query, subject, res, "either", "query")
#'
getLowestDists <- function(query = NULL, subject = NULL, res.nrst = NULL,
                           side = "either", relativeTo = "subject") {
    if (is.null(query) | is.null(subject) | is.null(res.nrst)) {
      stop("One of following is null: query, subject, res.nrst")
    }
    
    if (side == "either") {
      ## get lowest dist to either annot boundary from 5p side of the query
      dist.s <- start(query)[res.nrst$queryHits] -
        start(subject)[res.nrst$subjectHits]
      dist.e <- start(query)[res.nrst$queryHits] -
        end(subject)[res.nrst$subjectHits]
      dist5p <- ifelse(abs(dist.s) < abs(dist.e), dist.s, dist.e)
      
      ## get lowest dist to either annot boundary from 3p side of the query
      dist.s <- end(query)[res.nrst$queryHits] -
        start(subject)[res.nrst$subjectHits]
      dist.e <- end(query)[res.nrst$queryHits] -
        end(subject)[res.nrst$subjectHits]
      dist3p <- ifelse(abs(dist.s) < abs(dist.e), dist.s, dist.e)
    } else {
      ## no need to do calcs to start & end of subject since this clause
      ## assumes you have taken 5' or 3' of the subject!
      ## get lowest dist to annot boundary from 5p side of the query
      dist5p <- start(query)[res.nrst$queryHits] -
        start(subject)[res.nrst$subjectHits]
      
      ## get lowest dist to annot boundary from 3p side of the query
      dist3p <- end(query)[res.nrst$queryHits] -
        start(subject)[res.nrst$subjectHits]
    }
    
    ## get lowest distance from the lowest 5p or 3p of the query!
    dist.lowest <- ifelse(abs(dist5p) < abs(dist3p), dist5p, dist3p)
    
    ## fix signs to match biological upstream or downstream relative to
    ## query or subject!
    if (relativeTo == 'query') {
      bore <- as.character(strand(query))[res.nrst$query] == "+"
      dist.lowest2 <- ifelse(bore,-dist.lowest, dist.lowest)
    } else {
      bore <- as.character(strand(subject))[res.nrst$subjectHits] == "-"
      dist.lowest2 <- ifelse(bore,-dist.lowest, dist.lowest)
    }
    rm(bore)
    
    res.nrst$dist <- dist.lowest2
    
    ## fix cases where two nested features were returned by choosing
    ## the lowest absolute distances for both features.
    mins <- with(res.nrst, tapply(abs(dist), queryHits, min))
    res.nrst$lowest <- with(res.nrst, abs(dist) == mins[as.character(queryHits)])
    res.nrst <- droplevels(res.nrst[res.nrst$lowest,])
    res.nrst$lowest <- NULL
    
    res.nrst
  }

#' Generate a window size label.
#'
#' Function to generate aesthetically pleasing window size label given an
#' integer. This is one of the helper function used in
#' \code{\link{getFeatureCounts}} & \code{\link{getFeatureCountsBig}}.
#'
#' @param x vector of integers to generate the labels for.
#'
#' @return a character vector of length(x) which has x normalized and
#' suffixed by bp, Kb, Mb, or Gb depending on respective interval sizes.
#'
#' @seealso \code{\link{getFeatureCounts}}, \code{\link{makeGRanges}},
#' \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' getWindowLabel(c(0,1e7,1e3,1e6,2e9))
getWindowLabel <- function(x) {
  ind <- cut(abs(x), c(0, 1e3, 1e6, 1e9, 1e12),
             include.lowest = TRUE, right = FALSE, labels = FALSE)
  paste0(x / c(1, 1e3, 1e6, 1e9, 1e12)[ind],
         c("bp", "Kb", "Mb", "Gb")[ind])
}

#' Get counts of annotation within a defined window around each query
#' range or positions.
#'
#' Given a query object and window size(s), the function finds all the rows in
#' subject which are <= window size/2 distance away. If weights are assigned to
#' each positions in the subject, then tallied counts are multiplied
#' accordingly. For large annotations, use \code{\link{getFeatureCountsBig}}.
#'
#' @param sites.rd GRanges object to be used as the query.
#' @param features.rd GRanges object to be used as the subject or
#' the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated
#' annotation...serves as a prefix to windows sizes!
#' @param chromSizes named vector of chromosome/seqnames sizes to be used for
#' testing if a position is off the mappable region. DEPRECATED and will be
#' removed in future release.
#' @param widths a named/numeric vector of window sizes to be used for casting
#' a net around each position. Default: \code{c(1000,10000,1000000)}.
#' @param weightsColname if defined, weigh each row from features.rd when
#' tallying up the counts.
#' @param doInChunks break up sites.rd into small pieces of chunkSize to
#' perform the calculations. Default is FALSE. Useful if you are expecting
#' to find great deal of overlap between sites.rd and features.rd.
#' @param chunkSize number of rows to use per chunk of sites.rd.
#' Default to 10000. Only used if doInChunks=TRUE.
#' @param parallel use parallel backend to perform calculation with
#' \code{\link[foreach]{foreach}}. Defaults to FALSE. If no parallel backend is
#'  registered, then a serial version of foreach is ran using
#'  \code{\link[foreach]{registerDoSEQ}}.
#'
#' @return a GRanges object with new annotation columns appended at the end of
#' sites.rd. There will be a column for each width defined in widths parameter.
#' If widths was a named vector i.e. c("100bp"=100,"1K"=1000), then the colname
#' parameter will be pasted together with width name else default name will be
#' generated by the function.
#'
#' @note
#' \itemize{
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered
#'    before running the function. One can use any of the following libraries
#'    compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW,
#'    doMPI. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getNearestFeature}},
#' \code{\link{getSitesInFeature}}, \code{\link{getFeatureCountsBig}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene")
#' \dontrun{
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",
#' doInChunks=TRUE, chunkSize=200)
#' geneCounts
#' ## Parallel version of getFeatureCounts
#' # geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",
#' parallel=TRUE)
#' # geneCounts
#' }
getFeatureCounts <- function(sites.rd, features.rd,
                             colnam = NULL, chromSizes = NULL,
                             widths = c(1000,10000,1000000), 
                             weightsColname = NULL, doInChunks = FALSE, 
                             chunkSize = 10000, parallel = FALSE) {
  
  ## this is to avoid "no visible binding for global variable" in R CMD check
  query <- qID <- ok.chrs <- y <- freq <- NULL
  
  ## set global vars ##
  .checkArgsSetDefaults()
  
  if (!is.null(chromSizes)) {
    warning("decrepit option: chromSizes parameter is no longer required
            and will be ignored!")
  }
  
  if (doInChunks & chunkSize < length(sites.rd)) {
    rm("query","subject")
    
    # no need to execute all this if chunkSize is bigger than data size!!!
    total <- length(sites.rd)
    starts <- seq(1,total,by = chunkSize)
    stops <- unique(c(seq(chunkSize, total, by = chunkSize), total))
    stopifnot(length(starts) == length(stops))
    message("Breaking up sites.rd into chunks of ",chunkSize)
    res <- GRanges()
    for (x in 1:length(starts)) {
      res <- c(res,
               as(
                 getFeatureCounts(sites.rd[starts[x]:stops[x],],
                                  features.rd, colnam = colnam,
                                  widths = widths,
                                  weightsColname = weightsColname,
                                  parallel = parallel), 
                 "GRanges"
                 )
               )
    }
    
    return(res)
  } else {
    weighted <- ifelse(is.null(weightsColname),FALSE,TRUE)
    if (weighted) {
      mcols(subject)$weights <- mcols(features.rd)[,weightsColname]
    }
    rm(features.rd)
    
    # only get labels if not supplied
    if (is.null(names(widths))) {
      names(widths) <- getWindowLabel(widths)
    }
    
    ## chunkize the objects for parallel processing ##
    chunks <- if (parallel) {
      makeChunks(query, subject)
    } else {
      list(list("query" = query, "subject" = subject))
    }
    
    colnam <- cleanColname(colnam)
    
    ## perform overlap analysis in parallel by windows ##
    res <- foreach(
      x = iter(chunks), .inorder = FALSE, .combine = rbind,
      .export = c("widths","weighted","colnam"),
      .packages = c("GenomicRanges","dplyr")
    ) %dopar% {
      counts <- sapply(widths, function(y) {
        res.x <- findOverlaps(x$query, x$subject, select = 'all',
                              maxgap = (y/2), ignore.strand = TRUE)
        
        if (weighted) {
          res.x <- as.data.frame(res.x)
          res.x$weights <- mcols(x$subject)$weights[res.x$subjectHits]
          res.x$tempyID <- mcols(x$query)$tempyID[res.x$queryHits]
          res.x$counts <- with(res.x, ave(weights,tempyID,FUN = sum))
          res.x <- unique(res.x[,c("tempyID", "counts")])
          
          nohits <- setdiff(mcols(x$query)$tempyID, res.x$tempyID)
          res.x <- rbind(res.x, data.frame(tempyID = nohits, counts = 0))
          res.x[order(res.x$tempyID),"counts"]
        } else {
          countQueryHits(res.x)
        }
      })
      counts <- as.data.frame(counts)
      names(counts) <- paste(colnam,names(counts),sep = ".")
      counts$qID <- mcols(x$query)$tempyID
      counts
    }
    
    rm(chunks)
    
    ## merge results to the query object and return it ##
    .mergeAndReturn()
    
    sites.rd
  }
}

#' Get counts of annotation within a defined window around each query
#' range/position for large annotation objects spanning greater than
#' 1 billion rows.
#'
#' Given a query object and window size(s), the function finds all the rows in
#' subject which are <= window size/2 distance away. Note that here counting is
#' done using midpoint of the ranges in query instead of start-stop boundaries.
#' The counts will differ slightly when compared to
#' \code{\link{getFeatureCounts}}.
#'
#' @param sites.rd GRanges object to be used as the query.
#' @param features.rd GRanges object to be used as the subject or the
#' annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated
#' annotation...serves as a prefix to windows sizes!
#' @param widths a named/numeric vector of window sizes to be used for casting
#' a net around each position. Default: \code{c(1000,10000,1000000)}
#'
#' @return a GRanges object with new annotation columns appended at the end of
#' sites.rd. There will be a column for each width defined in widths parameter.
#' If widths was a named vector i.e. c("100bp"=100,"1K"=1000),
#' then the colname parameter will be pasted together with width name else
#' default name will be generated by the function.
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getNearestFeature}},
#' \code{\link{getSitesInFeature}}, \code{\link{getFeatureCounts}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' geneCounts1 <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene")
#' \dontrun{
#' geneCounts2 <- getFeatureCountsBig(alldata.rd,genes.rd,"NumOfGene")
#' identical(geneCounts1,geneCounts2)
#' }
getFeatureCountsBig <- function(sites.rd, features.rd,colnam = NULL, 
                                widths = c(1000,10000,1000000)) {
  
  ## this is to avoid "no visible binding for global variable" in R CMD check
  query <- qID <- ok.chrs <- y <- freq <- NULL
  
  ## set global vars ##
  .checkArgsSetDefaults()
  rm(features.rd)
  
  ranges(query) <- mid(ranges(query))
  query <- split(query, seqnames(query))
  subject <- split(subject,seqnames(subject))
  
  # only get labels if not supplied
  if (is.null(names(widths))) {
    names(widths) <- getWindowLabel(widths)
  }
  
  colnam <- cleanColname(colnam)
  
  ## get counts of midpoints using findInterval and add columns back to
  ## query object
  res <- lapply(names(widths), function(windowName) {
    message(".")
    columnName <- paste(colnam,names(widths[windowName]),sep = ".")
    
    res.i <- lapply(ok.chrs, function(x) {
      counts <- abs(findInterval(start(query[[x]]) - widths[windowName] / 2,
                                 sort(start(subject[[x]]))) -
                      findInterval(start(query[[x]]) + widths[windowName] / 2,
                                   sort(end(subject[[x]])))
      )
      res.x <- data.frame(qID = query[[x]]$tempyID)
      res.x[,columnName] <- counts
      res.x
    })
    
    do.call(rbind,res.i)
  })
  
  res <- do.call(cbind,res)
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()
  
  sites.rd
}

#' Find overlapping positions/ranges that match between the query and subject.
#'
#' When used in genomic context, the function annotates genomic positions of
#' interest with information like if they were in a gene or cpg island or
#' whatever annotation that was supplied in the subject.
#'
#' @param sites.rd GRanges object to be used as the query.
#' @param features.rd GRanges object to be used as the subject or the
#' annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated
#' annotation...serves a core! If allSubjectCols=TRUE, then this is used as a
#' prefix to all metadata column.
#' @param asBool Flag indicating whether to return results as TRUE/FALSE or the
#' property of an overlapping feature..namely feature name and orientation if
#' available. Defaults to FALSE.
#' @param feature.colnam column name from features.rd to be used for retrieving
#' the feature name. By default this is NULL assuming that features.rd has a
#' column that includes the word 'name' somewhere in it.
#' Not required if asBool=TRUE or allSubjectCols=TRUE
#' @param parallel use parallel backend to perform calculation with
#' \code{\link[foreach]{foreach}}. Defaults to FALSE. Not applicable when
#' asBool=T. If no parallel backend is registered, then a serial version of
#' foreach is ran using \code{\link[foreach]{registerDoSEQ}}.
#' @param allSubjectCols Flag indicating whether to return all annotations or
#' metadata columns from features.rd. Defaults to FALSE.
#' @param overlapType see \code{\link[IRanges]{findOverlaps}}. Defaults to 'any'
#'
#' @return a GRanges object with new annotation columns appended at the end
#' of sites.rd.
#'
#' @note
#' \itemize{
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered
#'   before running the function. One can use any of the following libraries
#'   compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW,
#'   doMPI. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}},
#' \code{\link{getNearestFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene")
#' InGenes
#' \dontrun{
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",asBool=TRUE)
#' InGenes
#' ## Parallel version of getSitesInFeature
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",asBool=TRUE,
#' parallel=TRUE)
#' InGenes
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",
#' allSubjectCols=TRUE, parallel=TRUE)
#' InGenes
#' }
getSitesInFeature <- function(sites.rd, features.rd, colnam = NULL,
                              asBool = FALSE, feature.colnam = NULL,
                              parallel = FALSE, allSubjectCols = FALSE,
                              overlapType = 'any') {
  ## this is to avoid "no visible binding for global variable" in R CMD check
  query <- qID <- ok.chrs <- y <- freq <- NULL
  
  ## set global vars ##
  .checkArgsSetDefaults()
  
  ## chunkize the objects for parallel processing ##
  mcols(subject)$featureName <- mcols(features.rd)[,feature.colnam]
  rm(features.rd)
  
  chunks <- if (parallel) {
    makeChunks(query, subject)
  } else {
    list(list("query" = query, "subject" = subject))
  }
  
  colnam <- cleanColname(colnam)
  
  ## perform overlap analysis in parallel by windows ##
  res <- foreach(
    x = iter(chunks), .inorder = FALSE, .combine = rbind,
    .export = c("colnam","asBool","allSubjectCols"),
    .packages = c("GenomicRanges","dplyr")
  ) %dopar% {
    if (asBool) {
      strand(x$subject) <- "*"
      bore <- overlapsAny(x$query, x$subject, ignore.strand = TRUE,
                          type = overlapType)
      res.x <- data.frame(qID = mcols(x$query)$tempyID, featureName = bore)
    } else if (allSubjectCols) {
      # remove artificially added featureName column else it will be duplicated! 
      mcols(x$subject)$featureName <- NULL 
      res.x <- as.data.frame(findOverlaps(x$query, x$subject, 
                                          select = 'all', ignore.strand = TRUE,
                                          type = overlapType))
      res.x$qID <- mcols(x$query)$tempyID[res.x$queryHits]
      allSubjCols <- mcols(x$subject[res.x$subjectHits])
      rownames(allSubjCols) <- NULL
      names(allSubjCols) <- paste("featureName", names(allSubjCols), sep = ".")
      res.x <- cbind(res.x, as.data.frame(allSubjCols))
    } else {
      res.x <- as.data.frame(
        findOverlaps(
          x$query, x$subject, select = 'all', ignore.strand = TRUE,
          type = overlapType
        )
      )
      res.x$qID <- mcols(x$query)$tempyID[res.x$queryHits]
      
      ## collapse rows where query returned two hits with the same featureNames
      ## due to alternative splicing or something else.
      res.x$featureName <- mcols(x$subject)$featureName[res.x$subjectHits]
      res.x$strand <- as.character(strand(x$subject))[res.x$subjectHits]
      
      # isolate non-singletons to save time & memory! #
      res.x <- res.x %>% group_by(queryHits) %>% mutate(freq = n())
      besties <- ungroup(res.x) %>% filter(freq == 1)
      res.x <- ungroup(res.x) %>% filter(freq > 1)
      
      # collapse multiple featureName #
      res.x <- res.x %>% group_by(queryHits) %>%
        mutate(featureName = paste(unique(featureName), collapse = ","),
               strand = paste(unique(strand), collapse = ",")) %>%
        ungroup %>% select(queryHits, qID, featureName, strand) %>%
        unique
      
      # put singletons & curated non-singletons back together #
      res.x <- rbind(besties[,names(res.x)], res.x)
      rm(besties)
      
      res.x <- arrange(res.x, qID)
      
      names(res.x)[grepl("strand",names(res.x))] <- paste0(colnam,"Ort")
    }
    
    ## change column names for swift merging by
    # .mergeAndReturn()
    if (allSubjectCols) {
      names(res.x)[grepl("featureName", names(res.x))] <-
        gsub("featureName", colnam,
             names(res.x)[grepl("featureName", names(res.x))])
    } else {
      names(res.x)[grepl("featureName", names(res.x))] <- colnam
    }
    
    res.x
  }
  
  rm(chunks)
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()
  
  ## for legacy code support change sites.rd not in feature.rd to FALSE
  ## instead of NA
  if (!allSubjectCols) {
    mcols(sites.rd)[,colnam][is.na(mcols(sites.rd)[,colnam])] <- FALSE
  }
  
  sites.rd
}

#' Annotate a GRanges object using one of annotation functions.
#'
#' This is a wrapper function which calls one of following functions depending
#' on annotType parameter: \code{\link{getFeatureCounts}},
#' \code{\link{getFeatureCountsBig}}, \code{\link{getNearestFeature}},
#' \code{\link{get2NearestFeature}}, \code{\link{getSitesInFeature}}
#'
#' @param annotType one of following: within, nearest, twoNearest, counts,
#' countsBig.
#' @param ... Additional parameters to be passed to the respective annotation
#' function.
#' @param postProcessFun function to call on the resulting object for any post
#' processing steps.
#' @param postProcessFunArgs additional arguments for postProcessFun as a list.
#'
#' @return a GRanges object with new annotation columns appended at the end
#' of sites.rd.
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}},
#' \code{\link{getFeatureCountsBig}}, \code{\link{getNearestFeature}},
#' \code{\link{get2NearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene",asBool=TRUE)
#' \dontrun{
#' doAnnotation(annotType="counts",alldata.rd,genes.rd,"NumOfGene")
#' doAnnotation(annotType="nearest",alldata.rd,genes.rd,"NearestGene")
#' doAnnotation(annotType="countsBig",alldata.rd,genes.rd,"ChipSeqCounts")
#' geneCheck <- function(x,wanted) { x$isWantedGene <- x$InGene %in% wanted;
#' return(x) }
#' doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene",
#' postProcessFun=geneCheck,
#' postProcessFunArgs=list("wanted"=c("FOXJ3","SEPT9","RPTOR")) )
#' }
doAnnotation <-
  function(annotType = NULL, ..., postProcessFun = NULL,
           postProcessFunArgs = list()) {
    if (is.null(annotType)) {
      stop(
        "Please define the annotType parameter to identify which type of
        annotation to perform: within, nearest, counts"
      )
    }
    
    res <- switch(
      match.arg(annotType, c("within", "nearest", "twoNearest",
                             "counts", "countsBig")),
      within = getSitesInFeature(...),
      nearest = getNearestFeature(...),
      twoNearest = get2NearestFeature(...),
      counts = getFeatureCounts(...),
      countsBig = getFeatureCountsBig(...),
      stop("Invalid annoType parameter")
    )
    
    if (!is.null(postProcessFun)) {
      res <- do.call(postProcessFun, append(postProcessFunArgs,
                                            list(res), after = 0))
    }
    
    res
    }

#' Check args and set defaults.
#'
#' This function checks all the arguments passed to an annotation
#' function and set default values for later use. Evaluation of this function
#' happens in the parent function.
#'
.checkArgsSetDefaults <- function() {
  checks <- expression(
    if (!is(sites.rd,"GRanges")) {
      sites.rd <- as(sites.rd,"GRanges")
    },
    
    if (!is(features.rd,"GRanges")) {
      features.rd <- as(features.rd,"GRanges")
    },
    
    if (!identical(class(sites.rd),class(features.rd))) {
      stop(
        "sites.rd & features.rd are of different classes.
        Please make them the same class: GRanges"
      )
    },
    
    stopifnot(length(sites.rd) > 0),
    stopifnot(length(features.rd) > 0),
    
    if ("parallel" %in% names(formals())) {
      if (!parallel) {
        registerDoSEQ()
      }
    },
    
    if (exists("colnam")) {
      if (is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s)
             to be appended.")
      }
    },
    
    if (!any(unique(as.character(seqnames(sites.rd))) %in%
             unique(as.character(seqnames(features.rd))))) {
      stop(
        "There are no seqnames/chromosomes that are shared between the
         query (sites.rd) and subject (features.rd)"
      )
    },
    
    ## get feature names column for adding feature name to sites.rd
    ## dont throw error if dists.only flag is TRUE from getNearestFeature
    if (exists("feature.colnam")) {
      if (is.null(feature.colnam)) {
        answer <- try(getRelevantCol(colnames(mcols(features.rd)),
                                     c("name","featureName"),
                                     "featureName",multiple.ok = TRUE),
                      silent = TRUE)
        feature.colnam <- colnames(mcols(features.rd))[answer][1]
      }
      
      if (exists("dists.only")) {
        if (!dists.only & is.na(feature.colnam)) {
          stop("No featureName based column found.")
        }
      } else if (exists("allSubjectCols")) {
        if (!allSubjectCols & is.na(feature.colnam)) {
          stop("No featureName based column found.")
        }
      } else {
        if (is.na(feature.colnam)) {
          stop("No featureName based column found.")
        }
      }
    },
    
    ## check strand column of the feature ##
    if (all(strand(features.rd) == "*")) {
      message("Setting strand to '+' for features.rd parameter")
      strand(features.rd) <- "+"
    },
    
    ## check strand column of the feature ##
    if (all(strand(sites.rd) == "*")) {
      message("Setting strand to '+' for sites.rd parameter")
      strand(sites.rd) <- "+"
    },
    
    ## use only chroms that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(as.character(seqnames(sites.rd)),
                         as.character(seqnames(features.rd))),
    features.rd <- keepSeqlevels(features.rd, ok.chrs),
    good.rows <- as.character(seqnames(sites.rd)) %in% ok.chrs,
    
    ## extract required objects to streamline downstream code/steps
    ## tag each row with tempyID for merging with original object
    ## this is crucial since objects are divided into chunks which resets
    ## the index from 1...n
    ## tempyID would preserve the original order for parallel processing ##
    query <- sites.rd,
    mcols(query) <- NULL,
    mcols(query)$tempyID <- 1:length(query),
    mcols(sites.rd)$tempyID <- mcols(query)$tempyID,
    
    subject <- features.rd,
    if (exists("allSubjectCols")) {
      if (!allSubjectCols) {
        mcols(subject) <- NULL
        mcols(subject)$tempyID <- 1:length(subject)
      }
    } else {
      mcols(subject) <- NULL
      mcols(subject)$tempyID <- 1:length(subject)
    }
    
  )
  
  eval.parent(checks)
}

#' Merge results back to the query object and perform additional post
#' processing steps.
#'
#' This function merges all the calculation results back to the query object.
#' Additionally, if any flags were set, the function does the necessary checks
#' and processing to format the return object as required. Evaluation of this
#' function happens in the parent function.
#'
.mergeAndReturn <- function() {
  toDo <- expression(
    ## make sure res object has the required fields ##
    stopifnot("qID" %in% names(res)),
    
    ## make sure we only have one resulting row per query unless
    ## allSubjectCols is TRUE ##
    if (exists("allSubjectCols")) {
      if (!allSubjectCols) {
        stopifnot(!any(table(res$qID) > 1))
      }
    } else {
      stopifnot(!any(table(res$qID) > 1))
    },
    
    res <- DataFrame(res),
    
    ## setup new columns to be added using NA and add the proper class ##
    newCols <- grep(colnam,names(res),value = TRUE),
    mcols(sites.rd)[newCols] <- NA,
    for (newCol in newCols) {
      mcols(sites.rd)[[newCol]] <- as(mcols(sites.rd)[[newCol]],
                                      class(res[[newCol]]))
    },
    
    ## merge back results in the same order as rows in query/sites.rd ##
    rows <- match(mcols(sites.rd)$tempyID[good.rows], res$qID),
    mcols(sites.rd)[good.rows,][!is.na(rows),newCols] <-
      res[rows[!is.na(rows)], newCols],
    
    ## clear up any temp columns ##
    mcols(sites.rd)$tempyID <- NULL
  )
  
  eval.parent(toDo)
}

#-------------------- clustering ------------------------

standardizeSites <- function(unstandardizedSites){
  if( ! length(unstandardizedSites) > 0){
  return(unstandardizedSites)
  }
  #Get called start values for clustering  
  unstandardizedSites$Position <- ifelse(strand(unstandardizedSites) == "+", start(unstandardizedSites), end(unstandardizedSites))
  unstandardizedSites$Break <- ifelse(strand(unstandardizedSites) == "+", end(unstandardizedSites), start(unstandardizedSites))
  unstandardizedSites$Score <- 95
  unstandardizedSites$qEnd <- width(unstandardizedSites)
  
  #Positions clustered by 5L window and best position is chosen for cluster
  standardized <- clusterSites(
    psl.rd = unstandardizedSites,
    weight = rep(1, length(unstandardizedSites)) 
    )

  start(standardized) <- ifelse(strand(standardized) == "+", 
                                standardized$clusteredPosition, standardized$Break)
  end(standardized) <- ifelse(strand(standardized) == "-", 
                              standardized$clusteredPosition, standardized$Break)
  
  standardized$Position <- NULL
  standardized$Break <- NULL
  standardized$Score <- NULL
  standardized$qEnd <- NULL
  standardized$clusteredPosition <- NULL
  standardized$clonecount <- NULL
  standardized$clusterTopHit <- NULL
  
  sort(standardized)
}  

dereplicateSites <- function(sites){
  #Reduce sites which have the same starts, but loose range info
  #(no need to add a gapwidth as sites are standardized)
  sites.reduced <- flank(sites, -1, start=TRUE)
  sites.reduced <- unlist(reduce(sites.reduced, min.gapwidth = 0L, with.revmap=TRUE))
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  
  #Order original sites by revmap  
  dereplicatedSites <- sites[unlist(sites.reduced$revmap)]
  
  #Skip this step and provide similar output if length(sites) = 0
  if(length(sites) > 0){
    dereplicatedSites <- split(dereplicatedSites, Rle(values = seq(length(sites.reduced)), lengths = sites.reduced$counts))
  }  
  
  #Dereplicate reads with same standardized starts and provide the longeset width
  dereplicatedSites <- unlist(reduce(dereplicatedSites, min.gapwidth = 0L))
  mcols(dereplicatedSites) <- mcols(sites.reduced)

  dereplicatedSites
}



#------------------  HiReads -------------------------


#' Functions to process LM-PCR reads from 454/Illumina data.
#'
#' hiReadsProcessor contains set of functions which allow users to process 
#' LM-PCR products sequenced using any platform. Given an excel/txt file 
#' containing parameters for demultiplexing and sample metadata, the functions 
#' automate trimming of adaptors and identification of the genomic product.
#' Genomic products are further processed for QC and abundance quantification.
#'
#' @import Biostrings GenomicAlignments hiAnnotator BiocParallel xlsx plyr
#' sonicLength BiocGenerics
#' @docType package
#' @name hiReadsProcessor
#' @author Nirav V Malani
NULL

#' Sample Integration Sites Sequencing Data
#' 
#' This is a processed data object containing raw sequences and respective 
#' alignments to UCSC freeze hg18 from 112 integration site samples. The object
#' is of SimpleList class and follows a certain structural hierarchy explained
#' by the Introductory vignette.
#' 
#' @docType data
#' @keywords datasets
#' @format A SimpleList object 
#' @name seqProps
NULL

#' PSL file output
#' 
#' Sample BLAT PSL file output from samples included Integration Sites 
#' Sequencing Data \code{\link{seqProps}}
#' 
#' @docType data
#' @keywords datasets
#' @format a data frame of 1000 rows and 21 columns
#' @name psl
NULL

#' Read contents of a sequencing folder and make a SimpleList object
#'
#' Given a sequencing folder path, sample information file path, and sequence 
#' file extension pattern, the function returns a list of variables required to 
#' process the data. The function also calls \code{\link{read.sampleInfo}} 
#' which reads in sample processing metadata and formats it if needed.
#'
#' @param sequencingFolderPath full or relative path to the sequencing folder
#' @param sampleInfoFilePath full or relative path to the sample information 
#' file, which holds samples to quadrant/lane associations along with other 
#' metadata required to trim sequences or process it. Default to NULL, where 
#' the function tries to find xls or tab deliminated txt file in the sequencing 
#' folder which sounds similar to 'sampleinfo' and present you with choices of 
#' file to select from.
#' @param seqfilePattern regex/string to describe sequence file endings. 
#' See examples. Default is NULL.
#' @param interactive whether to prompt each time the function encounters an 
#' issue or use the defaults. Default is TRUE.
#'
#' @return a SimpleList list which is used by other functions to process and 
#' decode the data.
#'
#' @note 
#' \itemize{
#'   \item One must make sure that each sequencing file has sector name/number 
#'   prefixed at the beginning, else \code{\link{findBarcodes}} will fail trying 
#'   to find the filename.
#'   \item For paired end Illumina runs, make sure the filenames include R1, R2, 
#'   and I1 somewhere in the name denoting pair1, pair2, and index/barcode 
#'   reads, respectively. 
#' }
#'
#' @seealso \code{\link{read.sampleInfo}}, \code{\link{findBarcodes}}, 
#' \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples
#' runData <- system.file("extdata/FLX_sample_run/", 
#' package = "hiReadsProcessor")
#' read.SeqFolder(runData, seqfilePattern=".+fna.gz$")
#' \dontrun{
#' read.SeqFolder(".", seqfilePattern="\\.TCA.454Reads.fna$")
#' read.SeqFolder(".", seqfilePattern=".+fastq$")
#' read.SeqFolder(".", seqfilePattern=".+sff$")
#' }
read.SeqFolder <- function(sequencingFolderPath=NULL, sampleInfoFilePath=NULL, 
                           seqfilePattern=NULL, interactive=TRUE) {
  if(is.null(sequencingFolderPath)) {
    stop("No Sequencing Folder Path provided.")
  }
  
  ## get the sequencingFolderPath right!
  sequencingFolderPath <- normalizePath(sequencingFolderPath, mustWork=TRUE)
  
  seqFilePaths <- list.files(path=sequencingFolderPath, recursive=TRUE, 
                             full.names=TRUE, pattern=seqfilePattern)
  if(length(seqFilePaths)==0) {
    stop("No files found in the folder ", sequencingFolderPath,
         "matching following pattern: ", seqfilePattern)
  }
  
  ## read the sample info file
  if(is.null(sampleInfoFilePath)) {
    possibleFiles <- list.files(path=sequencingFolderPath, 
                                recursive=TRUE, full.names=TRUE, 
                                pattern=".*sampleinfo.+", ignore.case=TRUE)
    if (length(possibleFiles)==0) {
      stop("No sample information file found in folder: ", sequencingFolderPath)
    } else {
      if(interactive & length(possibleFiles)>1) {
        message("Please choose a sample information file to read the meta data 
                from:\n",
                paste(1:length(possibleFiles), possibleFiles, 
                      sep=": ", collapse="\n"))
        choice <- scan(what=integer(0), n=1, quiet=TRUE, multi.line=FALSE)
      } else {
        choice <- 1            
      }
      message("Choosing ", possibleFiles[choice], 
              " as sample information file.")
    }
    sampleInfoFilePath <- possibleFiles[choice]
  }
  sampleInfo <- read.sampleInfo(sampleInfoFilePath, interactive=interactive)
  

  ## do a quick test of filenames if any samples are from paired end illumina
  if(any(sampleInfo$pairedend)) {    
    sectors <- unique(sampleInfo$sector[sampleInfo$pairedend])
    for(sector in sectors) {
      vals <- grep(gsub("I1|R1|R2","",sector), seqFilePaths, value=TRUE)
      if(length(vals)!=3) {
        stop("Sector ",sector," is missing one of the files: R1, R2, or I1.")
      }
    }
  }
  
  if(length(sampleInfo)!=length(unique(gsub("I1|R1|R2","",seqFilePaths)))) {
    warning("Number of sectors (", length(sampleInfo),
            ") in sample information file doesn't match # of sector files (", 
            length(unique(gsub(seqfilePattern,'',seqFilePaths))), 
            ") found in the folder.")
  }
  
  return(SimpleList("sequencingFolderPath"=sequencingFolderPath, 
                    "seqFilePaths"=seqFilePaths, 
                    "seqfilePattern"=seqfilePattern, 
                    "sampleInfoFilePath"=sampleInfoFilePath, 
                    "sectors"=sampleInfo, "callHistory"=match.call()))
}

#' Read a sample information file and format appropriate metadata.
#'
#' Given a sample information file, the function checks if it includes required 
#' information to process samples present on each sector/quadrant/region/lane. 
#' The function also adds other columns required for processing with default 
#' values if not already defined ahead of time.
#'
#' @param sampleInfoPath full or relative path to the sample information file, 
#' which holds samples to quadrant/lane associations along with other metadata 
#' required to trim sequences or process it. 
#' @param splitBySector split the data frame into a list by sector column. 
#' Default is TRUE.
#' @param interactive whether to prompt each time the function encounters an 
#' issue, or use the defaults. Default is TRUE.
#'
#' @details
#' \itemize{
#'  \item Required Column Description:
#'    \itemize{
#'  	\item sector => region/quadrant/lane of the sequencing plate the sample 
#'    comes from. If files have been split by samples apriori, then the filename
#'    associated per sample without the extension. If this is a filename, then 
#'    be sure to enable 'alreadyDecoded' parameter in \code{\link{findBarcodes}}, 
#'    since contents of this column is pasted together with 'seqfilePattern' 
#'    parameter in \code{\link{read.SeqFolder}} to find the appropriate file 
#'    needed. For paired end data, this is basename of the FASTA/Q file holding 
#'    the sample data from the LTR side. For example, files such as 
#'    Lib3_L001_R2_001.fastq.gz or Lib3_L001_R2_001.fastq would be 
#'    Lib3_L001_R2_001, and consequently Lib3_L001_R1_001 would be used as the 
#'    second pair!
#'  	\item barcode => unique 4-12bp DNA sequence which identifies the sample. 
#'    If providing filename as sector, then leave this blank since it is assumed 
#'    that the data is already demultiplexed.
#'  	\item primerltrsequence => DNA sequence of the viral LTR primer 
#'    with/without the viral LTR sequence following the primer landing site. 
#'    If already trimmed, then mark this as SKIP.
#'  	\item sampleName => Name of the sample associated with the barcode
#'  	\item sampleDescription => Detailed description of the sample
#'  	\item gender => sex of the sample: male or female or NA
#'  	\item species => species of the sample: homo sapien, mus musculus, etc.
#'  	\item freeze => UCSC freeze to which the sample should be aligned to.
#'  	\item linkerSequence => DNA sequence of the linker adaptor following the 
#'    genomic sequence. If already trimmed, then mark this as SKIP.
#'  	\item restrictionEnzyme => Restriction enzyme used for digestion and 
#'    sample recovery. Can also be one of: Fragmentase or Sonication!
#' 		}
#'  \item Metadata Parameter Column Description:
#'   \itemize{
#'  	\item ltrBitSequence => DNA sequence of the viral LTR following the 
#'    primer landing site. Default is last 7bps of the primerltrsequence.
#'  	\item ltrBitIdentity => percent of LTR bit sequence to match during the 
#'    alignment. Default is 1.
#'  	\item primerLTRidentity => percent of primer to match during the 
#'    alignment. Default is .85
#'  	\item linkerIdentity => percent of linker sequence to match during the 
#'    alignment. Default is 0.55. Only applies to non-primerID/random tag based 
#'    linker search.
#'  	\item primerIdInLinker => whether the linker adaptor used has 
#'    primerID/random tag in it? Default is FALSE.
#'  	\item primerIdInLinkerIdentity1 => percent of sequence to match before 
#'    the random tag. Default is 0.75. Only applies to primerID/random tag 
#'    based linker search and when primeridinlinker is TRUE.
#'  	\item primerIdInLinkerIdentity2 => percent of sequence to match after the 
#'    random tag. Default is 0.50. Only applies to primerID/random tag based 
#'    linker search and when primeridinlinker is TRUE.
#'  	\item celltype => celltype information associated with the sample
#'  	\item user => name of the user who prepared or processed the sample
#'  	\item pairedEnd => is the data paired end? Default is FALSE.
#'    \item vectorFile => fasta file containing the vector sequence
#' 		}
#'  \item Processing Parameter Column Description:
#'   \itemize{
#'  	\item startWithin => upper bound limit of where the alignment should 
#'    start within the query. Default is 3.
#'  	\item alignRatioThreshold => cuttoff for (alignment span/read length). 
#'    Default is 0.7.
#'  	\item genomicPercentIdentity => cuttoff for (1-(misMatches/matches)). 
#'    Default is 0.98.
#'  	\item clusterSitesWithin => cluster integration sites within a defined 
#'    window size based on frequency which corrects for any sequencing errors. 
#'    Default is 5.
#'  	\item keepMultiHits => whether to keep sequences/reads that return 
#'    multiple best hits, aka ambiguous locations. 
#'  	\item processingDate => the date of processing 
#' 		}
#' }
#'
#' @return if splitBySector=TRUE, then an object of SimpleList named by 
#' quadrant/lane information defined in sampleInfo file, else a dataframe.
#'
#' @seealso \code{\link{read.SeqFolder}}, \code{\link{findBarcodes}}, 
#' \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' runData <- system.file("extdata/FLX_sample_run", 
#' package = "hiReadsProcessor")
#' read.sampleInfo(file.path(runData,"sampleInfo.xls"))
read.sampleInfo <- function(sampleInfoPath=NULL, splitBySector=TRUE, 
                            interactive=TRUE) {
  ## read file and make sampleInfo object with sample to metadata associations
  if(is.null(sampleInfoPath)) {
    stop("No sample information file path provided.")
  }
  
  sampleInfoPath <- normalizePath(sampleInfoPath, mustWork=TRUE)
  
  requiredCols <- c('sector', 'barcode', 'primerltrsequence', 'samplename', 
                    'sampledescription', 'gender', 'species', 'freeze', 
                    'linkersequence', 'restrictionenzyme')
  
  metaDataCols <- c('ltrbitsequence'='', 'ltrbitidentity'=1, 
                    'primerltridentity'=.85, 'linkeridentity'=.55, 
                    'primeridinlinker'=FALSE, 'primeridinlinkeridentity1'=.75, 
                    'primeridinlinkeridentity2'=.50, 'celltype'='',
                    'user'=Sys.getenv("USER"), 'startwithin'=3, 
                    'alignratiothreshold'=.7, 'clustersiteswithin'=5, 
                    'keepmultihits'=TRUE , 'genomicpercentidentity'=.98, 
                    'processingdate'=format(Sys.time(), "%Y-%m-%d "), 
                    'pairedend'=FALSE, 'vectorFile'='')
  
  if(grepl('.xls.?$', sampleInfoPath)) {
    sampleInfo <- unique(read.xlsx(sampleInfoPath, 
                                   sheetIndex=1, stringsAsFactors=FALSE))
  } else {
    sampleInfo <- unique(read.delim(sampleInfoPath, stringsAsFactors=FALSE))
  }
  names(sampleInfo) <- tolower(gsub("\\.|-|_", "", names(sampleInfo)))
  
  # check for required columns
  ColsNotThere <- !requiredCols %in% names(sampleInfo)
  if (any(ColsNotThere)) {
    absentCols <- requiredCols[ColsNotThere]
    stop("Following required column(s) is absent from the Sample Info file: ",
         paste(absentCols,sep="", collapse=", "))
  }
  
  # add missing meta data columns
  metaColsNotThere <- !names(metaDataCols) %in% names(sampleInfo)
  if(any(metaColsNotThere)) {
    sampleInfo <- cbind(sampleInfo,
                        as.data.frame(t(metaDataCols[metaColsNotThere]),
                                      stringsAsFactors = FALSE))        
  }
  
  # do some formatting to avoid later hassels!
  for(column in c('sector', 'barcode', 'primerltrsequence', 'ltrbitsequence', 
                  'samplename', 'linkersequence', 'restrictionenzyme')) {
    sampleInfo[,column] <- gsub(" ", "", sampleInfo[,column])
    if(column %in% c('barcode', 'primerltrsequence', 'ltrbitsequence', 
                     'linkersequence', 'restrictionenzyme')) {
      sampleInfo[,column] <- toupper(sampleInfo[,column])
    }
  }
  
  for(column in c('pairedend', 'keepmultihits', 'primeridinlinker')) {
    sampleInfo[,column] <- as.logical(sampleInfo[,column])
  }
  
  # confirm ltr bit is correct
  ltrbitTest <- sampleInfo$primerltrsequence=="SKIP"
  if(any(ltrbitTest)) { 
    ## add SKIP to ltrbit as well if primerltrsequence has been trimmed
    tofix <- which(ltrbitTest)
    message("adding SKIP to ltrbitsequence to ",length(tofix),
            " sample since primerltrsequence has been trimmed.")
    sampleInfo$ltrbitsequence[tofix] <- "SKIP"
  }
  
  ltrbitTest <- nchar(sampleInfo$ltrbitsequence)==0 | sampleInfo$ltrbitsequence==""
  if(any(ltrbitTest)) {
    tofix <- which(ltrbitTest)
    if(interactive) {
      message("LTR bit not found for ",length(tofix)," samples. 
              Use last 7 bases of the LTR primer as the LTR bit? (y/n)")
      choice <- scan(what=character(0), n=1, quiet=TRUE, multi.line=FALSE)
    } else {
      message("LTR bit not found for ",length(tofix),
              " samples. Using last 7 bases of the LTR primer as the LTR bit.")
      choice <- "y"
    }
    
    if(tolower(choice)=="y") {
      sampleInfo$ltrbitsequence <- substr(sampleInfo$primerltrsequence,
                                          nchar(sampleInfo$primerltrsequence)-6, 
                                          nchar(sampleInfo$primerltrsequence))
      sampleInfo$primerltrsequence <- substr(sampleInfo$primerltrsequence, 1,
                                             nchar(sampleInfo$primerltrsequence)-7)
    } else {
      warning("No LTR bit sequence found for following samples: ",
              paste(sampleInfo$samplename[tofix], sep="", collapse=", "),
              immediate.=TRUE)
    }        
  }
  
  # check if samplenames are up to the expectations
  samplenametest <- nchar(sampleInfo$samplename)==0 | sampleInfo$samplename==""
  if(any(samplenametest)) {
    stop("No sample names found in following rows of the sample information file ",
         sampleInfoPath, " : ", paste(which(samplenametest), 
                                      sep="", collapse=", "))
  }
  
  # check for sectors and their usage
  sectortest <- nchar(sampleInfo$sector)==0 | sampleInfo$sector=="" | 
    is.na(sampleInfo$sector)
  if(any(sectortest)) {
    tofix <- which(sectortest)
    
    if(interactive) {
      message("Sector information not found for ", length(tofix),
              " samples. Which sector are they from? (1,2,4,etc)")
      choice <- scan(what=character(0), quiet=TRUE, multi.line=FALSE)
    } else {
      message("Sector information not found for ", length(tofix),
              " samples. Assuming they are from sector 1.")
      choice <- "1"
    }
    
    if(length(choice)>0) {
      sampleInfo$sector[tofix] <- unlist(strsplit(choice,","))
    } else {
      stop("No Sector information found for following samples: ",
           paste(sampleInfo$samplename[tofix], sep="", collapse=", "))
    }
  }
  
  ## excel sometimes converts integers to doubles...
  ## make sure to remove the trailing 0
  sampleInfo$sector <- gsub("\\.0$", "", as.character(sampleInfo$sector))
  
  sampleSectorTest <- table(paste(sampleInfo$samplename, sampleInfo$sector))
  if(any(sampleSectorTest>1)) {
    stop("Duplicate sample names found on same quadrant in the ",
         "sample information file ", sampleInfoPath, " : ", 
         paste(sampleSectorTest[sampleSectorTest>1], sep="",collapse=", "))
  }
  
  # prepare the sample info object!
  if(splitBySector) {            
    sampleInfo <- SimpleList(split(sampleInfo, sampleInfo$sector))
    for(sector in 1:length(sampleInfo)) { 
      sampleData <- SimpleList(split(sampleInfo[[sector]], 
                                     sampleInfo[[sector]]$samplename))
      for(sample.i in 1:length(sampleData)) { 
        sampleData[[sample.i]] <- SimpleList(as.list(sampleData[[sample.i]])) 
      }
      sampleInfo[[sector]] <- SimpleList("samples"=sampleData) 
    }
  }
  return(sampleInfo)
}

#' Removes duplicate sequences from DNAStringSet object.
#'
#' Given a DNAStringSet object, the function dereplicates reads and 
#' adds counts=X to the definition line to indicate replication. 
#'
#' @param dnaSet DNAStringSet object to dereplicate. 
#'
#' @return DNAStringSet object with names describing frequency of repeat.
#'
#' @seealso \code{\link{replicateReads}}, \code{\link{removeReadsWithNs}}, 
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples 
#' dnaSet <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC", 
#' "CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC") 
#' dereplicateReads(dnaSet)
dereplicateReads <- function(dnaSet) {
  if(!is(dnaSet,"DNAStringSet")) {
    dnaSet <- DNAStringSet(dnaSet)
  }
  if(is.null(names(dnaSet))) {
    message("No names attribute found in dnaSet object...", 
            "using artifically generated names")
    names(dnaSet) <- paste("read", 1:length(dnaSet), sep="-")
  }
  dnaSet <- dnaSet[order(dnaSet)]
  counts <- BiocGenerics::table(dnaSet)
  dnaSet <- unique(dnaSet)
  names(dnaSet) <- paste0(names(dnaSet), 
                          "counts=", 
                          as.integer(counts[names(counts)[names(dnaSet)]]))
  return(dnaSet)
}

#' Replicate sequences from DNAStringSet object using counts identifier or vector
#'
#' Given a DNAStringSet object, the function replicates reads using counts=X 
#' marker at the end of definition line. 
#'
#' @param dnaSet DNAStringSet object to replicate. 
#' @param counts an integer or a numeric vector of length length(dnaSet) 
#' indicating how many times to repeat each sequence. Default is NULL, 
#' in which it uses counts=X notation from the definition line to 
#' replicate reads.
#'
#' @return DNAStringSet object.
#'
#' @seealso \code{\link{dereplicateReads}}, \code{\link{removeReadsWithNs}},
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' dnaSet <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC", 
#' "CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC") 
#' dnaSet <- dereplicateReads(dnaSet)
#' replicateReads(dnaSet)
replicateReads <- function(dnaSet, counts=NULL) {
  stopifnot(is(dnaSet,"DNAStringSet"))
  if(is.null(counts)) {
    if(is.null(names(dnaSet))) {
      stop("No names attribute found in dnaSet object")
    }
    counts <- as.numeric(sub(".+counts=(\\d+)","\\1", names(dnaSet)))
    if(all(is.na(counts))) {
      stop("No counts=X marker found at the end of definition line or ",
           "names attribute in dnaSet object")
    }
  }
  if (length(counts)==1) {
    counts <- rep(counts, length(dnaSet))
  }
  
  ids <- unlist(sapply(counts, seq_len))
  deflines <- paste0(rep(sub("(.+)counts=.+", "\\1", 
                             names(dnaSet)), times=counts), "_", ids)
  dnaSet <- rep(dnaSet, times=counts)
  names(dnaSet) <- deflines
  return(dnaSet)
}

#' Remove sequences with ambiguous nucleotides.
#'
#' Given a DNAStringSet object, the function removes any reads that has either repeating or total Ns which is greater than to maxNs threshold
#'
#' @param dnaSet DNAStringSet object to evaluate. 
#' @param maxNs integer value denoting the threshold of maximum allowed Ns. 
#' Default is 5.
#' @param consecutive boolean flag denoting whether Ns to filter is consecutive or total . Default is TRUE.
#'
#' @return DNAStringSet object.
#'
#' @seealso \code{\link{dereplicateReads}}, \code{\link{replicateReads}},
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples  
#' dnaSet <- c("CCTGAATCCTNNCAATGTCATCATC", "ATCCTGGCNATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGNNTCTGCAATGTGNGGNCCTAN", "GAAGNNNNNNGTTGAAGTTCACAC") 
#' removeReadsWithNs(dnaSet)
#' removeReadsWithNs(dnaSet, maxNs=4, consecutive=FALSE)
removeReadsWithNs <- function(dnaSet, maxNs=5, consecutive=TRUE) {
  if(!is(dnaSet,"DNAStringSet")) {
    dnaSet <- DNAStringSet(dnaSet)
  }
  if(consecutive) {
    good.row <- !grepl(paste(rep("N",maxNs+1), collapse=""), dnaSet, fixed=TRUE)
  } else {
    res <- alphabetFrequency(dnaSet)
    good.row <- res[,"N"] <= maxNs
  }
  return(dnaSet[good.row])
}

#' Breaks an object into chunks of N size.
#'
#' Given a linear/vector like object, the function breaks it into equal sized chunks either by chunkSize. This is a helper function used by functions in 'See Also' section where each chunk is sent to a parallel node for processing.
#'
#' @param x a linear object.
#' @param chunkSize number of rows to use per chunk of query. Defaults to length(x)/detectCores() or length(query)/bpworkers() depending on parallel backend registered. 
#'
#' @return a list of object split into chunks.
#'
#' @seealso \code{\link{primerIDAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{pairwiseAlignSeqs}}
#'
#' @export
#'
#' @examples
#' x <- c("GAGGCTGTCACCGACAAGGTTCTGA", "AATAGCGTGGTGACAGCCCACATGC", 
#'        "GGTCTTCTAGGGAACCTACGCCACA", "TTTCCGGCGGCAGTCAGAGCCAAAG", 
#'        "TCCTGTCAACTCGTAGATCCAATCA", "ATGGTCACCTACACACAACGGCAAT", 
#'        "GTCAGGACACCTAATCACAAACGGA", "AGACGCAGGTTCAGGCTGGACTAGA", 
#'        "ATCGTTTCCGGAATTCGTGCACTGG", "CAATGCGGGCACACGCTCTTACAGT")
#' chunkize(DNAStringSet(x), 5)
chunkize <- function(x, chunkSize = NULL) {
  chunks <- breakInChunks(length(x), 
                          ifelse(!is.null(chunkSize),
                                 length(x)/chunkSize, 
                                 ifelse(!is.null(bpworkers()), 
                                        length(x)/bpworkers(), 
                                        length(x)/detectCores())))
  mapply(function(z, y) x[z:y], start(chunks), end(chunks), 
         SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

#' Split DNAStringSet object using first X number of bases defined by a vector.
#'
#' Given a character vector of barcodes/MID to sample association and a DNAStringSet object, the function splits/demultiplexes the DNAStringSet object by first few bases dictated by length of barcodes/MID supplied. This is an accessory function used by \code{\link{findBarcodes}}
#'
#' @param barcodesSample a character vector of barcodes to sample name associations. Ex: c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2",...)
#' @param dnaSet DNAStringSet object to evaluate. 
#' @param trimFrom integer value serving as start point to trim the sequences from. This is calculated internally length barcode+1. Default is NULL.
#' @param showStats boolean flag denoting whether to show decoding statistics per sample & barcode. Default is FALSE.
#' @param returnUnmatched boolean flag denoting whether to return unmatched reads. Default is FALSE.
#'
#' @return DNAStringSet object split by sample name found in barcodesSample.
#'
#' @seealso \code{\link{findBarcodes}}, \code{\link{dereplicateReads}},
#' \code{\link{replicateReads}}
#'
#' @export
#'
#' @examples 
#' dnaSet <- DNAStringSet(c("read1"="ACATCCATAGAGCTACGACGACATCGACATA",
#' "read2"="GAATGGATGACGACTACAGCACGACGAGCAGCTACT",
#' "read3"="GAATGGATGCGCTAAGAAGAGA", "read4"="ACATCCATTCTACACATCT"))
#' splitByBarcode(c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"), dnaSet, 
#' showStats=TRUE)
splitByBarcode <- function(barcodesSample, dnaSet, trimFrom=NULL, 
                           showStats=FALSE, returnUnmatched=FALSE) {
  if(is.null(barcodesSample) | length(barcodesSample)==0) {
    stop("No barcodes to samples association vector provided in parameter ",
         "barcodesSample.")
  }
  
  if(is.null(dnaSet) | length(dnaSet)==0) {
    stop("No sequences provided in parameter dnaSet.")
  }
  
  if(is.null(names(dnaSet))) {
    stop("No names attribute found in dnaSet object")
  }
  
  message("Using following schema for barcode to sample associations")
  print(as.data.frame(barcodesSample))
  
  ## subset barcode string from the rest of sequence ##
  barcodelen <- unique(nchar(names(barcodesSample)))
  seqbarcodes <- substr(dnaSet,1,barcodelen)
  
  ## index rows that match your list of barcodes ##
  good.row <- seqbarcodes %in% names(barcodesSample)
  if(!any(good.row)) {
    stop("No matching barcoded sequences found on this quadrant.")
  }
  
  sampleNames <- barcodesSample[seqbarcodes[good.row]]    
  deflines <- sub("^(\\S+) .+$", "\\1", names(dnaSet)[good.row], perl=TRUE)
  
  ## if primer bases were utilized for tiebreaker, use the original length 
  ## instead of modified for trimming.
  if(is.null(trimFrom)) {        
    trimFrom <- barcodelen+1
  }
  
  ## remove sequences with unknown barcode and trim barcode itself ##
  unmatched <- DNAStringSet(dnaSet[!good.row])
  dnaSet <- DNAStringSet(dnaSet[good.row], start=trimFrom)
  names(dnaSet) <- deflines
  
  if(showStats) {
    message("Number of Sequences with no matching barcode: ",
            as.numeric(table(good.row)['FALSE']))
    message("Number of Sequences decoded:")
    print(as.data.frame(table(sampleNames)))
  }
  
  dnaSet <- as.list(split(dnaSet, as.character(sampleNames)))
  
  if(returnUnmatched) {
    dnaSet <- c(dnaSet, "unDecodedSeqs"=unmatched)
  }
  
  return(dnaSet)
}

#' Demultiplex reads by their barcodes
#'
#' Given a sample information object, the function reads in the raw fasta/fastq file, demultiplexes reads by their barcodes, and appends it back to the sampleInfo object. Calls \code{\link{splitByBarcode}} to perform the actual splitting of file by barcode sequences. If supplied with a character vector and reads themselves, the function behaves a bit differently. See the examples.
#'
#' @param sampleInfo sample information SimpleList object created using 
#' \code{\link{read.sampleInfo}}, which holds barcodes and sample names per 
#' sector/quadrant/lane or a character vector of barcodes to sample name 
#' associations. Ex: c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2",...)
#' @param sector If sampleInfo is a SimpleList object, then a numeric/character 
#' value or vector representing sector(s) in sampleInfo. Optionally if on high 
#' memory machine sector='all' will decode/demultiplex sequences from all 
#' sectors/quadrants. This option is ignored if sampleInfo is a character vector. 
#' Default is NULL. 
#' @param dnaSet DNAStringSet object containing sequences to be decoded or 
#' demultiplexed. Default is NULL. If sampleInfo is a SimpleList object, then 
#' reads are automatically extracted using \code{\link{read.seqsFromSector}} 
#' and parameters defined in sampleInfo object.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param returnUnmatched return unmatched sequences. Returns results as a list 
#' where x[["unDecodedSeqs"]] has culprits. Default is FALSE.
#' @param dereplicate return dereplicated sequences. Calls 
#' \code{\link{dereplicateReads}}, which appends counts=X to sequence 
#' names/deflines. Default is FALSE. Not applicable for paired end data since 
#' it can cause insyncronicity.
#' @param alreadyDecoded if reads have be already decoded and split into 
#' respective files per sample and 'seqfilePattern' parameter in 
#' \code{\link{read.SeqFolder}} is set to reading sample files and not the 
#' sector files, then set this to TRUE. Default is FALSE. Enabling this 
#' parameter skips the barcode detection step and loads the sequence file as is 
#' into the sampleInfo object. 
#'
#' @return If sampleInfo is an object of SimpleList then decoded sequences are 
#' appeneded to respective sample slots, else a named list of DNAStringSet 
#' object. If returnUnmatched=TRUE, then x[["unDecodedSeqs"]] has the 
#' unmatched sequences.
#'
#' @seealso \code{\link{splitByBarcode}}, \code{\link{dereplicateReads}},
#' \code{\link{replicateReads}}
#'
#' @export
#'
#' @aliases decodeByBarcode
#' 
#' @examples 
#' dnaSet <- DNAStringSet(c("read1"="ACATCCATAGAGCTACGACGACATCGACATA",
#' "read2"="GAATGGATGACGACTACAGCACGACGAGCAGCTACT",
#' "read3"="GAATGGATGCGCTAAGAAGAGA", "read4"="ACATCCATTCTACACATCT"))
#' findBarcodes(sampleInfo=c("ACATCCAT"="Sample1", "GAATGGAT"="Sample2"), 
#' dnaSet=dnaSet, showStats=TRUE, returnUnmatched=TRUE)
#' \dontrun{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' findBarcodes(seqProps, sector="all", showStats=TRUE)
#' }
findBarcodes <- function(sampleInfo, sector=NULL, dnaSet=NULL, 
                         showStats=FALSE, returnUnmatched=FALSE, 
                         dereplicate=FALSE, alreadyDecoded=FALSE) {
  
  ## tried PDict()...and its slower than this awesome code! ##    
  if(is(sampleInfo,"SimpleList")) {
    if(is.null(sector)) {
      stop("No sector provided in parameter sector.")
    }
    
    sectors <- sector <- as.character(sector)
    if(length(sectors)==1 & tolower(sectors)=="all") {
      sectors <- names(sampleInfo$sectors)
    }
    
    if(any(!sectors %in% names(sampleInfo$sectors))) {
      stop("Following sectors not found in names(sampleInfo$sectors): ",
           sectors[!sectors %in% names(sampleInfo$sectors)])
    }
    
    for(sector in sectors) {
      ## check everything is cool with the provided barcodes first before 
      ## reading the sequences!
      message("Decoding sector: ",sector) 
      isPaired <- any(as.logical(extractFeature(sampleInfo, sector, 
                                                feature='pairedend')[[sector]]))
      
      ## prepare a vector of barcode to sample associations ##
      sampleBarcodes <- toupper(extractFeature(sampleInfo, sector=sector,
                                               feature="barcode")[[sector]])
      barcodesSample <- structure(names(sampleBarcodes), 
                                  names=as.character(sampleBarcodes))
      
      if (length(table(nchar(as.character(sampleBarcodes))))>1) {
        stop("Multiple barcode lengths found.")
      }
      
      ## length of barcodes before any modifications done later if any! ##
      realbarcodelen <- unique(nchar(as.character(sampleBarcodes)))
      
      if (any(table(as.character(sampleBarcodes))>1)) {
        message("Duplicate barcode found on this sector.\n",
                "Please choose from one of the options below:\n",
                "\t1: Pick first few bases of primer for tiebreaker? ",
                "(This could be dangerous if the run has too many errors!)\n",
                "\t2: Use the last sample associated with the duplicate ",
                "as the primary sample?\n", 
                "\t3: Do not do anything.")
        
        choice <- scan(what=integer(0), n=1, quiet=TRUE, multi.line=FALSE)
        
        if(choice==1) {  
          message("Enter # of bases to use from primer:")
          howmany <- scan(what=integer(0), n=1, quiet=TRUE, multi.line=FALSE)
          samplePrimers <- extractFeature(sampleInfo,
                                          sector=sector,
                                          feature="primerltrsequence")[[sector]]
          samplePrimers <- toupper(samplePrimers)
          newBarcodes <- toupper(paste0(sampleBarcodes,
                                        substr(samplePrimers,1,howmany))) 
          counts <- table(newBarcodes)
          ## only take 1 to 1 associations!
          rows <- newBarcodes %in% names(which(counts==1)) 
          
          if (any(counts>1)) {            
            message("Tie breaking failed...",
                    "try choosing high number of bases from primer possibly? ",
                    "Here are the failed barcodes: ",
                    paste(names(which(counts>1)),collapse=", "))
            
            message("Ignore samples associated with those barcodes and ",
                    "continue processing? (y/n)")
            whatsup <- scan(what=character(0), n=1, quiet=TRUE, 
                            multi.line=FALSE)            
            if(whatsup=='n') {
              stop("Aborting processing due to ambiguous barcode ",
                   "association for samples: ",
                   paste(names(sampleBarcodes[!rows]),collapse=", "))
            } else {              
              message("Ignoring following samples due to duplicate barcodes: ",
                      paste(names(sampleBarcodes[!rows]),collapse=", "))
            }
          }
          
          barcodesSample <- structure(names(sampleBarcodes[rows]), 
                                      names=newBarcodes[rows])
          
        } else if(choice==2) {
          message("Overwriting duplicate samples associated with the same barcode...")
        } else {
          stop("Aborting due to duplicate barcode found on this sector")
        }
      }
      
      dnaSet <- read.seqsFromSector(sampleInfo, sector, isPaired)
      
      if(alreadyDecoded) {
        if(length(barcodesSample)>1) {
          stop("alreadyDecoded parameter is set to TRUE. There shouldn't be more ",
               "than one sample associated to a sequence file.")
        }
        
        if(is.list(dnaSet)) {
          dnaSet <- sapply(dnaSet, function(x) {
            names(x) <- sub("^\\S+-(\\S+) .+$", "\\1", names(x), perl=TRUE)  
            x
          })          
        } else {
          names(dnaSet) <- sub("^\\S+-(\\S+) .+$", "\\1", 
                               names(dnaSet), perl=TRUE) 
        }        
        
        if(isPaired) {
          ## no need to store barcode/index reads if alreadyDecoded!
          dnaSet <- list(dnaSet[c("pair1","pair2")])
        } else {          
          if(dereplicate) {
            dnaSet <- list(dereplicateReads(dnaSet))
          } else {
            dnaSet <- list(dnaSet)
          }
        }
        names(dnaSet) <- as.character(barcodesSample)
      } else {
        if(isPaired) {
          bc <- splitByBarcode(barcodesSample, dnaSet[["barcode"]], 
                               trimFrom=realbarcodelen+1, showStats=showStats, 
                               returnUnmatched=returnUnmatched)
          p1 <- sapply(bc, function(x) dnaSet[['pair1']][names(x)])
          p2 <- sapply(bc, function(x) dnaSet[['pair2']][names(x)])
          stopifnot(identical(sapply(bc,length), sapply(p1,length)))
          stopifnot(identical(sapply(bc,length), sapply(p2,length)))
          dnaSet <- mapply(function(x,y) list("pair1"=x, "pair2"=y), p1, p2, 
                           SIMPLIFY=FALSE)
          rm("bc","p1","p2")
        } else {
          dnaSet <- splitByBarcode(barcodesSample, dnaSet, 
                                   trimFrom=realbarcodelen+1, 
                                   showStats=showStats, 
                                   returnUnmatched=returnUnmatched)
          if(dereplicate) {
            dnaSet <- dereplicateReads(dnaSet)
          }
        }
      }
      
      for(samplename in names(dnaSet)) {
        if(samplename=="unDecodedSeqs") {                                        
          metadata(sampleInfo$sectors[[sector]]) <- 
            append(metadata(sampleInfo$sectors[[sector]]), 
                   list("unDecodedSeqs"=dnaSet[[samplename]]))            
        } else {
          sampleInfo$sectors[[sector]]$samples[[samplename]]$decoded <- 
            dnaSet[[samplename]]
        }
      }
      metadata(sampleInfo$sectors[[sector]]) <- 
        append(metadata(sampleInfo$sectors[[sector]]),
               list("decodedBy"=barcodesSample))
    }
    sampleInfo$callHistory <- append(sampleInfo$callHistory, match.call())
    decoded <- sampleInfo
    cleanit <- gc()
  } else {
    decoded <- splitByBarcode(sampleInfo, dnaSet, trimFrom=NULL, 
                              showStats=showStats, 
                              returnUnmatched=returnUnmatched)
    cleanit <- gc()
  }
  return(decoded)
}

#' @export
decodeByBarcode <- findBarcodes

#' Align a short pattern to variable length target sequences.
#'
#' Align a fixed length short pattern sequence (i.e. primers or adaptors) to subject sequences using \code{\link{pairwiseAlignment}}. This function uses default of type="overlap", gapOpening=-1, and gapExtension=-1 to align the patternSeq against subjectSeqs. One can adjust these parameters if prefered, but not recommended. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT or see one of following \code{\link{blatSeqs}}, \code{\link{findAndRemoveVector}}
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. This is generally bigger than patternSeq, and cases where subjectSeqs is smaller than patternSeq will be ignored in the alignment.
#' @param patternSeq DNAString object or a sequence containing the query sequence to search. This is generally smaller than subjectSeqs. 
#' @param side which side of the sequence to perform the search: left, right or middle. Default is 'left'.
#' @param qualityThreshold percent of patternSeq to match. Default is 1, full match.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param bufferBases use x number of bases in addition to patternSeq length to perform the search. Beneficial in cases where the pattern has homopolymers or indels compared to the subject. Default is 5. Doesn't apply when side='middle'.
#' @param doRC perform reverse complement search of the defined pattern. Default is TRUE.
#' @param returnUnmatched return sequences which had no or less than 5\% match 
#' to the patternSeq. Default is FALSE.
#' @param returnLowScored return sequences which had quality score less than 
#' the defined qualityThreshold. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to FALSE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#' @param ... extra parameters for \code{\link{pairwiseAlignment}}
#' @note 
#' \itemize{
#'  \item For qualityThreshold, the alignment score is calculated by 
#'  (matches*2)-(mismatches+gaps) which programatically translates to 
#'  round(nchar(patternSeq)*qualityThreshold)*2.
#'  \item Gaps and mismatches are weighed equally with value of -1 which can 
#'  be overriden by defining extra parameters 'gapOpening' & 'gapExtension'.
#'  \item If qualityThreshold is 1, then it is a full match, if 0, then any 
#'  match is accepted which is useful in searching linker sequences at 3' end. 
#'  Beware, this function only searches for the pattern sequence in one 
#'  orientation. If you are expecting to find the pattern in both orientation, 
#'  you might be better off using BLAST/BLAT!
#'  \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'  before running the function. One can use any of the following 
#'  \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#'
#' @return 
#' \itemize{
#'  \item IRanges object with starts, stops, and names of the aligned sequences. 
#'  \item If returnLowScored or returnUnmatched = T, then a CompressedIRangesList 
#'  where x[["hits"]] has the good scoring hits, x[["Rejected"]] has the failed 
#'  to match qualityThreshold hits, and x[["Absent"]] has the hits where the 
#'  aligned bit is <=10\% match to the patternSeq.
#' }
#'
#' @seealso \code{\link{primerIDAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, 
#' \code{\link{doRCtest}}, \code{\link{findAndTrimSeq}}, \code{\link{blatSeqs}},
#' \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' subjectSeqs <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC")
#' subjectSeqs <- DNAStringSet(xscat("AAAAAAAAAA", subjectSeqs))
#' pairwiseAlignSeqs(subjectSeqs, "AAAAAAAAAA", showStats=TRUE)
#' pairwiseAlignSeqs(subjectSeqs, "AAATAATAAA", showStats=TRUE, 
#' qualityThreshold=0.5)
#'
pairwiseAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, side="left", 
                              qualityThreshold=1, showStats=FALSE, bufferBases=5, 
                              doRC=TRUE, returnUnmatched=FALSE, 
                              returnLowScored=FALSE, parallel=FALSE, ...) {
  dp <- NULL
  
  .checkArgs_SEQed()
  
  if(parallel) {
    subjectSeqs2 <- chunkize(subjectSeqs)
    hits <- bplapply(subjectSeqs2, function(x) 
      pairwiseAlignSeqs(x, patternSeq, side, qualityThreshold, showStats=FALSE, 
                        bufferBases, doRC, returnUnmatched,  returnLowScored, 
                        parallel=FALSE, ...), BPPARAM=dp)    
    hits <- do.call(c, hits)
    if(is(hits,"CompressedIRangesList")) {
      attrs <- unique(names(hits))
      hits <- sapply(attrs, 
                     function(x) unlist(hits[names(hits)==x],use.names=FALSE))
      IRangesList(hits)
    } else {
      hits
    }
  } else {
    qualityThreshold <- as.numeric(qualityThreshold)
    
    ## only get the relevant side of subject sequence with extra bufferBases to 
    ## account for indels/mismatches & save memory while searching and avoid 
    ## searching elsewhere in the sequence
    if(tolower(side)=="left") {
      badSeqs <- DNAStringSet()
      culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
      if(any(culprits)) {
        badSeqs <- subjectSeqs[culprits]
        message(length(badSeqs),
                " sequences were removed from aligning since they were",
                " shorter than pattern getting aligned: ",
                (nchar(patternSeq)+bufferBases),"bp")            
        subjectSeqs <- subjectSeqs[!culprits]            
      }
      subjectSeqs2 <- subseq(subjectSeqs, start=1, 
                             end=(nchar(patternSeq)+bufferBases))
      overFromLeft <- rep(0,length(subjectSeqs))
    } else if (tolower(side)=="right") { 
      overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
      overFromLeft[overFromLeft<1] <- 1
      subjectSeqs2 <- subseq(subjectSeqs, start=overFromLeft)
    } else {
      subjectSeqs2 <- subjectSeqs
      overFromLeft <- rep(0, length(subjectSeqs))
    }
    
    ## search both ways to test which side yields more hits!
    if(doRC) {
      patternSeq <- tryCatch(doRCtest(subjectSeqs2, patternSeq, 
                                      qualityThreshold),
                             error=function(e) patternSeq)
    }
    
    ## type=overlap is best for primer trimming...see Biostrings Alignment vignette
    if(any(names(match.call()) %in% c("type","gapOpening","gapExtension"))) {
      hits <- pairwiseAlignment(subjectSeqs2, patternSeq, ...)        
    } else {
      hits <- pairwiseAlignment(subjectSeqs2, patternSeq, type="overlap", 
                                gapOpening=-1, gapExtension=-1, ...)
    }
    stopifnot(length(hits)==length(subjectSeqs2))
    
    scores <- round(score(hits))
    highscored <- scores >= round(nchar(patternSeq)*qualityThreshold)*2
    
    if(!any(highscored)) {
      stop("No hits found which passed the qualityThreshold")
    }
    
    ## basically a small subset of highscored
    unmatched <- nchar(hits) <= round(nchar(patternSeq)*.1) 
    
    # no point in showing stats if all sequences are a potential match! #
    if(showStats & qualityThreshold!=0) {
      bore <- as.numeric(table(highscored)['FALSE'])
      message("Total of ",ifelse(is.na(bore),0,bore),
              " did not have the defined pattern sequence (", patternSeq,
              ") that passed qualityThreshold on the ", side, " side")
    }
    
    ## extract starts-stops of the entire pattern hit ##
    starts <- start(pattern(hits))
    ends <- end(pattern(hits))
    namesq <- names(subjectSeqs)
    hits <- IRanges(start=starts+overFromLeft-ifelse(side=="right",2,0), 
                    end=ends+overFromLeft-ifelse(side=="right",2,0), 
                    names=namesq)
    rm("scores","subjectSeqs2","subjectSeqs","starts","ends","namesq")
    
    ## no need to test if there were any multiple hits since pairwiseAlignment 
    ## will only output one optimal alignment...see the man page.
    if(!returnLowScored & !returnUnmatched) {      
      hits <- hits[highscored]
    } else {
      hitstoreturn <- IRangesList("hits"=hits[highscored])
      if(returnLowScored & length(hits[!highscored])>0) {
        hitstoreturn <- append(hitstoreturn, 
                               IRangesList("Rejected"=hits[!highscored]))
      }
      
      if(returnUnmatched & length(hits[unmatched])>0) {
        hitstoreturn <- append(hitstoreturn, 
                               IRangesList("Absent"=hits[unmatched]))
      }
      hits <- hitstoreturn
      rm(hitstoreturn)
    }             
    cleanit <- gc()
    return(hits)
  }
}

#' Align a short pattern with PrimerID to variable length target sequences.
#'
#' Align a fixed length short pattern sequence containing primerID to variable length subject sequences using \code{\link{pairwiseAlignment}}. This function uses default of type="overlap", gapOpening=-1, and gapExtension=-1 to align the patterSeq against subjectSeqs. The search is broken up into as many pieces +1 as there are primerID and then compared against subjectSeqs. For example, patternSeq="AGCATCAGCANNNNNNNNNACGATCTACGCC" will launch two search jobs one per either side of Ns. For each search, qualityThreshold is used to filter out candidate alignments and the area in between is chosen to be the primerID. This strategy is benefical because of Indels introduced through homopolymer errors. Most likely the length of primerID(s) wont the same as you expected!
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. 
#' @param patternSeq DNAString object or a sequence containing the query sequence to search with the primerID.
#' @param qualityThreshold1 percent of first part of patternSeq to match. Default is 0.75.
#' @param qualityThreshold2 percent of second part of patternSeq to match. Default is 0.50.
#' @param doAnchored for primerID based patternSeq, use the base before and after primer ID in patternSeq as anchors?. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern. Default is TRUE.
#' @param returnUnmatched return sequences if it had no or less than 5\% match to the first part of patternSeq before the primerID. Default is FALSE.
#' @param returnRejected return sequences if it only has a match to one side of patternSeq or primerID length does not match # of Ns +/-2 in the pattern. Default is FALSE.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param ... extra parameters for \code{\link{pairwiseAlignment}}
#'
#' @note 
#' \itemize{
#'  \item For qualityThreshold1 & qualityThreshold2, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2
#'  \item Gaps and mismatches are weighed equally with value of -1 which can be overriden by defining extra parameters 'gapOpening' & 'gapExtension'.
#'  \item If qualityThreshold is 1, then it is a full match, if 0, then any match is accepted which is useful in searching linker sequences at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST/BLAT!
#' }
#'
#' @return 
#' \itemize{
#'  \item A CompressedIRangesList of length two, where x[["hits"]] is hits covering the entire patternSeq, and x[["primerIDs"]] is the potential primerID region. 
#'  \item If returnUnmatched = T, then x[["Absent"]] is appended which includes reads not matching the first part of patternSeq. 
#'  \item If returnRejected=TRUE, then x[["Rejected"]] includes reads that only matched one part of patternSeq or places where no primerID was found in between two part of patternSeq, and x[["RejectedprimerIDs"]] includes primerIDs that didn't match the correct length. 
#'  \item If doAnchored=TRUE, then x[["unAnchoredprimerIDs"]] includes reads that didn't match the base before and after primer ID on patternSeq.
#' }
#'
#' @seealso \code{\link{vpairwiseAlignSeqs}}, \code{\link{pairwiseAlignSeqs}},
#' \code{\link{doRCtest}}, \code{\link{blatSeqs}}, 
#' \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' subjectSeqs <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC")
#' ids <- c("GGTTCTACGT", "AGGAGTATGA", "TGTCGGTATA", "GTTATAAAAC", 
#' "AGGCTATATC", "ATGGTTTGTT")
#' subjectSeqs <- xscat(subjectSeqs, xscat("AAGCGGAGCCC",ids,"TTTTTTTTTTT"))
#' patternSeq <- "AAGCGGAGCCCNNNNNNNNNNTTTTTTTTTTT"
#' primerIDAlignSeqs(DNAStringSet(subjectSeqs), patternSeq, doAnchored = TRUE)
primerIDAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, 
                              qualityThreshold1=0.75, qualityThreshold2=0.50, 
                              doAnchored=FALSE, doRC=TRUE, returnUnmatched=FALSE, 
                              returnRejected=FALSE, showStats=FALSE, ...) {
  
  .checkArgs_SEQed()
  
  qualityThreshold1 <- as.numeric(qualityThreshold1)
  qualityThreshold2 <- as.numeric(qualityThreshold2)
  
  ## make sure there are Ns in the patternSeq for considering primerIDs
  if(length(unlist(gregexpr("N",patternSeq)))<4) {
    stop("There should be minimum of atleast 4 Ns in patternSeq to be ",
         "considered as a primerID sequence.")
  }
  
  ## Get the right orientation of the supplied patternSeq to peform proper search 
  ## at later two step search phase. 
  if(doRC) {
    patternSeq <- tryCatch(doRCtest(subjectSeqs, patternSeq),
                           error=function(e) patternSeq)
  }
  
  primerIDpos <- unlist(gregexpr("N", patternSeq))
  
  ## Perform primerID extraction by breaking the pattern into two parts...
  ## for sanity sakes due to homopolymers ##
  pattern1 <- as.character(subseq(DNAString(patternSeq), 1, primerIDpos[1]-1))
  pattern2 <- as.character(subseq(DNAString(patternSeq), 
                                  primerIDpos[length(primerIDpos)]+1))
  
  pattern1.hits <- pairwiseAlignSeqs(subjectSeqs, pattern1, "middle", 
                                     qualityThreshold=qualityThreshold1, 
                                     doRC=FALSE, returnUnmatched=TRUE, ...)
  pattern2.hits <- pairwiseAlignSeqs(subjectSeqs, pattern2, "middle", 
                                     qualityThreshold=qualityThreshold2,
                                     doRC=FALSE, ...)
  
  ## Set aside reads which has no match to the pattern1...
  ## most likely mispriming if primerID is on 5' or read was too loong if on 3'
  ## Hits returned from pairwiseAlignSeqs will be filtered for low scored hits...
  ## so no need to check for those from subjectSeqs
  unmatched <- pattern1.hits[["Absent"]]
  pattern1.hits <- pattern1.hits[["hits"]]
  
  ## remove reads which only have a match to one of the patterns...crossover most likely!
  rejected1 <- setdiff(names(pattern1.hits), names(pattern2.hits))
  rejected2 <- setdiff(names(pattern2.hits), names(pattern1.hits))
  rejected <- c(pattern1.hits[names(pattern1.hits) %in% rejected1], 
                pattern2.hits[names(pattern2.hits) %in% rejected2])
  rm("rejected1","rejected2")
  if(showStats) { 
    message("Removed ", length(rejected),
            " read(s) for only matching one of pattern1 or pattern2") 
  }
  
  ## use only reads which match to both sides of the patterns.
  good.rows <- intersect(names(pattern1.hits), names(pattern2.hits))
  pattern1.hits <- pattern1.hits[names(pattern1.hits) %in% good.rows]
  pattern2.hits <- pattern2.hits[names(pattern2.hits) %in% good.rows]
  
  stopifnot(identical(names(pattern1.hits), names(pattern2.hits)))
  stopifnot(identical(names(pattern1.hits), good.rows))
  
  ## make sure there is no overlap of ranges between pattern1.hits & 
  ## pattern2.hits if there is, then no primerID was found...remove it
  badAss <- end(pattern1.hits) >= start(pattern2.hits)
  if(any(badAss)) {
    rejected <- c(rejected, pattern1.hits[badAss], pattern2.hits[badAss])
    pattern1.hits <- pattern1.hits[!badAss]
    pattern2.hits <- pattern2.hits[!badAss]
    good.rows <- good.rows[!badAss]
    bore <- table(badAss)["TRUE"]
    message("Removed ", ifelse(is.na(bore),0,bore),
            " read(s) for not having primerID present between pattern 1 & 2")
  }
  
  hits <- IRanges(start=start(pattern1.hits), 
                  end=end(pattern2.hits), 
                  names=good.rows)        
  
  primerIDs <- IRanges(start=end(pattern1.hits)+1, 
                       end=start(pattern2.hits)-1, 
                       names=good.rows)        
  
  if(length(hits)==0) {
    stop("No hits found that matched both sides of patternSeq with ",
         "primerID in the middle.")
  }    
  
  ## do anchored search for only sequences that matched both sides of patternSeq
  if(doAnchored) {
    message("Found ", length(primerIDs), 
            " total primerIDs before anchored filter.")
    
    ## get anchors of bases flanking Ns
    anchorBase.s <- substr(patternSeq, primerIDpos[1]-1, primerIDpos[1]-1)
    anchorBase.e <- substr(patternSeq, primerIDpos[length(primerIDpos)]+1, 
                           primerIDpos[length(primerIDpos)]+1)
    
    anchorBase.s.i <- trimSeqs(subjectSeqs, 
                               resize(pattern1.hits,width=1,fix="end"))
    anchorBase.e.i <- trimSeqs(subjectSeqs, 
                               resize(pattern2.hits,width=1,fix="start"))
    rows <- anchorBase.s==as.character(anchorBase.s.i) & 
      anchorBase.e==as.character(anchorBase.e.i)
    
    unAnchored <- hits[!rows]
    primerIDs <- primerIDs[rows]
    hits <- hits[rows]
    message("Found ", length(primerIDs), 
            " total primerIDs after anchored filter.")
    rm("rows","anchorBase.s.i","anchorBase.e.i")
  }
  rm("good.rows","pattern1.hits","pattern2.hits")
  cleanit <- gc()
  
  ## also remove any primerIDs that were too short or big than intended.
  nSize <- 2
  badAss <- !width(primerIDs) %in% 
    (length(primerIDpos)-nSize):(length(primerIDpos)+nSize)
  if(any(badAss)) {
    rejectedprimerIDs <- hits[badAss]
    hits <- hits[!badAss]
    primerIDs <- primerIDs[!badAss]
    message("Removed ", table(badAss)["TRUE"],
            " read(s) for not having right primerID length")
  }
  
  hits <- IRangesList("hits"=hits, "primerIDs"=primerIDs)
  
  if(exists("unAnchored")) {
    if(length(unAnchored)>0) { 
      hits <- append(hits, IRangesList("unAnchoredprimerIDs"=unAnchored)) 
    }
  }
  
  if(returnUnmatched & length(unmatched)>0) {
    hits <- append(hits, IRangesList("Absent"=unmatched))
  }
  
  if(returnRejected) {
    if(length(rejected)>0) { 
      hits <- append(hits, IRangesList("Rejected"=rejected)) 
    }
    if(exists("rejectedprimerIDs")) { 
      if(length(rejectedprimerIDs)>0) { 
        hits <- append(hits, IRangesList("RejectedprimerIDs"=rejectedprimerIDs)) 
      } 
    }
  }    
  return(hits)
}

#' Align a short pattern to variable length target sequences.
#'
#' Align a fixed length short pattern sequence to subject sequences using \code{\link{vmatchPattern}}. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT.
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. This is generally bigger than patternSeq, and cases where subjectSeqs is smaller than patternSeq will be ignored in the alignment.
#' @param patternSeq DNAString object or a sequence containing the query sequence to search. This is generally smaller than subjectSeqs. 
#' @param side which side of the sequence to perform the search: left, right, or middle. Default is 'left'.
#' @param qualityThreshold percent of patternSeq to match. Default is 1, full match. This is supplied to max.mismatch parameter of \code{\link{vmatchPattern}} as round(nchar(patternSeq)*(1-qualityThreshold)).
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param bufferBases use x number of bases in addition to patternSeq length to perform the search. Beneficial in cases where the pattern has homopolymers or indels compared to the subject. Default is 5. Doesn't apply when side='middle'.
#' @param doRC perform reverse complement search of the defined pattern. Default is TRUE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to FALSE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#' @param ... extra parameters for \code{\link{vmatchPattern}} except for 'max.mismatch' since it's calculated internally using the 'qualityThreshold' parameter.
#'
#'
#' @note 
#' \itemize{
#'  \item For qualityThreshold, the alignment score is calculated by (matches*2)-(mismatches+gaps) which programatically translates to round(nchar(patternSeq)*qualityThreshold)*2.
#'  \item No indels are allowed in the function, if expecting indels then use \code{\link{pairwiseAlignSeqs}}.
#'  \item If qualityThreshold is 1, then it is a full match, if 0, then any match is accepted which is useful in searching linker sequences at 3' end. Beware, this function only searches for the pattern sequence in one orientation. If you are expecting to find the pattern in both orientation, you might be better off using BLAST/BLAT!
#'  \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'  before running the function. One can use any of the following 
#'  \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#' 
#' @return IRanges object with starts, stops, and names of the aligned sequences.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}, 
#' \code{\link{doRCtest}}, \code{\link{findAndTrimSeq}}, \code{\link{blatSeqs}}, 
#' \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' subjectSeqs <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC")
#' subjectSeqs <- DNAStringSet(xscat("AAAAAAAAAA", subjectSeqs))
#' vpairwiseAlignSeqs(subjectSeqs, "AAAAAAAAAA", showStats=TRUE)
#' vpairwiseAlignSeqs(subjectSeqs, "AAAAAAAAAA", showStats=TRUE,
#' qualityThreshold=0.5)
vpairwiseAlignSeqs <- function(subjectSeqs=NULL, patternSeq=NULL, side="left", 
                               qualityThreshold=1, showStats=FALSE, 
                               bufferBases=5, doRC=TRUE, parallel=FALSE, ...) {
  dp <- NULL
  
  .checkArgs_SEQed()
  
  qualityThreshold <- as.numeric(qualityThreshold)
  
  ## Only get the relevant side of subject sequence with extra bufferBases to...
  ## account for indels/mismatches & save memory while searching & avoid searching...
  ## elsewhere in the sequence
  if(tolower(side)=="left") {
    badSeqs <- DNAStringSet()
    culprits <- width(subjectSeqs) < (nchar(patternSeq)+bufferBases)
    if(any(culprits)) {
      badSeqs <- subjectSeqs[culprits]
      message(length(badSeqs),
              " sequences were removed from aligning since they were",
              " shorter than pattern getting aligned: ",
              (nchar(patternSeq)+bufferBases),"bp")
      subjectSeqs <- subjectSeqs[!culprits]            
    }
    subjectSeqs2 <- subseq(subjectSeqs, start=1, 
                           end=(nchar(patternSeq)+bufferBases))
    overFromLeft <- rep(0, length(subjectSeqs))
  } else if (tolower(side)=="right") { 
    overFromLeft <- width(subjectSeqs)-(nchar(patternSeq)+bufferBases)
    overFromLeft[overFromLeft<1] <- 1
    subjectSeqs2 <- subseq(subjectSeqs, start=overFromLeft)
  } else {
    subjectSeqs2 <- subjectSeqs
    overFromLeft <- rep(0, length(subjectSeqs))
  }
  
  ## search both ways to test which side yields more hits!        
  if(doRC) {        
    patternSeq <- tryCatch(doRCtest(subjectSeqs2, patternSeq,
                                    qualityThreshold),
                           error=function(e) patternSeq)
  }
  
  if(parallel) {
    subjectSeqs2 <- chunkize(subjectSeqs2)
    maxmis <- round(nchar(patternSeq)*(1-qualityThreshold))
    hits <- bplapply(subjectSeqs2, function(x) {
      bore <- vmatchPattern(patternSeq, x, max.mismatch=maxmis, ...)
      unlist(bore, recursive=TRUE, use.names=TRUE)
    }, BPPARAM=dp) 
    hits <- do.call(c, hits)
  } else {
    max.mm <- round(nchar(patternSeq)*(1-qualityThreshold))
    hits <- vmatchPattern(patternSeq, subjectSeqs2, max.mismatch=max.mm, ...)
    hits <- unlist(hits, recursive=TRUE, use.names=TRUE)
  }
  
  ## test if there were any multiple hits which are overlapping and if they are 
  ## reduce them
  counts <- Rle(names(hits))
  if(any(runLength(counts)>1)) {
    reduced <- reduce(GRanges(seqnames=names(hits),
                              IRanges(start=start(hits), end=end(hits))))
    counts <- seqnames(reduced)
    hits <- ranges(reduced)
    names(hits) <- as.character(seqnames(reduced))
    rm(reduced)
    if(any(runLength(counts)>1)) {
      message("More than 1 pattern (", patternSeq,") match found for:",
              paste(runValue(counts)[runLength(counts)>1],collapse=","),
              "\nUsing the latter occuring hit as the dominant for each read.")
      toremove <- c()
      for(culprits in as.character(runValue(counts)[runLength(counts)>1])) {
        rows <- which(names(hits) %in% culprits)
        toremove <- c(toremove, rows[1:length(rows)-1])
      }
      hits <- hits[-toremove]
      counts <- Rle(names(hits))
      if(any(runLength(counts)>1)) {
        stop("More than 1 pattern unresolved (",patternSeq,") match found for:",
             paste(runValue(counts)[runLength(counts)>1],collapse=","))
      }
    }
  }
  
  good.row <- names(subjectSeqs2) %in% names(hits)
  
  if(showStats) {
    bore <- as.numeric(table(good.row)['FALSE'])    
    message("Total of ", ifelse(is.na(bore),0,bore),
            " did not have the defined pattern sequence (", patternSeq,
            ") that passed qualityThreshold on the ", side, " side")
  }  
  
  starts <- start(hits)
  ends <- end(hits)
  namesq <- names(hits)
  rm("hits","subjectSeqs2")     
  
  hits <- IRanges(start=starts+overFromLeft[good.row]-ifelse(side=="right",2,0), 
                  end=ends+overFromLeft[good.row]-ifelse(side=="right",2,0), 
                  names=namesq)
  
  cleanit <- gc()
  return(hits)
}

#' Test if pattern aligns better in +/- orientation.
#'
#' Given a fixed length pattern sequence and variable length subject sequences, the function roughly finds which orientation of pattern yields the most hits. The function doing the heavylifting is \code{\link{vcountPattern}}. This is an accessory function used in function listed under See Also section below. 
#'
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern. 
#' @param patternSeq DNAString object or a sequence containing the query sequence to search.
#' @param qualityThreshold percent of patternSeq to match. Default is 0.50, half match. This is supplied to max.mismatch parameter of \code{\link{vcountPattern}} as round(nchar(patternSeq)*(1-qualityThreshold)).
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to FALSE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'
#' @return patternSeq that aligned the best. 
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{primerIDAlignSeqs}}
#'
#' @export
#'
#' @examples 
#' subjectSeqs <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC")
#' subjectSeqs <- xscat("AAAAAAAAAA", subjectSeqs)
#' doRCtest(subjectSeqs, "TTTTTTTTT")
doRCtest <- function(subjectSeqs=NULL, patternSeq=NULL, 
                     qualityThreshold=0.5, parallel=TRUE) {
  
  dp <- NULL
  
  .checkArgs_SEQed()
    
  patternSeq.rc <- as.character(reverseComplement(DNAString(patternSeq)))
  hits <- bplapply(c(patternSeq,patternSeq.rc), function(x) {
    max.mm <- round(nchar(x)*(1-qualityThreshold))
    counts <- pmin(vcountPattern(x, subjectSeqs, max.mismatch=max.mm),1)
    sum(counts)
  }, BPPARAM=dp)
  
  if(all(hits==0)) {
    stop("No hits found")
  }
  
  if(hits[[1]] < hits[[2]]) {
    message("There were less/no good hits found for pattern ", patternSeq,
            " than its reverse complement...", patternSeq.rc,
            "\nUsing the latter to perform the searching.")
    patternSeq <- patternSeq.rc
  }
  
  return(patternSeq)
}

# This function generates alignment statistics for findXXXXXX functions. 
# Evaluation of this function happens in the parent function.
.showFindStats <- function() {
  checks <- expression(
    isFindLinker <- grepl("linker",trimmedObj,ignore.case=TRUE) |
      grepl("linker",featureTrim,ignore.case=TRUE),
    listed <- sapply(get(trimmedObj), is.list),
    if(any(listed)) {
      totals <- if(isFindLinker) {        
        sapply(get(trimmedObj)[listed], sapply, function(x) length(x$hits))
      } else {
        sapply(get(trimmedObj)[listed], sapply, length)
      }
      counts <- sapply(get(rawObj)[names(get(trimmedObj)[listed])], 
                       sapply, length)
      stopifnot(identical(colnames(totals), colnames(counts)))
      toprint <- cbind(data.frame("Total"=t(totals)),
                       data.frame("valueTempCol"=t(100*(totals/counts))))
      names(toprint)[3:4] <- gsub("valueTempCol", valueColname, 
                                  names(toprint)[3:4])
      mean.test <- mean(toprint[,paste0(valueColname,
                                        ifelse(isFindLinker,
                                               ".pair2",".pair1"))])<=5      
    } else {
      if(isFindLinker) {
        toprint <- as.data.frame(sapply(sapply(linkerTrimmed[!listed], 
                                               "[[", "hits"), length))
      } else {
        toprint <- as.data.frame(sapply(get(trimmedObj)[!listed], length))
      }
      names(toprint) <- "Total"            
      counts <- sapply(get(rawObj)[names(get(trimmedObj)[!listed])], length)
      toprint[,valueColname] <- 100*(toprint$Total/counts[rownames(toprint)])
      mean.test <- mean(toprint[,valueColname])<=5
    },
    
    toprint$SampleName <- rownames(toprint),
    rownames(toprint) <- NULL,
    
    if(showStats) {            
      message("Sequence Stats after ",featureTrim," alignment:")
      print(toprint)        
    },
    
    ## if <= 5% of sequences found primers...then something is wrong with the 
    ## primer sequences provided???
    if(doTest & mean.test) {
      stop("Something seems to be wrong with the ",featureTrim,
           "s provided for each sample. On average <= 5% of sequences found ",
           ,featureTrim," match for the entire run!!!")
    }
  )
  
  eval.parent(checks)
}

#' Find the 5' primers and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 5' primers for each sample per 
#' sector and adds the results back to the object. This is a specialized 
#' function which depends on many other functions shown in 'see also section' 
#' to perform specialized trimming of 5' primer/adaptor found in the sampleInfo 
#' object. The sequence itself is never trimmed but rather coordinates of primer
#' portion is recorded back to the object and used subsequently by 
#' \code{\link{extractSeqs}} function to perform the trimming.
#'
#' @param sampleInfo sample information SimpleList object outputted from 
#' \code{\link{findBarcodes}}, which holds decoded sequences for samples 
#' per sector/quadrant along with information of sample to primer associations.
#' @param alignWay method to utilize for detecting the primers. One of 
#' following: "slow" (Default), or "fast". Fast, calls 
#' \code{\link{vpairwiseAlignSeqs}} and uses \code{\link{vmatchPattern}} at its 
#' core, which is less accurate with indels and mismatches but much faster. 
#' Slow, calls \code{\link{pairwiseAlignSeqs}} and uses 
#' \code{\link{pairwiseAlignment}} at its core, which is accurate with indels 
#' and mismatches but slower.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern/primer. 
#' Default is FALSE.
#' @param parallel use parallel backend to perform calculation . Defaults to TRUE. 
#' If no parallel backend is registered, then a serial version is ran using
#' \code{\link{SerialParam}}. Parllelization is done at sample level 
#' per sector. Use parallel2 for parallelization at sequence level.
#' @param samplenames a vector of samplenames to process. Default is NULL, which
#' processes all samples from sampleInfo object.
#' @param bypassChecks skip checkpoints which detect if something was odd with 
#' the data? Default is FALSE.
#' @param parallel2 perform parallelization is sequence level. Default is FALSE.
#' Useful in cases where each sector has only one sample with numerous sequences.
#' @param ... extra parameters to be passed to either \code{\link{vmatchPattern}} 
#' or \code{\link{pairwiseAlignment}} depending on 'alignWay' parameter.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new 
#' data added under each sector and sample. New data attributes include: primed
#'
#' @note 
#' \itemize{
#'  \item For paired end data, qualityThreshold for pair 2 is decreased by 
#'  0.10 to increase chances of matching primer sequence.
#'  \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'  before running the function. One can use any of the following 
#'  \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#' \code{\link{extractFeature}}, \code{\link{extractSeqs}}, 
#' \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, 
#' \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples 
#' 
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' findPrimers(seqProps, showStats=TRUE)
#' }
findPrimers <- function(sampleInfo, alignWay="slow", showStats=FALSE, 
                        doRC=FALSE, parallel=TRUE, samplenames=NULL, 
                        bypassChecks=FALSE, parallel2=FALSE, ...) {    
  dp <- NULL
  
  .checkArgs_SEQed()
  
  ## test if there are decoded sequences in the sampleinfo object ##
  decoded <- extractFeature(sampleInfo, feature="decoded")
  samplesDecoded <- sapply(decoded, names, simplify=FALSE)
  sectorsDecoded <- names(which(sapply(sapply(decoded,length), ">", 0)))
  rm(decoded)
  cleanit <- gc()
  
  if(length(sectorsDecoded)==0) {
    stop("No decoded information found in sampleInfo...",
         "did you run findBarcodes()?")
  }
  
  for(sector in sectorsDecoded) {
    message("Processing sector ", sector)
    
    ## refine sample list if specific samples are supplied ##
    samplesToProcess <- samplesDecoded[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    ## prepare sample to primer association ##
    ltrPrimers <- toupper(extractFeature(sampleInfo, sector, samplesToProcess,
                                         feature="primerltrsequence")[[sector]])
    
    ## find any samples which need to be skipped
    skippers <- ltrPrimers=="SKIP"
    if(!all(skippers) & (length(ltrPrimers[!skippers])==0 | 
                           mean(nchar(ltrPrimers[!skippers]))<=5 | 
                           all(is.na(ltrPrimers[!skippers])))) {
      stop("Either the primer size is too short (<=5) or ",
           "no primers are found in sample information object.")
    }

    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
                               feature="primed", 
                               value=structure(rep("SKIPPED",
                                                   length(skip.samples)), 
                                               names=skip.samples))
    }

    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {			
      ## get the decoded reads ##
      decoded <- extractSeqs(sampleInfo, sector, samplesToProcess,
                             feature="decoded")[[sector]]
      
      ## find paired end samples
      isPaired <- extractFeature(sampleInfo, sector, 
                                 feature='pairedend')[[sector]]
      
      ## trim the primers ##
      message("\tFinding Primers.")
      primerIdentity <- extractFeature(sampleInfo, sector,
                                       feature="primerltridentity")[[sector]]
      stopifnot(length(primerIdentity)>0)
      
      primerTrimmed <- bpmapply(function(subjectSeqs, patternSeq, paired,
                                         qualT, ...) {
        if(paired) {
          if(alignWay=="fast") {
            p1 <- vpairwiseAlignSeqs(subjectSeqs$pair1, patternSeq, "left", 
                                     (qualT-0.05), doRC=doRC, 
                                     parallel=parallel2, ...)
          } else if (alignWay=="slow") {
            p1 <- tryCatch(pairwiseAlignSeqs(subjectSeqs$pair1, patternSeq, 
                                             "left", qualT, doRC=doRC, 
                                             parallel=parallel2, ...),
                           error=function(e) geterrmessage())
          }                                   
          
          ## side needs to be middle since we dont 
          ## know where in pair2 the primer can be
          if(alignWay=="fast") {
            p2 <- vpairwiseAlignSeqs(subjectSeqs$pair2, patternSeq, "middle", 
                                     (qualT-0.15), doRC=doRC,  
                                     parallel=parallel2, ...)
          } else if (alignWay=="slow") {
            p2 <- tryCatch(pairwiseAlignSeqs(subjectSeqs$pair2, patternSeq, 
                                             "middle", (qualT-0.1), doRC=doRC, 
                                             parallel=parallel2, ...),
                           error=function(e) geterrmessage())
          }
          list("pair1"=p1, "pair2"=p2)
        } else {
          if(alignWay=="fast") {
            vpairwiseAlignSeqs(subjectSeqs, patternSeq, "left", (qualT-0.05), 
                               doRC=doRC, parallel=parallel2, ...)
          } else if (alignWay=="slow") {
            tryCatch(pairwiseAlignSeqs(subjectSeqs, patternSeq, "left", qualT, 
                                       doRC=doRC, parallel=parallel2, ...),
                     error=function(e) geterrmessage())
          }                                   
        }
      }, decoded[samplesToProcess], ltrPrimers[samplesToProcess],
      isPaired[samplesToProcess], primerIdentity[samplesToProcess],
      SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)
      names(primerTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",primerTrimmed))) {
        stop("Error encountered in LTR Trimming function",
             paste(names(primerTrimmed[grepl("simpleError",primerTrimmed)]),
                   collapse=", "))
      }        
      
      ## remove samples with no primer hits from further processing ##
      culprits <- grep("No hits found",primerTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for primer alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        primerTrimmed <- primerTrimmed[-c(culprits)]
      }
      
      cleanit <- gc()
      
      if(length(primerTrimmed)>0) {
        if(!bypassChecks | showStats) {
          eval(expression(trimmedObj <- "primerTrimmed", rawObj <- "decoded", 
                          featureTrim <- "Primer", 
                          valueColname <- "PercOfDecoded",
                          doTest <- TRUE))
          .showFindStats()
        }
        
        ## modify metadata attribute, write primer coordinates back to 
        ## sampleInfo object & trim
        message("Adding primer info back to the object")
        sampleInfo <- addFeature(sampleInfo, sector, names(primerTrimmed),
                                 feature="primed", value=primerTrimmed)
      }
      
      rm(primerTrimmed)
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory, match.call())
  return(sampleInfo)
}

#' Find the 5' LTRs and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 5' LTR following the primer for
#' each sample per sector and adds the results back to the object. This is a 
#' specialized function which depends on many other functions shown in 'see also
#' section' to perform specialized trimming of 5' viral LTRs found in the 
#' sampleInfo object. The sequence itself is never trimmed but rather 
#' coordinates of LTR portion is added to primer coordinates and recorded back 
#' to the object and used subsequently by \code{\link{extractSeqs}} function to 
#' perform the trimming. This function heavily relies on 
#' \code{\link{pairwiseAlignSeqs}}.
#'
#' @param sampleInfo sample information SimpleList object outputted from 
#' \code{\link{findPrimers}}, which holds decoded and primed sequences for 
#' samples per sector/quadrant along with information of sample to LTR 
#' associations.
#' @param showStats toggle output of search statistics. Default is FALSE. 
#' For paired end data, stats for "pair2" is relative to decoded and/or 
#' primed reads.
#' @param doRC perform reverse complement search of the defined pattern/LTR 
#' sequence. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using 
#' \code{\link{SerialParam}}. Parllelization is done at sample level per 
#' sector.
#' @param samplenames a vector of samplenames to process. Default is NULL, 
#' which processes all samples from sampleInfo object.
#' @param bypassChecks skip checkpoints which detect if something was odd with 
#' the data? Default is FALSE.
#' @param parallel2 perform parallelization is sequence level. Default is FALSE. 
#' Useful in cases where each sector has only one sample with numerous sequences.
#'
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#' @return a SimpleList object similar to sampleInfo paramter supplied with new 
#' data added under each sector and sample. New data attributes include: LTRed
#'
#' @note 
#' \itemize{
#'  \item For paired end data, qualityThreshold for pair 2 is decreased by 
#'  0.05 to increase chances of matching LTR sequence.
#'  \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'  before running the function. One can use any of the following 
#'  \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#' \code{\link{extractFeature}}, \code{\link{extractSeqs}}, 
#' \code{\link{primerIDAlignSeqs}}, \code{\link{findPrimers}}, 
#' \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' findLTRs(seqProps, showStats=TRUE)
#' }
findLTRs <- function(sampleInfo, showStats=FALSE, doRC=FALSE, 
                     parallel=TRUE, samplenames=NULL, bypassChecks=FALSE, 
                     parallel2=FALSE, ...) {    
  dp <- NULL
  
  .checkArgs_SEQed()
  
  ## test if there are primed sequences in the sampleinfo object ##   
  primed <- extractFeature(sampleInfo, feature="primed")
  samplesprimed <- sapply(primed, names, simplify=FALSE)
  sectorsPrimed <- names(which(sapply(sapply(primed, length), ">", 0)))
  rm(primed)
  cleanit <- gc()
  
  if(length(sectorsPrimed)==0) {
    stop("No primed information found in sampleInfo. 
         Did you run findPrimers()?")
  }
  
  for(sector in sectorsPrimed) {
    message("Processing sector ",sector)
    
    ## refine sample list if specific samples are supplied ##
    samplesToProcess <- samplesprimed[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    ## prepare sample to LTR bit associations ##
    sampleLTRbits <- toupper(extractFeature(sampleInfo,sector,samplesToProcess,
                                            feature="ltrbitsequence")[[sector]])
    
    # find any samples which need to be skipped
    skippers <- sampleLTRbits=="SKIP"
    if(!all(skippers) & (length(sampleLTRbits[!skippers])==0 | 
                           mean(nchar(sampleLTRbits[!skippers]))<=1 | 
                           all(is.na(sampleLTRbits[!skippers])))) {
      stop("Either LTR bit sequence is too short (<=1) or no LTR bits ",
           "found in sample information object.")
    }

    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
                               feature="LTRed", 
                               value=structure(rep("SKIPPED",
                                                   length(skip.samples)), 
                                               names=skip.samples))
    }

    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {
      ## get the primer trimmed reads ##
      primerTrimmed <- extractSeqs(sampleInfo, sector, samplesToProcess,
                                   feature="primed")[[sector]]
      
      ## find paired end samples...
      ## add reads which aren't primed in pair2 for cases where reads were long
      isPaired <- extractFeature(sampleInfo, sector, samplesToProcess,
                                 feature='pairedend')[[sector]]
      if(any(isPaired)) {
        message("Getting reads from pair2 which weren't primed...")
        decoded <- extractSeqs(sampleInfo, sector, names(which(isPaired)),
                               feature="decoded", pairReturn='pair2')[[sector]]
        rows <- intersect(names(decoded), names(primerTrimmed))
        
        bore <- mapply(function(x,y) {
          x$pair2 <- c(x$pair2, y[!names(y) %in% names(x$pair2)])
          x
        }, primerTrimmed[rows], decoded[rows])
        primerTrimmed[rows] <- as(bore, "DataFrame")        
        rm(bore)
      }
      
      ## trim LTRbits using slow method since its the best! ##
      message("\tFinding LTR bits.")
      ltrBitIdentity <- extractFeature(sampleInfo,sector=sector,
                                       feature="ltrbitidentity")[[sector]]
      
      ltrTrimmed <- bpmapply(function(subjectSeqs, patternSeq, paired,
                                      qualT, ...) {
        if(paired) {
          p1 <- tryCatch(pairwiseAlignSeqs(subjectSeqs$pair1, patternSeq, 
                                           "left", qualityThreshold=qualT, 
                                           doRC=doRC, parallel=parallel2, ...),
                         error=function(e) geterrmessage())
          
          p2 <- tryCatch(pairwiseAlignSeqs(subjectSeqs$pair2, patternSeq, 
                                           "middle", qualityThreshold=(qualT-0.05), 
                                           doRC=doRC, parallel=parallel2, ...),
                         error=function(e) geterrmessage())
          
          list("pair1"=p1, "pair2"=p2)          
        } else {
          tryCatch(pairwiseAlignSeqs(subjectSeqs, patternSeq, "left", 
                                     qualityThreshold=qualT, doRC=doRC, 
                                     parallel=parallel2, ...),
                   error=function(e) geterrmessage())
        }        
      }, primerTrimmed[samplesToProcess], sampleLTRbits[samplesToProcess],
      isPaired[samplesToProcess], ltrBitIdentity[samplesToProcess],
      SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)
      names(ltrTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",ltrTrimmed))) {
        stop("Error encountered in LTR Trimming function",
             paste(names(ltrTrimmed[grepl("simpleError",ltrTrimmed)]),
                   collapse=", "))
      }
      
      ## remove samples with no LTR hits from further processing ##
      culprits <- grep("No hits found",ltrTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for LTR bit alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        ltrTrimmed <- ltrTrimmed[-c(culprits)]
      }    
      
      cleanit <- gc()
      
      if(length(ltrTrimmed)>0) {
        if(!bypassChecks | showStats) {
          eval(expression(trimmedObj <- "ltrTrimmed", rawObj <- "primerTrimmed", 
                          featureTrim <- "LTR", valueColname <- "PercOfPrimed",
                          doTest <- TRUE))
          .showFindStats()
        }
        
        ## modify metadata attribute, add LTR bit coordinates to primer and write 
        ## back to sampleInfo object & trim...remember everything is relative to 
        ## the entire read length!
        message("Adding LTR info back to the object")
        for(x in names(ltrTrimmed)) {
          cat(".")
          if(isPaired[[x]]) {
            worked <- sapply(ltrTrimmed[[x]],length)>0
            if(any(worked)) {
              for(pair in names(which(worked))) {
                primed <- sampleInfo$sectors[[sector]]$samples[[x]]$primed[[pair]]
                rows <- names(ltrTrimmed[[x]][[pair]]) %in% names(primed)
                primed.end <- end(primed[names(ltrTrimmed[[x]][[pair]][rows])])
                rm(primed)
                
                end(ltrTrimmed[[x]][[pair]][rows]) <- 
                  end(ltrTrimmed[[x]][[pair]][rows]) + primed.end              
                start(ltrTrimmed[[x]][[pair]][rows]) <- 
                  start(ltrTrimmed[[x]][[pair]][rows]) + primed.end              
                rm(primed.end)
              }
            }
            sampleInfo$sectors[[sector]]$samples[[x]]$LTRed <- ltrTrimmed[[x]]
          } else {
            if(length(ltrTrimmed[[x]])>0) {
              primed <- sampleInfo$sectors[[sector]]$samples[[x]]$primed
              primed.end <- end(primed[names(ltrTrimmed[[x]])])
              rm(primed)
              
              end(ltrTrimmed[[x]]) <- end(ltrTrimmed[[x]]) + primed.end
              start(ltrTrimmed[[x]]) <- start(ltrTrimmed[[x]]) + primed.end
              rm(primed.end)
              
              sampleInfo$sectors[[sector]]$samples[[x]]$LTRed <- ltrTrimmed[[x]]
            }
          }        
        }
      }
      rm("ltrTrimmed","primerTrimmed")
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Find the 3' linkers and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 3' linkers for each sample per 
#' sector and adds the results back to the object. This is a specialized 
#' function which depends on many other functions shown in 'see also section' to
#' perform specialized trimming of 3' primer/linker adaptor sequence found in 
#' the sampleInfo object. The sequence itself is never trimmed but rather 
#' coordinates of linker portion is recorded back to the object and used 
#' subsequently by \code{\link{extractSeqs}} function to perform the trimming. 
#' This function heavily relies on either \code{\link{pairwiseAlignSeqs}} or 
#' \code{\link{primerIDAlignSeqs}} depending upon whether linkers getting 
#' aligned have primerID in it or not.
#'
#' @param sampleInfo sample information SimpleList object outputted from 
#' \code{\link{findPrimers}} or \code{\link{findLTRs}}, which holds decoded 
#' sequences for samples per sector/quadrant along with information of sample 
#' to primer associations.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param doRC perform reverse complement search of the defined pattern/linker 
#' sequence. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with
#'  \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#'  registered, then a serial version is ran using 
#'  \code{\link{SerialParam}}. Parllelization is done at sample level 
#'  per sector.
#' @param samplenames a vector of samplenames to process. Default is NULL, 
#' which processes all samples from sampleInfo object.
#' @param bypassChecks skip checkpoints which detect if something was odd 
#' with the data? Default is FALSE.
#' @param parallel2 perform parallelization is sequence level. Default is FALSE. 
#' Useful in cases where each sector has only one sample with numerous sequences.
#'
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#' @return a SimpleList object similar to sampleInfo paramter supplied with new 
#' data added under each sector and sample. New data attributes include: 
#' linkered. If linkers have primerID then, primerIDs attribute is appended 
#' as well. 
#'
#' @note 
#' \itemize{
#'  \item For paired end data, qualityThreshold for pair 2 is increased by 
#'  0.25 or set to 1 whichever is lower to increase quality & full match to 
#'  linker sequence.
#'  \item If no linker matches are found with default options, then try 
#'  doRC=TRUE. 
#'  \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'  before running the function. One can use any of the following 
#'  \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#' \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, 
#' \code{\link{findPrimers}}, \code{\link{extractFeature}}, 
#' \code{\link{extractSeqs}}, \code{\link{findAndTrimSeq}},
#' \code{\link{findIntegrations}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' findLinkers(seqProps, showStats=TRUE, doRC=TRUE)
#' }
findLinkers <- function(sampleInfo, showStats=FALSE, doRC=FALSE, parallel=TRUE, 
                        samplenames=NULL, bypassChecks=FALSE, 
                        parallel2=FALSE, ...) {    
  dp <- NULL
  
  .checkArgs_SEQed()
  
  ## test if there are decoded sequences in the sampleinfo object ##
  decoded <- extractFeature(sampleInfo, feature="decoded")
  toProcessSamples <- sapply(decoded, names, simplify=FALSE)
  sectorsDecoded <- names(which(sapply(sapply(decoded,length), ">", 0)))
  rm(decoded)
  cleanit <- gc()
  
  if(length(sectorsDecoded)==0) {
    stop("No decoded information found in sampleInfo...",
         "did you run findBarcodes()?")
  }
  
  for(sector in sectorsDecoded) {
    message("Processing sector ",sector)
    
    ## subset samplenames from all samples if defined ##
    samplesToProcess <- toProcessSamples[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    ## prepare sample to linker association ##
    sampleLinkers <- toupper(extractFeature(sampleInfo,sector,samplesToProcess,
                                            feature="linkersequence")[[sector]])
    
    # find any samples which need to be skipped
    skippers <- sampleLinkers=="SKIP"
    if(!all(skippers) & (length(sampleLinkers[!skippers])==0 | 
                           mean(nchar(sampleLinkers[!skippers]))<=10 | 
                           all(is.na(sampleLinkers[!skippers])))) {
      stop("Either Linker sequence is too short (<=10) or ",
           "no Linkers found in sample information object.")
    }

    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
                               feature="linkered", 
                               value=structure(rep("SKIPPED",
                                                   length(skip.samples)), 
                                               names=skip.samples))
    }

    ## get the linker quality thresholds for non primerID based samples ##
    linkerIdentity <- extractFeature(sampleInfo, sector, samplesToProcess,
                                     feature="linkeridentity")[[sector]]
    stopifnot(length(linkerIdentity)>0)
    
    ## check if any are primerIDed and get their identity thresholds ##
    primerIded <- extractFeature(sampleInfo, sector, samplesToProcess,
                                 feature="primeridinlinker")[[sector]]        
    primerIded.threshold1 <- 
      extractFeature(sampleInfo, sector, samplesToProcess,
                     feature="primeridinlinkeridentity1")[[sector]]        
    primerIded.threshold2 <- 
      extractFeature(sampleInfo, sector, samplesToProcess,
                     feature="primeridinlinkeridentity2")[[sector]]        

    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {
      
      toProcess <- extractSeqs(sampleInfo, sector, samplesToProcess,
                               feature="decoded")[[sector]]
      
      ## find paired end samples
      isPaired <- extractFeature(sampleInfo,sector,feature='pairedend')[[sector]]
      
      ## trim the Linkers ##
      message("\tFinding Linkers.")
      linkerTrimmed <- bpmapply(function(subjectSeqs, patternSeq, qualT,
                                         paired, pIded, qualT1, qualT2, ...) {        
        if(pIded) {
          if(paired) {
            p1 <- tryCatch(primerIDAlignSeqs(subjectSeqs$pair1, patternSeq, 
                                             doAnchored=TRUE, returnUnmatched=TRUE, 
                                             returnRejected=TRUE, doRC=doRC, 
                                             qualityThreshold1=qualT1, 
                                             qualityThreshold2=qualT2, 
                                             parallel=parallel2, ...),
                           error=function(e) geterrmessage())
            
            p2 <- tryCatch(primerIDAlignSeqs(subjectSeqs$pair2, patternSeq, 
                                             doAnchored=TRUE, returnUnmatched=TRUE, 
                                             returnRejected=TRUE, doRC=doRC, 
                                             qualityThreshold1=pmin(qualT1+.25,1), 
                                             qualityThreshold2=pmin(qualT2+.25,1), 
                                             parallel=parallel2, ...),
                           error=function(e) geterrmessage())
            
            list("pair1"=p1, "pair2"=p2)
          } else {
            tryCatch(primerIDAlignSeqs(subjectSeqs, patternSeq, doAnchored=TRUE, 
                                       returnUnmatched=TRUE, returnRejected=TRUE,
                                       doRC=doRC, qualityThreshold1=qualT1,
                                       qualityThreshold2=qualT2, 
                                       parallel=parallel2, ...),
                     error=function(e) geterrmessage())
          }          
        } else {
          ## use side="middle" since more junk sequence can be present after 
          ## linker which would fail pairwiseAlignSeqs if side='right' for
          ## single end reads or "pair1"
          if(paired) {
            p1 <- tryCatch(pairwiseAlignSeqs(subjectSeqs$pair1, patternSeq, 
                                             "middle", qualityThreshold=qualT, 
                                             returnUnmatched=TRUE, 
                                             returnLowScored=TRUE, doRC=doRC, 
                                             parallel=parallel2, ...),
                           error=function(e) geterrmessage())
            
            # pair2 should end with linker, hence side='right'            
            p2 <- tryCatch(pairwiseAlignSeqs(subjectSeqs$pair2, patternSeq, 
                                             "right", 
                                             qualityThreshold=pmin(qualT+.25,1), 
                                             returnUnmatched=TRUE, 
                                             returnLowScored=TRUE, doRC=doRC, 
                                             parallel=parallel2, ...),
                           error=function(e) geterrmessage())
            
            list("pair1"=p1, "pair2"=p2)
          } else {
            tryCatch(pairwiseAlignSeqs(subjectSeqs, patternSeq, "middle", 
                                       qualityThreshold=qualT, returnUnmatched=TRUE, 
                                       returnLowScored=TRUE, doRC=doRC, 
                                       parallel=parallel2, ...),
                     error=function(e) geterrmessage())
          }
        }        
      },
      toProcess[samplesToProcess], sampleLinkers[samplesToProcess],
      linkerIdentity[samplesToProcess], isPaired[samplesToProcess],
      primerIded[samplesToProcess], primerIded.threshold1[samplesToProcess],
      primerIded.threshold2[samplesToProcess], SIMPLIFY=FALSE, USE.NAMES=FALSE, 
      BPPARAM=dp)      
      names(linkerTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",linkerTrimmed))) {
        stop("Error encountered in Linker Trimming functions",
             paste(names(linkerTrimmed[grepl("simpleError",linkerTrimmed)]),
                   collapse=", "))
      }
      
      ## remove samples with no linker hits from further processing ##
      culprits <- grep("No hits found",linkerTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for Linker alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        linkerTrimmed <- linkerTrimmed[-c(culprits)]
      }        
      
      cleanit <- gc()
      if(length(linkerTrimmed)>0) {
        if(!bypassChecks | showStats) {
          eval(expression(trimmedObj <- "linkerTrimmed", rawObj <- "toProcess", 
                          featureTrim <- "Linker", 
                          valueColname <- "PercOfDecoded",
                          doTest <- TRUE))
          .showFindStats()
        }
        
        message("Adding linker info back to the object")
        ## modify metadata attribute and write back to sampleInfo object
        ## for primerID based samples...write back all the returned attibutes
        for(x in names(linkerTrimmed)) {
          cat(".") 
          if(isPaired[[x]]) {
            for(y in unique(unlist(sapply(linkerTrimmed[[x]], names)))) {
              bore <- sapply(linkerTrimmed[[x]], "[[", y)
              newAttrName <- paste0(ifelse(y=="hits","",y),"linkered")
              sampleInfo$sectors[[sector]]$samples[[x]][[newAttrName]] <- bore
              rm(bore)
            }
          } else {
            for(y in names(linkerTrimmed[[x]])) {
              newAttrName <- paste0(ifelse(y=="hits","",y),"linkered")
              sampleInfo$sectors[[sector]]$samples[[x]][[newAttrName]] <- 
                linkerTrimmed[[x]][[y]]
            }
          }
        }
      }
      rm(linkerTrimmed)
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Find vector DNA in reads and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds vector fragments following the 
#' LTR piece for each sample per sector and adds the results back to the object. 
#' This is a specialized function which depends on many other functions shown
#' in 'see also section' to perform specialized trimming of 5' viral LTRs found
#' in the sampleInfo object. The sequence itself is never trimmed but rather 
#' coordinates of vector portion is added to LTR coordinates and recorded back 
#' to the object and used subsequently by \code{\link{extractSeqs}} function to 
#' perform the trimming. This function heavily relies on \code{\link{blatSeqs}}.
#' In order for this function to work, it needs vector sequence which is read in
#' using 'vectorFile' metadata supplied in the sample information file in 
#' \code{\link{read.sampleInfo}}
#'
#' @param sampleInfo sample information SimpleList object outputted from 
#' \code{\link{findLTRs}}, which holds decoded, primed, and LTRed sequences for 
#' samples per sector/quadrant.
#' @param showStats toggle output of search statistics. Default is FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is  
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#' Parllelization is done at sample level per sector. Use parallel2 for 
#' parallelization at sequence level.
#' @param samplenames a vector of samplenames to process. Default is NULL, which
#' processes all samples from sampleInfo object.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new 
#' data added under each sector and sample. New data attributes include: vectored
#'
#' @note 
#' \itemize{
#'  \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'  before running the function. One can use any of the following 
#'  \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{blatSeqs}},
#' \code{\link{extractFeature}}, \code{\link{extractSeqs}}, 
#' \code{\link{findPrimers}}, \code{\link{findLTRs}}, 
#' \code{\link{findLinkers}}, \code{\link{findAndTrimSeq}},
#' \code{\link{findAndRemoveVector}}
#'
#' @export
#'
#' @examples 
#' 
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' findVector(seqProps, showStats=TRUE)
#' }
findVector <- function(sampleInfo, showStats=FALSE, parallel=TRUE, 
                       samplenames=NULL) {    
  
  dp <- NULL

  .checkArgs_SEQed()

  ## test if there are primed sequences in the sampleinfo object ##   
  primed <- extractFeature(sampleInfo, feature="primed")
  samplesprimed <- sapply(primed, names, simplify=FALSE)
  sectorsPrimed <- names(which(sapply(sapply(primed, length), ">", 0)))
  rm(primed)
  cleanit <- gc()
  
  if(length(sectorsPrimed)==0) {
    stop("No primed information found in sampleInfo. Did you run findPrimers()?")
  }
  
  for(sector in sectorsPrimed) {
    message("Processing sector ",sector)
    
    ## refine sample list if specific samples are supplied ##
    samplesToProcess <- samplesprimed[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    ## get the vector files ##
    vectorFiles <- extractFeature(sampleInfo, sector, samplesToProcess,
                                  feature='vectorfile')[[sector]]
    vectorFiles <- structure(file.path(sampleInfo$sequencingFolderPath, 
                                       vectorFiles), 
                             names=names(vectorFiles))
    
    # find any samples which need to be skipped #
    # gotta use grep here because it could be a file path #
    skippers <- grepl("SKIP$",vectorFiles)
    names(skippers) <- names(vectorFiles)
    if(!all(skippers) & (length(vectorFiles[!skippers])==0 | 
                           mean(nchar(vectorFiles[!skippers]))<=1 | 
                           all(is.na(vectorFiles[!skippers])))) {
      stop("Problem detecting 'vectorFile' association to samples in the
           supplied sampleInfo object.")
    }
    
    skip.samples <- names(which(skippers))
    if(length(skip.samples)>0) {
      samplesToProcess <- samplesToProcess[!samplesToProcess %in% skip.samples]
      message("Skipping samples ", paste(skip.samples,collapse=","))
      sampleInfo <- addFeature(sampleInfo, sector, skip.samples, 
                               feature="vectored", 
                               value=structure(rep("SKIPPED",
                                                   length(skip.samples)), 
                                               names=skip.samples))
    }

    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {
      ## get the LTR trimmed reads ##
      ltrTrimmed <- extractSeqs(sampleInfo, sector, samplesToProcess,
                                feature="LTRed")[[sector]]
      samplesToProcess <- names(ltrTrimmed)
      
      ## find paired end samples...
      ## add reads which aren't primed in pair2 for cases where reads were long
      isPaired <- extractFeature(sampleInfo, sector, samplesToProcess,
                                 feature='pairedend')[[sector]]
      if(any(isPaired)) {
        message("Getting reads from pair2 which weren't LTRed...")
        decoded <- extractSeqs(sampleInfo, sector, names(which(isPaired)),
                               feature="decoded", pairReturn='pair2')[[sector]]
        rows <- intersect(names(decoded), names(ltrTrimmed))
        
        bore <- mapply(function(x,y) {
          x$pair2 <- c(x$pair2, y[!names(y) %in% names(x$pair2)])
          x
        }, ltrTrimmed[rows], decoded[rows])
        ltrTrimmed[rows] <- as(bore, "DataFrame")        
        rm(bore)
      }

      ## find vector bits using BLAT! ##
      message("\tFinding Vector bits.")
      
      vecTrimmed <- mapply(function(subjectSeqs, paired, vectorFile) {
        Vector <- readDNAStringSet(vectorFile)
        if(paired) {
          reads <- subjectSeqs$pair1
          p1 <- findAndRemoveVector(reads, Vector, 10, TRUE, parallel)$hits
          p1 <- if(is(p1,"GRanges")) 
          { IRanges(ranges(p1), names=p1$qNames) } else { p1 }
          
          reads <- subjectSeqs$pair2
          p2 <- findAndRemoveVector(reads, Vector, 10, TRUE, parallel)$hits
          p2 <- if(is(p2,"GRanges"))
          { IRanges(ranges(p2), names=p2$qNames) } else { p2 }
          
          list("pair1"=p1, "pair2"=p2)          
        } else {
          p <- findAndRemoveVector(subjectSeqs, Vector, 10, TRUE, parallel)$hits
          if(is(p,"GRanges"))
          { IRanges(ranges(p), names=p$qNames) } else { p }
        }        
      }, ltrTrimmed[samplesToProcess], isPaired[samplesToProcess], 
      vectorFiles[samplesToProcess], SIMPLIFY=FALSE, USE.NAMES=FALSE)
      names(vecTrimmed) <- samplesToProcess
      
      ## check if any error occured during alignments ##
      if(any(grepl("simpleError",vecTrimmed))) {
        stop("Error encountered in LTR Trimming function",
             paste(names(vecTrimmed[grepl("simpleError",vecTrimmed)]),
                   collapse=", "))
      }
      
      ## remove samples with no LTR hits from further processing ##
      culprits <- grep("No hits found",vecTrimmed)
      if(length(culprits)>0) {
        message("Following sample(s) had no hits for LTR bit alignment: ",
                paste(samplesToProcess[culprits],collapse=", "))
        samplesToProcess <- samplesToProcess[-c(culprits)]
        vecTrimmed <- vecTrimmed[-c(culprits)]
      }    
      
      cleanit <- gc()
      
      if(length(vecTrimmed)>0) {
        if(showStats) {
          eval(expression(trimmedObj <- "vecTrimmed", 
                          rawObj <- "ltrTrimmed", 
                          featureTrim <- "Vector", 
                          valueColname <- "PercOfLTRed",
                          doTest <- FALSE))
          .showFindStats()
        }
        
        ## modify metadata attribute, add vector bit coordinates to LTR 
        # and write back to sampleInfo object & trim...
        # remember everything is relative to the entire read length!
        message("Adding vector info back to the object")
        for(x in names(vecTrimmed)) {
          cat(".")
          if(isPaired[[x]]) {
            worked <- sapply(vecTrimmed[[x]],length)>0
            if(any(worked)) {
              for(pair in names(which(worked))) {
                ltred <- 
                  sampleInfo$sectors[[sector]]$samples[[x]]$LTRed[[pair]]
                rows <- 
                  names(vecTrimmed[[x]][[pair]]) %in% names(ltred)
                ltred.end <- 
                  end(ltred[names(vecTrimmed[[x]][[pair]][rows])])
                rm(ltred)
                
                end(vecTrimmed[[x]][[pair]][rows]) <- 
                  end(vecTrimmed[[x]][[pair]][rows]) + ltred.end              
                start(vecTrimmed[[x]][[pair]][rows]) <- 
                  start(vecTrimmed[[x]][[pair]][rows]) + ltred.end              
                rm(ltred.end)
              }
            }
            sampleInfo$sectors[[sector]]$samples[[x]]$vectored <- 
              vecTrimmed[[x]]
          } else {
            if(length(vecTrimmed[[x]])>0) {
              ltred <- sampleInfo$sectors[[sector]]$samples[[x]]$LTRed
              ltred.end <- end(ltred[names(vecTrimmed[[x]])])
              rm(ltred)
              
              end(vecTrimmed[[x]]) <- 
                end(vecTrimmed[[x]]) + ltred.end
              start(vecTrimmed[[x]]) <- 
                start(vecTrimmed[[x]]) + ltred.end
              rm(ltred.end)
              
              sampleInfo$sectors[[sector]]$samples[[x]]$vectored <-
                vecTrimmed[[x]]
            }
          }        
        }
      }
      rm("vecTrimmed")
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Find the integration sites and add results to SampleInfo object. 
#'
#' Given a SampleInfo object, the function finds integration sites for each 
#' sample using their respective settings and adds the results back to the 
#' object. This is an all-in-one function which aligns, finds best hit per 
#' read per sample, cluster sites, and assign ISU IDs. It calls 
#' \code{\link{blatSeqs}}, \code{\link{read.psl}}, \code{\link{getIntegrationSites}}, 
#' \code{\link{clusterSites}}, \code{\link{otuSites}}. here must be linkered 
#' reads within the sampleInfo object in order to use this function using the 
#' default parameters. If you are planning on BLATing non-linkered reads, 
#' then change the seqType to one of accepted options for the 'feature' 
#' parameter of \code{\link{extractSeqs}}, except for '!' based features.
#'
#' @param sampleInfo sample information SimpleList object outputted from 
#' \code{\link{findLinkers}}, which holds decoded, primed, LTRed, and Linkered 
#' sequences for samples per sector/quadrant along with metadata.
#' @param seqType which type of sequence to align and find integration sites. 
#' Default is NULL and determined automatically based on type of restriction 
#' enzyme or isolation method used. If restriction enzyme is Fragmentase, MuA, 
#' Sonication, or Sheared then this parameter is set to genomicLinkered, else 
#' it is genomic. Any one of following options are valid: genomic, 
#' genomicLinkered, decoded, primed, LTRed, linkered.
#' @param genomeIndices an associative character vector of freeze to full or 
#' relative path of respective of indexed genomes from BLAT(.nib or .2bit files).
#' For example: c("hg18"="/usr/local/blatSuite34/hg18.2bit", "mm8"="/usr/local/blatSuite34/mm8.2bit"). Be sure to supply an index per freeze supplied 
#' in the sampleInfo object. Default is NULL.
#' @param samplenames a vector of samplenames to process. Default is NULL, 
#' which processes all samples from sampleInfo object.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#' @param autoOptimize if aligner='BLAT', then should the blatParameters be 
#' automatically optimized based on the reads? Default is FALSE. When TRUE, 
#' following parameters are adjusted within the supplied blatParameters vector: 
#' stepSize, tileSize, minScore, minIdentity. This parameter is useful 
#' when aligning reads of various lengths to the genome. 
#' Optimization is done using only read lengths. In beta phase!
#' @param doSonic calculate integration sites abundance using breakpoints. 
#' See \code{\link{getSonicAbund}} for more details. Default is FALSE.
#' @param doISU calculate integration site unit for multihits. 
#' See \code{\link{isuSites}} for more details. Default is FALSE.
#' @param ... additional parameters to be passed to \code{\link{blatSeqs}}.
#'
#' @return a SimpleList object similar to sampleInfo parameter supplied with 
#' new data added under each sector and sample. New data attributes include: 
#' psl, and sites. The psl attributes holds the genomic hits per read along 
#' with QC information. The sites attribute holds the condensed integration 
#' sites where genomic hits have been clustered by the Position column and 
#' cherry picked to have each site pass all the QC steps. 
#'
#' @note If parallel=TRUE, then be sure to have a parallel backend registered 
#' before running the function. One can use any of the following 
#' \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{findLTRs}}, 
#' \code{\link{findLinkers}}, \code{\link{startgfServer}}, 
#' \code{\link{read.psl}}, \code{\link{blatSeqs}}, \code{\link{blatListedSet}}, 
#' \code{\link{pslToRangedObject}}, \code{\link{clusterSites}}, 
#' \code{\link{isuSites}}, \code{\link{crossOverCheck}}, 
#' \code{\link{getIntegrationSites}}, \code{\link{getSonicAbund}},
#' \code{\link{annotateSites}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' findIntegrations(seqProps, genomeIndices=c("hg18"="/usr/local/genomeIndexes/hg18.noRandom.2bit"), numServers=2)
#' }
findIntegrations <- function(sampleInfo, seqType=NULL,
                             genomeIndices=NULL, samplenames=NULL,
                             parallel=TRUE, autoOptimize=FALSE, 
                             doSonic=FALSE, doISU=FALSE, ...) {
  
  ## to avoid 'no visible binding for global variable' during checks ##
  freeze <- restrictionenzyme <- startwithin <- alignratiothreshold <- NULL
  genomicpercentidentity <- keepmultihits <- clustersiteswithin <- dp <- NULL
  
  .checkArgs_SEQed()
  
  if(is.null(genomeIndices)) {
    stop("No genomeIndices provided")
  }
  
  ## test if there are linkered sequences in the sampleinfo object if 
  ## specific feature/seqType is not defined ##   
  feature <- ifelse(is.null(seqType), "linkered", seqType)
  message("Checking for ", feature, " reads.")  	
  featured <- extractFeature(sampleInfo, feature=feature)
  samplesfeatured <- sapply(featured, names, simplify=FALSE)
  sectorsfeatured <- names(which(sapply(sapply(featured,length),">",0)))
  rm(featured)
  cleanit <- gc()
  
  if(length(sectorsfeatured)==0) {
    stop("No ", feature, " reads found in sampleInfo object provided.")
  }
  
  ## subset specific samples if defined ##
  samplesToProcess <- unlist(samplesfeatured,use.names=FALSE)
  if(!is.null(samplenames)) {
    samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
  }
  
  message("Creating hashes of settings for aligning and processing.")
  
  ## setup settings hashes for aligning the genomic seqs by species & 
  ## processing hits later
  for(setting in c("restrictionenzyme", "freeze", "startwithin", 
                   "alignratiothreshold", "genomicpercentidentity", 
                   "clustersiteswithin", "keepmultihits")) {
    setter <- extractFeature(sampleInfo, samplename=samplesToProcess, 
                             feature=setting)
    names(setter) <- NULL
    setter <- unlist(setter)
    assign(setting,setter)
  }
  
  minReadLength <- 10

  ## Align by respective species ##
  alignedFiles <- c()        
  for(f in unique(freeze)) {
    message("Aligning to: ", f)
    
    message("Getting sequences to align")
    # get sequences to align #
    if(is.null(seqType)) {
      wanted <- names(restrictionenzyme[!grepl("FRAG|SONIC|MU|SHEAR",
                                               restrictionenzyme,
                                               ignore.case=TRUE)])
      wanted <- wanted[wanted %in% names(freeze[freeze==f])]
      seqs <- extractSeqs(sampleInfo, samplename=wanted, 
                          feature="genomic", strict=TRUE,
                          minReadLength=minReadLength)
      if(any(as.numeric(sapply(seqs,length))>0)) {
        write.listedDNAStringSet(seqs, filePrefix=paste0("processed",f))
      }
      
      wanted <- names(restrictionenzyme[grepl("FRAG|SONIC|MU|SHEAR",
                                              restrictionenzyme,
                                              ignore.case=TRUE)])
      wanted <- wanted[wanted %in% names(freeze[freeze==f])]
      seqs <- extractSeqs(sampleInfo, samplename=wanted, 
                          feature="genomicLinkered", strict=TRUE,
                          minReadLength=minReadLength)
      if(any(as.numeric(sapply(seqs,length))>0)) {
        write.listedDNAStringSet(seqs, filePrefix=paste0("processed",f))
      }
    } else {
      seqs <- extractSeqs(sampleInfo, samplename=names(freeze[freeze==f]), 
                          feature=seqType, minReadLength=minReadLength)
      if(any(as.numeric(sapply(seqs,length))>0)) {
        write.listedDNAStringSet(seqs, filePrefix=paste0("processed",f))
      }
    }

    ## Align seqs ##
    
    queryFiles <- list.files(pattern = paste0("^processed",f))
    
    # merge the extra args with blatSeq() params #
    dots <- list(...)
    blatArgs <- as.list(args(blatSeqs))
    blatArgs <- head(blatArgs, -1) ## remove the last NULL
    
    if(length(dots)>0) {
      for(a in names(dots)) {
        blatArgs[[a]] <- dots[[a]]
      }
    }
    
    blatParams <- eval(blatArgs$blatParameters)
    
    ## autoOptimize BLAT settings if enabled ##
    if(autoOptimize) {        
      # BLAT formula: 2*stepSize+tileSize-1
      queryLengths <- summary(fasta.info(queryFiles, use.names=FALSE))   
      blatParams[["tileSize"]] <- pmax(pmin(queryLengths[["Min."]],15), 8)
      blatParams[["stepSize"]] <- pmax(round((queryLengths[["1st Qu."]]/4)),5)
      blatParams[["minScore"]] <- queryLengths[["Min."]]
      blatParams[["minIdentity"]] <- round(queryLengths[["Median"]]*.95)
      blatArgs[["blatParameters"]] <- blatParams
    }
    
    ## start the requested number of gfServers ##
    numServers <- blatArgs$numServers
    blatArgs[["numServers"]] <- 1L
    gfServerOpts <- c("tileSize","stepSize","minMatch","maxGap", "trans","log",
                      "seqLog","syslog","logFacility","mask","repMatch",
                      "maxDnaHits", "maxTransHits","maxNtSize","maxAsSize",
                      "canStop")
    port <- blatArgs$port + 0:(numServers-1)
    for(n in 1:numServers) {
      searchCMD <- sprintf("gfServer status %s %s", blatArgs$host, port[n])
      if(system(searchCMD, ignore.stderr=TRUE)!=0) {
        message(sprintf("Starting gfServer # %s.", n))
        startgfServer(seqDir=genomeIndices[[f]], host=blatArgs$host, 
                      port=port[n], gfServerOpts=blatParams[names(blatParams) 
                                                            %in% gfServerOpts])
      }
    }
    
    ## align read files in parallel! ##
    blatArgs[["standaloneBlat"]] <- FALSE
    aFile <- bpmapply(function(port, x) {
      blatArgs[["query"]] <- x
      blatArgs[["subject"]] <- genomeIndices[[f]]
      blatArgs[["port"]] <- port
      do.call("blatSeqs", blatArgs)
    }, rep(port, length=length(queryFiles)), queryFiles, 
    SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)
    
    aFile <- unlist(aFile)
    
    ## stop the requested number of gfServers ##
    sapply(port, function(x) stopgfServer(port=x))

    alignedFiles <- c(alignedFiles, aFile)
    message("Cleaning!")
    cleanit <- gc()
    system(paste("rm",paste0("processed",f,".*.fa")))    
  }
  
  ## read all hits and split by samples ##
  message("Reading aligned files.")  
  psl <- read.psl(alignedFiles, bestScoring=TRUE, asGRanges=TRUE, 
                  removeFile=TRUE, parallel=FALSE)    

  mcols(psl)$samplename <- sub("^(.+)-(.+)$","\\1", mcols(psl)$qName)
  psl <- split(psl, mcols(psl)$samplename)
  cleanit <- gc()
  
  ## pair up alignments if sample==paired end ##
  isPaired <- extractFeature(sampleInfo, samplename=samplesToProcess,
                             feature='pairedend')
  names(isPaired) <- NULL
  isPaired <- unlist(isPaired)
  if(any(isPaired)) {
    for(x in names(which(isPaired))) {
      psl[[x]] <- pairUpAlignments(psl[[x]], parallel=parallel)  
    }    
  }
  
  ## begin processing hits ##
  psl.hits <- bplapply(names(psl), function(x) {
    message("Processing ", x)
    
    # add qc info for bestscoring hits #
    psl.x <- 
      getIntegrationSites(psl[[x]], startWithin=startwithin[[x]], 
                          alignRatioThreshold=alignratiothreshold[[x]], 
                          genomicPercentIdentity=genomicpercentidentity[[x]],
                          oneBased=TRUE)
    
    # filter multihits if applicable #
    if(!as.logical(keepmultihits[[x]])) {
      psl.x <- psl.x[!mcols(psl.x)$isMultiHit, ]
    }
    
    # cluster sites by positions #
    psl.x <- clusterSites(psl.rd=psl.x, windowSize=clustersiteswithin[[x]])
    
    # get sonicAbund #
    if(doSonic) {
      psl.x <- getSonicAbund(psl.rd=psl.x)
    }

    # get sites ISU for tagging multihits #
    if(doISU) {        
      psl.x <- isuSites(psl.rd=psl.x)
    }
    
    psl.x
  }, BPPARAM=dp)
  names(psl.hits) <- names(psl)
  
  message("Adding PSL hits back to the object.")
  sampleInfo <- addFeature(sampleInfo, sector=NULL, samplename=names(psl.hits),
                           feature="psl", value=psl.hits)
  
  message("Adding sites back to the object.")
  psl.hits <- sapply(psl.hits, function(x) {
    x <- subset(x, mcols(x)$clusterTopHit & mcols(x)$pass.allQC)
    ranges(x) <- IRanges(mcols(x)$clusteredPosition, width=1) 
    x[, setdiff(colnames(mcols(x)), 
                c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert',
                  'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'tSize', 
                  'blockCount', 'blockSizes', 'qStarts', 'score', 'tStarts', 
                  'pass.startWithin', 'alignRatio', 'pass.alignRatio', 
                  'percIdentity', 'pass.percIdentity', 'pass.allQC', 
                  'clusterTopHit', 'width', 'Position', 'clusteredPosition'))]
  })
  
  sampleInfo <- addFeature(sampleInfo, sector=NULL, samplename=names(psl.hits),
                           feature="sites", value=psl.hits)
  
  cleanit <- gc()
  
  sampleInfo$callHistory <- append(sampleInfo$callHistory,match.call())
  return(sampleInfo)
}

#' Find the 5' primers and add results to SampleInfo object. 
#'
#' Given a sampleInfo object, the function finds 5' primers for each sample per 
#' sector and adds the results back to the object. This is a specialized 
#' function which depends on many other functions shown in 'see also section' 
#' to perform specialized trimming of 5' primer/adaptor found in the sampleInfo 
#' object. The sequence itself is never trimmed but rather coordinates of primer
#' portion is recorded back to the object and used subsequently by 
#' \code{\link{extractSeqs}} function to perform the trimming.
#'
#' @param sampleInfo sample information SimpleList object outputted from 
#' \code{\link{findIntegrations}}, which holds genomic integration sites.
#' @param annots a named list of GRanges object holding features for annotating 
#' integration sites. The name attribute of the list is used as column name.
#' @param samplenames a vector of samplenames to process. Default is NULL, which
#' processes all samples from sampleInfo object.
#' @param parallel use parallel backend to perform calculation. Defaults to TRUE. 
#' If no parallel backend is registered, then a serial version is ran using
#' \code{\link{SerialParam}}. Parllelization is done at sample level 
#' per sector. Use parallel2 for parallelization at sequence level.
#' @param ... additional parameters to be passed to \code{\link{doAnnotation}}
#' except for sites.rd, features.rd, and colnam.
#'
#' @return a SimpleList object similar to sampleInfo paramter supplied with new 
#' data added under each sector and sample. New data attributes include: primed
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{isuSites}}, 
#' \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, 
#' \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' data(genes)
#' genes <- makeGRanges(genes)
#' cpgs <- getUCSCtable("cpgIslandExt","CpG Islands")
#' cpgs <- makeGRanges(cbind(cpgs,strand="*"), chromCol = "chrom")
#' annots <- list("RefGenes"=genes,"CpG"=cpgs)
#' annotateSites(seqProps, annots, annotType="nearest", side="5p")
#' }
annotateSites <- function(sampleInfo, annots=NULL, samplenames=NULL, 
                          parallel=TRUE, ...) {    
  dp <- NULL
  
  .checkArgs_SEQed()
  if(is.null(annots)) {
    stop("annots parameter is empty. Please provide an annotation object.")
  }
  
  ## test if there are sites in the sampleinfo object ##
  sites <- extractFeature(sampleInfo, feature="sites")
  samplesDone <- sapply(sites, names, simplify=FALSE)
  sectorsDone <- names(which(sapply(sapply(sites,length), ">", 0)))
  rm(sites)
  cleanit <- gc()
  
  if(length(sectorsDone)==0) {
    stop("No sites information found in sampleInfo...",
         "did you run findIntegrations()?")
  }
  
  for(sector in sectorsDone) {
    message("Processing sector ", sector)
    
    ## refine sample list if specific samples are supplied ##
    samplesToProcess <- samplesDone[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }

    ## dont bother searching if no samples are to be processed! ##
    if(length(samplesToProcess)>0) {
      
      ## get the sites ##
      sites <- extractFeature(sampleInfo, sector, samplesToProcess,
                              feature="sites")[[sector]]
      
      ## do Annotations ##
      message("\tAnnotating Sites.")
      
      annotated <- bplapply(sites[samplesToProcess], function(x) {
        for(f in names(annots)) {
          x <- doAnnotation(sites.rd=x, features.rd=annots[[f]], colnam=f, ...)
        }
        x
      }, BPPARAM=dp)
      stopifnot(identical(names(annotated),samplesToProcess))
      
      if(length(annotated)>0) {                  
        ## Store annotated sites back to the object ##
        sampleInfo <- addFeature(sampleInfo, sector, names(annotated),
                                 feature="sites", value=annotated)
      }
      
      rm(annotated)
      cleanit <- gc()
    }
  }
  sampleInfo$callHistory <- append(sampleInfo$callHistory, match.call())
  return(sampleInfo)
}

#' Compare LTRed/Primed sequences to all linkers. 
#'
#' Given a SampleInfo object, the function compares LTRed sequences from each sample per sector to all the linker sequences present in the run. The output is a summary table of counts of good matches to all the linkers per sample. 
#'
#' @param sampleInfo sample information SimpleList object outputted from \code{\link{findPrimers}} or \code{\link{findLTRs}}, which holds decoded sequences for samples per sector/quadrant along with information of sample to primer associations.
#' @param qualityThreshold percent of linker length to match, round(nchar(linker)*qualityThreshold). Default is 0.55. Only applies to non-primerID based linkers
#' @param qualityThreshold1 percent of first part of patternSeq to match. Default is 0.75. Only applies to primerID based linker search.
#' @param qualityThreshold2 percent of second part of patternSeq to match. Default is 0.50. Only applies to primerID based linker search.
#' @param doRC perform reverse complement search of the linker sequence. Default is TRUE. Highly recommended!
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'  Parllelization is done at sample level per sector.
#' @param samplenames a vector of samplenames to process. Default is NULL, which processes all samples from sampleInfo object.
#' @param ... extra parameters to be passed to \code{\link{pairwiseAlignment}}.
#'
#' @return a dataframe of counts. 
#'
#' @note If parallel=TRUE, then be sure to have a parallel backend registered 
#' before running the function. One can use any of the following 
#' \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#' \code{\link{primerIDAlignSeqs}}, \code{\link{findLTRs}}, 
#' \code{\link{findPrimers}}, \code{\link{findAndTrimSeq}}
#'
#' @export
troubleshootLinkers <- function(sampleInfo, qualityThreshold=0.55, 
                                qualityThreshold1=0.75, qualityThreshold2=0.50, 
                                doRC=TRUE, parallel=TRUE, samplenames=NULL, 
                                ...) {    
  dp <- NULL
  
  .checkArgs_SEQed()
  
  ## test if there are decoded sequences in the sampleinfo object ##
  toProcess <- extractFeature(sampleInfo, feature="decoded")
  toProcessSamples <- sapply(toProcess, names, simplify=FALSE)
  sectors <- names(toProcess)   
  rm(toProcess)
  cleanit <- gc()
  
  results <- data.frame()
  for(sector in sectors) {
    message("Processing sector ",sector)
    
    ## prepare sample to linker association ##
    sampleLinkers <- toupper(extractFeature(sampleInfo, sector=sector,
                                            feature="linkersequence")[[sector]])
    if(length(sampleLinkers)==0 | mean(nchar(sampleLinkers))<=10 | 
         all(is.na(sampleLinkers))) {
      stop("Either Linker sequence is too short (<=10) or ",
           "no Linkers found in sample information object.")
    }
    
    samplesToProcess <- toProcessSamples[[sector]]
    if(!is.null(samplenames)) {
      samplesToProcess <- samplesToProcess[samplesToProcess %in% samplenames]
    }
    
    toProcess.seqs <- extractSeqs(sampleInfo,sector,samplesToProcess,
                                  feature="decoded")[[sector]]    
    
    ## find paired end samples
    isPaired <- extractFeature(sampleInfo, sector, feature='pairedend')[[sector]]
    
    ## do all by all comparison of linkers ##
    message("\tFinding Linkers.")
    for(linkerSeq in unique(as.character(sampleLinkers))) {
      message("Checking ",linkerSeq)
      
      linkerTrimmed <- bplapply(samplesToProcess, function(x) {
        if(length(unlist(gregexpr("N",linkerSeq)))>3) {
          if(isPaired[[x]]) {
            p1 <- try_default(length(primerIDAlignSeqs(
              toProcess.seqs[[x]]$pair1, linkerSeq, qualityThreshold1, 
              qualityThreshold2, doRC=doRC, ...)$hits), 0, quiet=TRUE)
            
            p2 <- try_default(length(primerIDAlignSeqs( 
              toProcess.seqs[[x]]$pair2, linkerSeq, qualityThreshold1, 
              qualityThreshold2, doRC=doRC, ...)$hits), 0, quiet=TRUE)
            
            list("pair1"=p1, "pair2"=p2)
          } else {
            try_default(length(primerIDAlignSeqs(
              toProcess.seqs[[x]], linkerSeq, qualityThreshold1, 
              qualityThreshold2, doRC=doRC, ...)$hits), 0, quiet=TRUE)
          }
        } else {
          if(isPaired[[x]]) {
            p1 <- try_default(length(pairwiseAlignSeqs(
              toProcess.seqs[[x]]$pair1, linkerSeq, "middle", qualityThreshold, 
              doRC=doRC, ...)), 0, quiet=TRUE)
            
            p2 <- try_default(length(pairwiseAlignSeqs(
              toProcess.seqs[[x]]$pair2, linkerSeq, "right", qualityThreshold, 
              doRC=doRC, ...)), 0, quiet=TRUE)
            
            list("pair1"=p1, "pair2"=p2)                       
          } else {
            try_default(length(pairwiseAlignSeqs(
              toProcess.seqs[[x]], linkerSeq, "middle", qualityThreshold,
              doRC=doRC, ...)), 0, quiet=TRUE)
          }
        }
      }, BPPARAM=dp)      
      names(linkerTrimmed) <- samplesToProcess
      
      isPaired <- sapply(linkerTrimmed, is.list)
      if(any(isPaired)) {
        totalSeqs <- sapply(toProcess.seqs[isPaired], sapply, length)
        linkerhits <- sapply(linkerTrimmed[isPaired], unlist)
        PercentOfTotal <- linkerhits/totalSeqs[,names(linkerTrimmed[isPaired])]        
      } else {
        totalSeqs <- sapply(toProcess.seqs[!isPaired], length)
        linkerhits <- as.numeric(unlist(linkerTrimmed[!isPaired]))
        PercentOfTotal<- linkerhits/totalSeqs[names(linkerTrimmed[!isPaired])]
      }
      
      bore <- data.frame("linkerSeq"=linkerSeq, 
                         "samplename"=names(linkerTrimmed), 
                         "linkerhits"=linkerhits, 
                         "PercentOfTotal"=PercentOfTotal,
                         stringsAsFactors=FALSE)
      rownames(bore) <- NULL
      results <- rbind(results, bore)
      
      cleanit <- gc()
    }        
  }  
  sampleLinkers <- extractFeature(sampleInfo, feature="linkersequence")
  names(sampleLinkers) <- NULL
  sampleLinkers <- unlist(sampleLinkers)
  linkersample <- as.data.frame(sampleLinkers)
  linkersample <- tapply(rownames(linkersample), linkersample$sampleLinkers, 
                         paste, collapse=",")
  
  results$CorrectLinker <- with(results,
                                sampleLinkers[as.character(samplename)]==
                                  as.character(linkerSeq))
  results$CorrectSample <- with(results, linkersample[as.character(linkerSeq)])
  return(results)
}

#' Find and trim a short pattern sequence from the subject. 
#'
#' This function facilitates finding and trimming of a short pattern sequence from a collection of subject sequences. The trimming is dictated by side parameter. For more information on the trimming process see the 'side' parameter documentation in \code{\link{trimSeqs}}. For information regarding the pattern alignment see the documentation for \code{\link{pairwiseAlignSeqs}}. This function is meant for aligning a short pattern onto large collection of subjects. If you are looking to align a vector sequence to subjects, then please use BLAT.
#'
#' @param patternSeq DNAString object or a sequence containing the query sequence to search.
#' @param subjectSeqs DNAStringSet object containing sequences to be searched for the pattern.
#' @param side which side of the sequence to perform the search & trimming: left, right or middle. Default is 'left'.
#' @param offBy integer value dictating if the trimming base should be offset by X number of bases. Default is 0.
#' @param alignWay method to utilize for detecting the primers. One of following: "slow" (Default), "fast", or "blat". Fast, calls \code{\link{vpairwiseAlignSeqs}} and uses \code{\link{vmatchPattern}} at its core, which is less accurate with indels and mismatches but much faster. Slow, calls \code{\link{pairwiseAlignSeqs}} and uses \code{\link{pairwiseAlignment}} at its core, which is accurate with indels and mismatches but slower. Blat will use \code{\link{blatSeqs}}.
#' @param ... parameters to be passed to \code{\link{pairwiseAlignment}}, \code{\link{vpairwiseAlignSeqs}} or \code{\link{blatSeqs}} depending on which method is defined in 'alignWay' parameter.
#'
#' @return DNAStringSet object with pattern sequence removed from the subject sequences. 
#'
#' @note If parallel=TRUE, then be sure to have a parallel backend registered 
#' before running the function. One can use any of the following 
#' \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{extractFeature}}, \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}, \code{\link{findPrimers}}, \code{\link{findLinkers}}
#'
#' @export
#'
#' @examples 
#' findAndTrimSeq(patternSeq="AGACCCTTTT",
#' subjectSeqs=DNAStringSet(c("AGACCCTTTTGAGCAGCAT","AGACCCTTGGTCGACTCA",
#' "AGACCCTTTTGACGAGCTAG")), qualityThreshold=.85, doRC=FALSE, side="left", 
#' offBy=1, alignWay = "slow")
#'
findAndTrimSeq <- function(patternSeq, subjectSeqs, side = "left", offBy = 0, 
                           alignWay = "slow", ...) {
  
  # to avoid NOTE during R cmd check #
  removeSubjectNamesAfter <- NULL
  
  .checkArgs_SEQed()
  
  coords <- switch(alignWay,
                   fast = vpairwiseAlignSeqs(subjectSeqs, patternSeq, side,...),
                   slow = pairwiseAlignSeqs(subjectSeqs, patternSeq, side, ...),
                   blat = with(read.psl(blatSeqs(subjectSeqs, patternSeq, ...)),
                               IRanges(qStart,qEnd))
  )
  
  res <- trimSeqs(subjectSeqs, coords, side, offBy)
  if(removeSubjectNamesAfter) {
    names(res) <- NULL
  }
  res
}

#' Find and trim vector sequence from reads. 
#'
#' This function facilitates finding and trimming of long/short fragments of 
#' vector present in LM-PCR products. The algorithm looks for vector sequence 
#' present anywhere within the read and trims according longest contiguous 
#' match on either end of the read. Alignment is doing using BLAT
#'
#' @param reads DNAStringSet object containing sequences to be trimmed for vector.
#' @param Vector DNAString object containing vector sequence to be searched in reads.
#' @param minLength integer value dictating minimum length of trimmed sequences to return. Default is 10.
#' @param returnCoords return the coordinates of vector start-stop for the matching reads. Defaults to FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'
#' @return DNAStringSet object with Vector sequence removed from the reads. If returnCoords=TRUE, then a list of two named elements "hits" & "reads". The first element, "hits" is a GRanges object with properties of matched region and whether it was considered valid denoted by 'good.row'. The second element, "reads" is a DNAStringSet object with Vector sequence removed from the reads.
#'
#' @note If parallel=TRUE, then be sure to have a parallel backend registered 
#' before running the function. One can use any of the following 
#' \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, \code{\link{pslToRangedObject}}, \code{\link{blatSeqs}}, \code{\link{read.blast8}}, \code{\link{findAndTrimSeq}}
#' @export
findAndRemoveVector <- function(reads, Vector, minLength=10, 
                                returnCoords=FALSE, parallel=TRUE) {
  
  .checkArgs_SEQed()
  
  ## use tryCatch to be safe incase not reads were aligned in read.psl stops
  ## with "No hits found"
  hits <- tryCatch(read.psl(blatSeqs(query=reads, subject=Vector, 
                                     parallel=parallel,
                                     blatParameters=c(stepSize=3,tileSize=8,
                                                      minIdentity=70,minScore=5,
                                                      repMatch=112312,out="psl")
                                     ), bestScoring=FALSE),
                   error = function(z) NULL)
  
  ## delete temporary files like subjectFile.fa.*.tempyS generated and not
  ## removed after the command above due to an abrupt stop ##
  toremove <- list.files(pattern="subjectFile.fa.+.tempyS")
  if(length(toremove)>0) {
    file.remove(toremove)
  }
  
  if(!is.null(hits)) {
    hits <- reduce(pslToRangedObject(hits, useTargetAsRef=FALSE),
                   min.gapwidth=10, ignore.strand=TRUE)
    hits$qSize <- width(reads[as.character(seqnames(hits))])
    
    ## only consider hits that start or end with vector...
    ## others would be wierdos! ##
    tocheck <- start(hits) <= 10 | hits$qSize-end(hits) <= 11
    hits <- hits[tocheck]
  }
  
  if(length(hits)>0) {
    ## queries where genomic sequence is surrounded by vector sequence (counts>1)...
    ## we will fix these later!
    hits$qNames <- as.character(seqnames(hits))
    counts <- table(hits$qNames)
    hits$toCure <- hits$qNames %in% names(counts[counts>1])
    rm(counts)
    
    ### queries composed of recombined forms of vector
    widthSum <- tapply(width(hits), hits$qNames, sum)
    hits$widthSum <- widthSum[hits$qNames]
    rm(widthSum)    
    hits$allCovered <- hits$qSize >= (hits$widthSum-3) & 
      hits$qSize <= (hits$widthSum+3)
    
    ### seqs w/ vector fragments which are not 'allCovered' with vector
    toTrim <- reads[names(reads) %in% hits$qNames[!hits$allCovered]]
    
    ### seqs w/o any vector sequences
    reads <- reads[!names(reads) %in% hits$qNames] 
    
    ### queries where genomic sequence is either at the beginning or the end
    hits$qEndDiff <- hits$qSize-end(hits)
    hits$StartDiff <- start(hits)-0
    cured <- hits[!hits$toCure]
    
    # take the side which has the higher unmatched bases to the vector sequence #
    cured$trueStart <- ifelse(cured$StartDiff > cured$qEndDiff, 0, end(cured)+1)
    cured$trueEnd <- ifelse(cured$StartDiff > cured$qEndDiff, 
                            start(cured)-1, cured$qSize)
    cured$type <- "genomic either side"
    
    ### fix queries where genomic sequence is surrounded by vector sequence
    if(any(hits$toCure)) {
      hits.list <- split(hits[hits$toCure], hits$qNames[hits$toCure])
      bores <- bplapply(hits.list, function(bore) {        
        trueRange <- gaps(ranges(bore))
        bore$trueStart <- start(trueRange[width(trueRange)== 
                                            max(width(trueRange))
                                          ])
        bore$trueEnd <- end(trueRange[width(trueRange)==max(width(trueRange))])
        bore$type <- "genomic in middle"
        bore$good.row <- TRUE
        bore[1]        
      })  
      bores <- unlist(GRangesList(bores), use.names=FALSE)
      cured <- c(cured, bores[,names(mcols(cured))])
      rm("hits.list","bores")
    }
    rm("hits")
    
    atStart <- cured$StartDiff<=10
    atEnd <- cured$qEndDiff<=10
    
    ### remove sequences where vector sequence is in the middle 
    cured$toCure[!atStart & !atEnd] <- FALSE 
    cured$type[!atStart & !atEnd] <- "vector in middle" 
    
    ## extract sequence with trueStart & trueEnd 
    good.row <- !cured$toCure & !cured$allCovered
    trimmed <- subseq(toTrim[as.character(seqnames(cured))[good.row]],
                      start=cured$trueStart[good.row]+1,
                      end=cured$trueEnd[good.row])
    
    ### combined seqs.good with trimmed ###
    reads <- c(reads, trimmed[width(trimmed)>=minLength])
    
    if(returnCoords) {
      list("hits"=cured, "reads"=reads)
    } else {
      reads
    }
  } else {
    message("No vector hits found")
    if(returnCoords) {
      list("hits"="No hits found", "reads"=reads)
    } else {
      reads
    }
  }
}
#' Trim sequences from a specific side.
#'
#' This function trims a DNAStringSet object using the ranges from left, right, or middle of the sequence. This is a helper function utilized in \code{\link{primerIDAlignSeqs}} and \code{\link{extractSeqs}}. If dnaSet and coords are not the same length, then they are required to have a names attribute to perform the matched trimming. 
#'
#' @param dnaSet DNAStringSet object containing sequences to be trimmed.
#' @param coords IRanges object containing boundaries.
#' @param side either 'left','right',or the Default 'middle'.
#' @param offBy integer value dictating if the supplied coordinates should be offset by X number of bases. Default is 0.
#'
#' @return a DNAStringSet object with trimmed sequences. 
#'
#' @note If side is left, then any sequence following end of coords+offBy is returned. If side is right, then sequence preceding start of coords-offBy is returned. If side is middle, then sequence contained in coords is returned where offBy is added to start and subtracted from end in coords.
#'
#' @seealso \code{\link{extractSeqs}}, \code{\link{primerIDAlignSeqs}}
#'
#' @export
#'
#' @examples 
#' dnaSet <- DNAStringSet(c("AAAAAAAAAACCTGAATCCTGGCAATGTCATCATC", 
#' "AAAAAAAAAAATCCTGGCAATGTCATCATCAATGG", "AAAAAAAAAAATCAGTTGTCAACGGCTAATACGCG",
#' "AAAAAAAAAAATCAATGGCGATTGCCGCGTCTGCA", "AAAAAAAAAACCGCGTCTGCAATGTGAGGGCCTAA",
#' "AAAAAAAAAAGAAGGATGCCAGTTGAAGTTCACAC"))
#' coords <- IRanges(start=1, width=rep(10,6))
#' trimSeqs(dnaSet, coords, side="left", offBy=1)
#' trimSeqs(dnaSet, coords, side="middle")
trimSeqs <- function(dnaSet, coords, side="middle", offBy=0) {
  stopifnot(class(dnaSet) %in% c("DNAStringSet", "DNAString"))
  stopifnot(class(coords)=="IRanges")
  
  if(length(dnaSet)==0 | length(coords)==0) {
    stop("dnaSet/coords is empty. Please supply reads/coords to be trimmed.")
  }
  
  # check if both dnaSet and coords has 'names' attribute, 
  # if yes then check if they have matching names, else check lengths. 
  if(is.null(names(dnaSet)) | is.null(names(coords))) {
    stopifnot(length(dnaSet)==length(coords))
  } else {
    rows <- match(names(coords), names(dnaSet))
    if(any(is.na(rows))) {
      stop("Some of the names in coords are not present in dnaSet")
    }
    if(!is.ordered(rows)) {
      dnaSet <- dnaSet[rows]
      if(!identical(names(dnaSet), names(coords))) {
        stop("Names are not identical between dnaSet and coords parameter")
      }
    }
  }        
  
  # temp helper function to show messages #
  .showMessage <- function(x) {
    message("Following sequences were removed from trimming since their ",
            "coordinates+offBy were out of sequence length: ", 
            paste(x,collapse=", "))
  }
  
  # trim by side and check if any of the coords are off the sequence length in dnaSet
  if(tolower(side)=="left") {
    test <- end(coords)+offBy > width(dnaSet) | end(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], start=end(coords[!test])+offBy)
  } else if (tolower(side)=="right") {
    test <- start(coords)-offBy > width(dnaSet) | end(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], end=start(coords[!test])-offBy)
  } else {
    test <- start(coords)+offBy > width(dnaSet) | 
      end(coords)-offBy > width(dnaSet) | start(coords)+offBy < 1
    if(any(test)) {
      .showMessage(names(dnaSet)[test])
    }
    subseq(dnaSet[!test], 
           start=start(coords[!test])+offBy, 
           end=end(coords[!test])-offBy)
  }    
}

#' Extract sequences for a feature in the sampleInfo object.
#'
#' Given a sampleInfo object, the function extracts sequences for a defined 
#' feature.
#'
#' @param sampleInfo sample information SimpleList object, which samples per 
#' sector/quadrant information along with other metadata.
#' @param sector specific sector to extract sequences from. Default is NULL, 
#' which extracts all sectors. 
#' @param samplename specific sample to extract sequences from. Default is NULL, 
#' which extracts all samples. 
#' @param feature which part of sequence to extract (case sensitive). Options 
#' include: primed, !primed, LTRed, !LTRed, linkered, !linkered, primerIDs, 
#' genomic, genomicLinkered, decoded, and unDecoded. If a sample was primerIDed 
#' and processed by \code{\link{primerIDAlignSeqs}}, then all the rejected and
#' unmatched attributes can be prepended to the feature. Example: vectored,
#' Rejectedlinkered, RejectedprimerIDslinkered, Absentlinkered, or 
#' unAnchoredprimerIDslinkered. When feature is genomic, it includes sequences 
#' which are primed, LTRed, linkered, and !linkered. The genomicLinkered is same 
#' as genomic minus the !linkered. When feature is decoded, it includes 
#' everything that demultiplexed. The '!' in front of a feature extracts the 
#' inverse. One can only get unDecoded sequences if returnUnmatched was TRUE in
#' \code{\link{findBarcodes}}. If \code{\link{findVector}} was run and 
#' "vectored" feature was found in the sampleInfo object, then genomic & 
#' genomicLinkered output will have vectored reads removed.
#' @param trim whether to trim the given feature from sequences or keep it. 
#' Default is TRUE. This option is ignored for feature with '!'.
#' @param minReadLength threshold for minimum length of trimmed sequences to return.
#' @param sideReturn if trim=TRUE, which side of the sequence to return: left, 
#' middle, or right. Defaults to NULL and determined automatically. Doesn't 
#' apply to features: decoded, genomic or genomicLinkered.
#' @param pairReturn if the data is paired end, then from which pair to return 
#' the feature. Options are "pair1", "pair2", or defaults to "both". Ignored if 
#' data is single end. 
#' @param strict this option is used when feature is either 'genomic' or
#' 'genomicLinkered'. When a sample has no LTRed reads, primer ends are used as 
#' starting points by default to extract the genomic part. Enabling this option
#' will strictly ensure that only reads with primer and LTR are trimmed for the
#' 'genomic' or 'genomicLinkered' feature. Default is FALSE.
#'
#' @return a listed DNAStringSet object structed by sector then sample. 
#' Note: when feature='genomic' or 'genomicLinkered' and when data is paired end, 
#' then "pair2" includes union of reads from both pairs which found LTR.
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{findLTRs}}, 
#' \code{\link{findLinkers}}, \code{\link{trimSeqs}}, \code{\link{extractFeature}},
#' \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' samples <- c('Roth-MLV3p-CD4TMLVWell6-Tsp509I', 
#' 'Roth-MLV3p-CD4TMLVWell6-MseI', 'Roth-MLV3p-CD4TMLVwell5-MuA')
#' extractSeqs(seqProps, sector='2', samplename=samples, feature="primed")
#' extractSeqs(seqProps, sector='2', samplename=samples, feature="!primed")
#' extractSeqs(seqProps, sector='2', samplename=samples, feature="linkered")
#' extractSeqs(seqProps, sector='2', samplename=samples, feature="genomic")
extractSeqs <- function(sampleInfo, sector=NULL, samplename=NULL, 
                        feature="genomic", trim=TRUE, minReadLength=1, 
                        sideReturn=NULL, pairReturn="both", strict=FALSE) {
  
  .checkArgs_SEQed()
  
  # get all sectors and samplenames in each sector or a specific sector
  res <- getSectorsForSamples(sampleInfo, sector, samplename)
  sectors <- res[["sectors"]]
  samplenames <- res[["samplenames"]]
  
  if(tolower(feature)=="undecoded") {
    sapply(sectors, function(y) { 
      allmetadata <- metadata(sampleInfo$sectors[[y]])
      if("unDecodedSeqs" %in% names(allmetadata)) {
        allmetadata$unDecodedSeqs 
      } else {
        message("No unDecoded attribute found for the supplied ",
                "sampleInfo object and sector ", y)
      }
    })
  } else {        
    res <- sapply(sectors, function(y) {
      sapply(samplenames[[y]], function(x,y) { 
        decoded <- sampleInfo$sectors[[y]]$samples[[x]]$decoded
        isPaired <- sampleInfo$sectors[[y]]$samples[[x]]$pairedend
        if(isPaired & pairReturn!='both') {
          decoded <- decoded[[pairReturn]]
        }
        
        if(tolower(feature)!="decoded") {
          # get ride of ! from feature, else R wont know what to do when making a... 
          # new object with ! in the front.
          assign(gsub("!","",feature), 
                 sampleInfo$sectors[[y]]$samples[[x]][[gsub("!","",feature)]])
        }
        
        if(tolower(feature)=="decoded") {
          decoded
        } else if (tolower(feature) %in% c("genomic","genomiclinkered")) {
          primed <- sampleInfo$sectors[[y]]$samples[[x]]$primed
          LTRed <- sampleInfo$sectors[[y]]$samples[[x]]$LTRed
          linkered <- sampleInfo$sectors[[y]]$samples[[x]]$linkered
          vectored <- sampleInfo$sectors[[y]]$samples[[x]]$vectored
          
          if(isPaired & pairReturn!='both') {
            primed <- primed[[pairReturn]]
            LTRed <- LTRed[[pairReturn]]
            linkered <- linkered[[pairReturn]]
            vectored <- vectored[[pairReturn]]
          }
          
          # if reads dont have LTRs...use primers instead...but throw a warning!
          LTRed.test <- if(is.list(LTRed)) {
            any(sapply(LTRed,is.null))
          } else {
            is.null(LTRed)
          }
          
          if(LTRed.test & strict) {
            warning("LTRed information not found for ",x, " skipping...", 
                    immediate.=TRUE)
            return()
          } else if (LTRed.test & !strict) {
            warning("LTRed information not found for ",x,
                    "...using primer end as starting boundary.", 
                    immediate.=TRUE)
          }
          
          # check vectored attribute is not null #
          vec.test <- if(is.list(vectored)) {
            any(sapply(vectored,is.null))
          } else {
            is.null(vectored)
          }
          if(vec.test) {
            rm(vectored)                                
          }
          
          # test if there are any primed and/or linkered reads
          p.l.test <- if(is.list(primed)) {
            any(sapply(primed,is.null), sapply(linkered,is.null))
          } else {
            any(is.null(primed), is.null(linkered))
          }
          
          if(p.l.test) { 
            message("No sequences found for requested feature (", feature,
                    ") for sample: ", x,"...skipping.") 
          } else {
            if(trim) {
              # get all LTRed reads and make ends = size of each read
              if(isPaired & pairReturn=='both') {
                if(LTRed.test & !strict) {
                  LTRed <- primed
                }
                
                if(exists("vectored")) {                  
                  toremove <- unlist(sapply(vectored, names), use.names=FALSE)
                  message("Removing ",length(toremove)," vectored reads...")
                  LTRed <- lapply(LTRed, function(p) p[!names(p) %in% toremove])
                }
                
                starts <- lapply(LTRed, function(p) structure(end(p), 
                                                              names=names(p)))
                ## add reads present pair1 to pair2 where LTR wasn't found which
                ## makes things consistent to pair1 atleast!
                loners <- setdiff(names(starts$pair1), names(starts$pair2))
                starts$pair2 <- c(starts$pair2, 
                                  structure(rep(0L,length(loners)), names=loners))
                ends <- lapply(decoded, function(p) structure(width(p), 
                                                              names=names(p)))
                
                stopifnot(identical(names(starts), names(ends)))
                stopifnot(mapply(function(s,e) all(names(s) %in% names(e)),
                                 starts, ends))
                
                coords <- mapply(function(s,e) {
                  IRanges(start=s+1, end=e[names(s)], names=names(s))
                }, starts, ends)
                
                # alter ends for reads where linker was present
                # fix cases where start of linker earlier than edge of LTR/primer end
                ends <- lapply(linkered, function(p) structure(start(p), 
                                                               names=names(p)))
                stopifnot(identical(names(coords), names(ends)))
                
                coords <- mapply(function(p, e) {
                  there <- intersect(names(p), names(e))
                  e <- as.numeric(e[there]-1)
                  culprits <- start(p[there]) > e
                  end(p[there][!culprits]) <- e[!culprits]
                  end(p[there][culprits]) <- start(p[there][culprits])
                  p
                }, coords, ends)
                
                if(tolower(feature)=="genomiclinkered") { 
                  stopifnot(identical(names(coords), names(linkered)))
                  coords <- mapply(function(p, l) p[names(p) %in% names(l)], 
                                   coords, linkered)
                }
                
                # trim it and return non zero length sequences
                if(any(sapply(coords, length)>0)) {
                  stopifnot(identical(names(coords), names(decoded)))
                  mapply(function(d, p) {
                    seqs <- trimSeqs(d,p)
                    seqs[width(seqs)>=minReadLength]
                  }, decoded, coords)
                } else {
                  message("No linkered reads found for sample: ",x,"...skipping.")
                }
              } else {
                if(LTRed.test & !strict) {
                  LTRed <- primed
                }
                
                if(exists("vectored")) {                  
                  message("Removing ",length(vectored)," vectored reads...")
                  LTRed <- LTRed[!names(LTRed) %in% names(vectored)]
                }
                
                starts <- structure(end(LTRed), names=names(LTRed))
                ends <- structure(width(decoded), names=names(decoded))
                
                stopifnot(identical(names(starts), names(ends[names(starts)])))
                coords <- IRanges(start=starts+1,
                                  end=ends[names(starts)],
                                  names=names(starts))
                
                # alter ends for reads where linker was present
                # fix cases where start of linker earlier than edge of LTR/primer end
                ends <- structure(start(linkered), names=names(linkered))
                there <- intersect(names(coords), names(ends))
                ends <- as.numeric(ends[there]-1)
                culprits <- start(coords[there]) > ends
                end(coords[there][!culprits]) <- ends[!culprits]
                end(coords[there][culprits]) <- start(coords[there][culprits])
                
                if(tolower(feature)=="genomiclinkered") { 
                  coords <- coords[names(coords) %in% names(linkered)] 
                }
                
                # trim it and return non zero length sequences
                if(length(coords)>0) {
                  seqs <- trimSeqs(decoded,coords)
                  seqs[width(seqs)>=minReadLength]
                } else {
                  message("No linkered reads found for sample: ",x,"...skipping.")
                }                
              }
            } else {
              # just returning ranges...simply match names from decoded to the request
              if(isPaired & pairReturn=='both') {
                if(tolower(feature)=="genomiclinkered") { 
                  mapply(function(d,l) d[names(d) %in% names(l)],
                         decoded, linkered)
                } else { ## everything past primer and LTR
                  if(LTRed.test & !strict) {
                    mapply(function(d,l) d[names(d) %in% names(l)],
                           decoded, primed)
                  } else {
                    mapply(function(d,l) d[names(d) %in% names(l)],
                           decoded, LTRed)
                  }
                }
              } else {
                if(tolower(feature)=="genomiclinkered") { 
                  decoded[names(decoded) %in% names(linkered)]
                } else { ## everything past primer and LTR
                  if(is.null(LTRed)) {
                    decoded[names(decoded) %in% names(primed)]
                  } else {
                    decoded[names(decoded) %in% names(LTRed)]
                  }
                }
              }              
            }
          }
        } else {
          notFeature <- grepl("!",feature)
          ## no need to trim if looking for "not" based feature since there are 
          ## no coordinates for it
          if(notFeature) {
            if(isPaired & pairReturn=='both') {
              stopifnot(identical(names(decoded), 
                                  names(get(gsub("!","",feature)))))
              toreturn <- mapply(function(d,g) d[!names(d) %in% names(g)], 
                                 decoded, get(gsub("!","",feature)))
              
              toreturn <- list()
              toreturn[["pair1"]] <- 
                decoded$pair1[!names(decoded$pair1) %in% 
                                names(get(gsub("!","", feature))$pair1)]
              toreturn[["pair2"]] <- 
                decoded$pair2[!names(decoded$pair2) %in% 
                                names(get(gsub("!","", feature))$pair2)]
            } else {
              toreturn <- decoded[!names(decoded) %in% 
                                    names(get(gsub("!","",feature)))]
            }
            
            if(length(toreturn)) {
              toreturn
            } else {
              message("No sequences found for requested feature (", feature,
                      ") for sample: ", x,"...skipping.") 
            }
          } else {                    
            if(is.null(get(feature))) { 
              message("No sequences found for requested feature (", feature,
                      ") for sample: ", x,"...skipping.") 
            } else {
              if(isPaired & pairReturn=='both') {
                stopifnot(identical(names(decoded), names(get(feature))))                
                res.seq <- list("pair1"=
                                  decoded$pair1[names(decoded$pair1) %in% 
                                                  names(get(feature)$pair1)],
                                "pair2"=
                                  decoded$pair2[names(decoded$pair2) %in% 
                                                  names(get(feature)$pair2)])
              } else {
                res.seq <- decoded[names(decoded) %in% names(get(feature))]  
              }
              
              if(trim) {
                if(is.null(sideReturn)) {
                  sidetype <- ifelse(grepl("primerID",feature,ignore.case=TRUE), 
                                     "middle", 
                                     ifelse(grepl("primed|LTRed",feature,
                                                  ignore.case=TRUE), 
                                            "left", 
                                            ifelse(grepl("linkered",feature,
                                                         ignore.case=TRUE), 
                                                   "right", 
                                                   "middle")))
                } else {
                  sidetype <- tolower(sideReturn)
                }                            
                
                offByLength <- ifelse(sidetype=="middle",0,1)
                
                if(isPaired & pairReturn=='both') {
                  seqs.p1 <- trimSeqs(res.seq$pair1, get(feature)$pair1, 
                                      side=sidetype, offBy=offByLength)
                  seqs.p1 <- seqs.p1[width(seqs.p1)>=minReadLength]
                  
                  seqs.p2 <- trimSeqs(res.seq$pair2, get(feature)$pair2, 
                                      side=sidetype, offBy=offByLength)
                  seqs.p2 <- seqs.p2[width(seqs.p2)>=minReadLength]
                  
                  list("pair1"=seqs.p1, "pair2"=seqs.p2)
                } else {
                  seqs <- trimSeqs(res.seq, get(feature), 
                                   side=sidetype, offBy=offByLength)
                  seqs[width(seqs)>=minReadLength]
                }                
              } else {
                res.seq
              }
            }
          }
        }
      }, y=y)
    }, simplify=FALSE)
    
    simpletons <- !sapply(res,class)=="matrix"
    if(any(simpletons)) {
      lengthTest <- lapply(lapply(lapply(res[simpletons], 
                                         function(x) sapply(x, length)),">",0),
                           which)
      res <- mapply(function(x,y) x[y], res[simpletons], lengthTest, 
                    SIMPLIFY=FALSE)  
    }
    
    if(any(!simpletons)) {
      res[!simpletons] <- sapply(res[!simpletons], as, 'DataFrame')
    }
    
    res
  }
}

#' Extract a specific feature/attribute of the sampleInfo object.
#'
#' Given a sampleInfo object, the function extracts a defined feature(s) for given sample or sector.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param sector a vector or specific sector to extract the feature from. Default is NULL, which extracts all sectors. 
#' @param samplename a character vector or a specific sample to extract feature from. Default is NULL, which extracts all samples. 
#' @param feature Options include: primed, LTRed, linkered, decoded, and any of the metadata. Default is NULL. When feature='metadata', then it returns names of all the metadata elements associated with the sample as a comma separated list.
#'
#' @return a list or list of lists depending upon which parameters were supplied.
#'
#' @seealso \code{\link{addFeature}}, \code{\link{findPrimers}},
#'  \code{\link{findLTRs}}, \code{\link{findLinkers}}, 
#'  \code{\link{extractSeqs}}, \code{\link{trimSeqs}}, 
#'  \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' samples <- c('Roth-MLV3p-CD4TMLVWell6-Tsp509I', 
#' 'Roth-MLV3p-CD4TMLVWell6-MseI', 'Roth-MLV3p-CD4TMLVwell5-MuA')
#' extractFeature(seqProps, sector='2', samplename=samples, feature="primed")
#' extractFeature(seqProps, sector='2', samplename=samples, feature="linkered")
#' extractFeature(seqProps, sector='2', samplename=samples, feature="metadata")
extractFeature <- function(sampleInfo, sector=NULL, samplename=NULL, 
                           feature=NULL) {
  
  .checkArgs_SEQed()
  
  # get all sectors and samplenames in each sector or a specific sector
  res <- getSectorsForSamples(sampleInfo, sector, samplename)
  sectors <- res[["sectors"]]
  samplenames <- res[["samplenames"]]
  
  res <- sapply(sectors, function(y) {
    sapply(samplenames[[y]], function(x,y) { 
      if(tolower(feature)=="metadata") {
        paste(names(sampleInfo$sectors[[y]]$samples[[x]]), collapse=", ")
      } else {
        res <- sampleInfo$sectors[[y]]$samples[[x]][[feature]]
        ## convert any factor based vector to the appropriate regular vector
        if(class(res)=="factor") { 
          if(!any(is.na(suppressWarnings(as.numeric(levels(res)))))) { 
            as.numeric(as.character(res)) } 
          else { 
            as.character(res) 
          }    
        } else if (class(res)=="character"){ 
          ## if a numeric vector is stored as character, convert it back to numeric
          if(!any(is.na(suppressWarnings(as.numeric(res))))) { 
            as.numeric(res) 
          } else { 
            res 
          }
        } else {
          res
        }
      }
    }, y=y)
  }, simplify=FALSE)
  
  simpletons <- !sapply(res,class)=="matrix"
  if(any(simpletons)) {
    lengthTest <- lapply(lapply(lapply(res[simpletons], 
                                       function(x) sapply(x, length)),">",0),
                         which)
    res <- mapply(function(x,y) x[y], res[simpletons], lengthTest, 
                  SIMPLIFY=FALSE)  
  }
  
  if(any(!simpletons)) {
    res[!simpletons] <- sapply(res[!simpletons], as, 'DataFrame')
  }
  
  res
}

#' Add a specific feature/attribute to the sampleInfo object.
#'
#' Given a sampleInfo object, the function adds a new feature for the 
#' given samples & sectors.
#'
#' @param sampleInfo sample information SimpleList object, which samples per 
#' sector/quadrant information along with other metadata.
#' @param sector a vector or a specific sector to add the new feature(s) to. 
#' Default is NULL, in which case the sectors are searched from 
#' samplename parameter.
#' @param samplename a character vector or a specific sample to add the new 
#' feature(s) to. Default is NULL.
#' @param feature a string of naming the new feature to add for the defined 
#' samplename and sector.
#' @param value named vector of samplenames & values which is assigned for the 
#' defined sector, samplename, and feature. Example: c("Sample1"="ACDTDASD")
#'
#' @return modified sampleInfo object with new feature(s) added.
#'
#' @seealso \code{\link{findPrimers}}, \code{\link{extractSeqs}}, 
#' \code{\link{trimSeqs}}, \code{\link{extractFeature}}, 
#' \code{\link{getSectorsForSamples}}
#'
#' @export
#'
#' @examples 
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' extractFeature(seqProps, sector="2", 
#' samplename="Roth-MLV3p-CD4TMLVWell6-MseI", feature="metadata")
#' seqProps <- addFeature(seqProps, sector="2", 
#' samplename="Roth-MLV3p-CD4TMLVWell6-MseI", feature="foo", 
#' value=c("Roth-MLV3p-CD4TMLVWell6-MseI"="woo"))
#' extractFeature(seqProps, sector="2", 
#' samplename="Roth-MLV3p-CD4TMLVWell6-MseI", feature="metadata")
addFeature <- function(sampleInfo, sector=NULL, samplename=NULL, feature=NULL, 
                       value=NULL) {
  
  .checkArgs_SEQed()
  
  if(is.null(value)) {
    stop("Please define the value parameter.")
  }
  
  if(is.null(samplename)) {
    stop("Please define the samplename to add the features to.")
  }
  
  if(!all(samplename %in% names(value))) {
    stop("Not all samplename(s) are found in value parameter")
  }

  # get all sectors and samplenames in each sector or a specific sector
  res <- getSectorsForSamples(sampleInfo, sector, samplename)
  sectors <- res[["sectors"]]
  samplenames <- res[["samplenames"]]
  
  for(y in sectors) {
    for(x in samplenames[[y]]) {
      cat(".")
      sampleInfo$sectors[[y]]$samples[[x]][[feature]] <- value[[x]]
    }
  }
  
  sampleInfo$callHistory <- append(sampleInfo$callHistory, match.call())
  return(sampleInfo)
}

#' Get sectors for samples defined in the sampleInfo object.
#'
#' Given a sampleInfo object, the function gets the sectors for each samplename. This is an accessory function utilized by other functions of this package to aid sector retrieval.
#'
#' @param sampleInfo sample information SimpleList object, which samples per sector/quadrant information along with other metadata.
#' @param sector a specific sector or vector of sectors if known ahead of time. Default is NULL, which extracts all sectors. 
#' @param samplename a specific sample or vector of samplenames to get sectors for. Default is NULL, which extracts all samples. 
#' @param returnDf return results in a dataframe. Default is FALSE.
#'
#' @return If returnDf=TRUE, then a dataframe of sector associated with each samplename, else a named list of length two: x[["sectors"]] and x[["samplenames"]]
#'
#' @seealso \code{\link{extractSeqs}}, \code{\link{extractFeature}}, 
#' \code{\link{addFeature}}
#'
#' @export
#'
#' @examples 
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' samples <- c('Roth-MLV3p-CD4TMLVWell6-Tsp509I', 
#' 'Roth-MLV3p-CD4TMLVWell6-MseI', 'Roth-MLV3p-CD4TMLVwell5-MuA')
#' getSectorsForSamples(seqProps, samplename=samples)
#' getSectorsForSamples(seqProps, samplename=samples, returnDf=TRUE)
getSectorsForSamples <- function(sampleInfo, sector=NULL, samplename=NULL,
                                 returnDf=FALSE) {
  
  .checkArgs_SEQed()
  
  # get all sectors and samplenames in each sector or a specific sector if defined
  if(is.null(sector)) {
    sectors <- names(sampleInfo$sectors)        
  } else {
    sectors <- sector
  }
  samplenames <- sapply(sectors,
                        function(x) names(sampleInfo$sectors[[x]]$samples),
                        simplify=FALSE)
  
  # if specific samplename(s) is defined, then search to find where they are
  if(is.null(samplename)) {
    samplename <- unlist(samplenames)     
  }
  
  if(!all(samplename %in% unlist(samplenames))) {            
    stop("Following sample(s) do not exist on given sector(s) (",
         paste(sectors,collapse=", "),") in the supplied sampleInfo object: ",
         paste(samplename[!samplename %in% unlist(samplenames)],collapse=", "))
  }
  sectors <- names(which(unlist(lapply(lapply(samplenames,"%in%",samplename), 
                                       any)))) 
  if(returnDf) {
    return(do.call(rbind, lapply(sectors, function(x) { 
      data.frame(samplename=samplenames[[x]][samplenames[[x]] %in% samplename], 
                 sector=x, stringsAsFactors=FALSE)
    })))
  } else {
    samplenames <- sapply(sectors, 
                          function(x) 
                            samplenames[[x]][samplenames[[x]] %in% samplename],
                          simplify=FALSE)
    return(list("sectors"=sectors, "samplenames"=samplenames))
  }
}

#' Read fasta/fastq/sff given the path or sampleInfo object.
#'
#' Given a sequence reads file path, the function returns a DNAStringSet object.
#'
#' @param seqFilePath a path to fasta/fastq/sff reads file or a sampleInfo 
#' object returned by \code{\link{read.SeqFolder}}
#' @param sector specific sector to reads sequences from. Default is 1, and not 
#' required if seqFilePath is a direct file path rather than sampleInfo object.
#' @param isPaired does the sector contain paired end reads? Default is FALSE
#'
#' @return if isPaired is FALSE, then a DNAStringSet object, else a list of 
#' DNAStringSet objects of three elements corresponding to reads from 
#' "barcode", "pair1", and "pair2". Note: "pair2" is reverse complemented!
#'
#' @seealso \code{\link{findBarcodes}}, \code{\link{read.SeqFolder}}, 
#' \code{\link{extractSeqs}}
#'
#' @export
#'
#' @examples 
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' read.seqsFromSector(seqProps, sector="2")
#' }
read.seqsFromSector <- function(seqFilePath=NULL, sector=1, isPaired=FALSE) {
  if(is.null(seqFilePath)) {
    stop("Missing seqFilePath!")
  }
  
  if(class(seqFilePath)=="SimpleList") {
    seqfilePattern <- seqFilePath$seqfilePattern
    seqFilePaths <- seqFilePath$seqFilePaths
    
    if(isPaired) {
      filePath <- seqFilePaths[grep(paste0("^",gsub("R1|R2|I1",".*",sector),
                                           seqfilePattern),
                                    basename(seqFilePaths))]
      filePath <- normalizePath(filePath, mustWork=TRUE)
      
      if(length(filePath)==0) {
        stop("No sequence file found for sector: ", sector,
             " in seqFilePath variable (",
             paste(seqFilePaths,collapse=" \n"), ") using pattern (",
             paste0(sector,seqfilePattern),")")
      }
      
      nobarcodes <- FALSE ## incase the data is already demultiplexed
      if(!any(grepl("I1", basename(filePath)))) {
        warning("No index/barcode file (I1) found.")
        nobarcodes <- TRUE
      }
      
      pair2 <- paste0(sub("(.*)R\\d.*","\\1",sector),
                      ifelse(sub(".*(R\\d).*","\\1",sector)=="R2","R1","R2"))
      if(!any(grepl(pair2, basename(filePath)))) {
        stop("Pair #2 or Linker end read file not present!")
      }
      seqFilePath <- filePath
      
    } else {
      filePath <- seqFilePaths[grep(paste0("^",sector,seqfilePattern),
                                    basename(seqFilePaths))]
      filePath <- normalizePath(filePath, mustWork=TRUE)
      
      if(length(filePath)==0) {
        stop("No sequence file found for sector: ", sector,
             " in seqFilePath variable (",
             paste(seqFilePaths,collapse=" \n"), ") using pattern (",
             paste0(sector,seqfilePattern),")")
      }
      
      if(length(filePath)>1) {
        stop("Multiple sequence file found for sector: ", sector,
             " in seqFilePath variable (", paste(seqFilePaths,collapse=" \n"),
             ") using pattern (", paste0(sector,seqfilePattern),")")
      }
      seqFilePath <- filePath
    }
  }
  
  doSingleEndCheck <- FALSE
  
  message("Reading:\n", paste(seqFilePath, collapse="\n"))
  if(any(grepl("fastq",seqFilePath,ignore.case=TRUE))) {    
    dnaSet <- sapply(seqFilePath, function(x) {
      dnaSet <- readDNAStringSet(x, format='fastq')      
      names(dnaSet) <- sub("^\\S+-(\\S+) .+$", "\\1", names(dnaSet), perl=TRUE)
      if(any(duplicated(names(dnaSet)))) {
        stop("Duplicate definition lines found in file: ", x)
      }
      if(length(dnaSet)==0) {
        stop("No sequences found in file: ", x)
      }
      dnaSet
    })
    
    if(isPaired & length(seqFilePath)>1) {
      LTRside <- grep(sector, basename(names(dnaSet)))
      #names(dnaSet[[LTRside]]) <- paste0("@pair1side@",names(dnaSet[[LTRside]]))
      
      linkerSide <- grep(pair2, basename(names(dnaSet)))      
      #names(dnaSet[[linkerSide]]) <- paste0("@pair2side@",names(dnaSet[[linkerSide]]))
      
      if(!nobarcodes) {
        barcodes <- grep("I1", basename(names(dnaSet)))
        dnaSet <- list("barcode"=dnaSet[[barcodes]],
                       "pair1"=dnaSet[[LTRside]],
                       "pair2"=reverseComplement(dnaSet[[linkerSide]]))
      } else {
        dnaSet <- list("pair1"=dnaSet[[LTRside]],
                       "pair2"=reverseComplement(dnaSet[[linkerSide]]))
      }
      
    } else {
      ## for single-end data!
      dnaSet <- dnaSet[[1]]
    }
  } else if (any(grepl("sff",seqFilePath,ignore.case=TRUE))) {    
    dnaSet <- rSFFreader::sread(rSFFreader::readSff(seqFilePath, 
                                                    use.qualities=FALSE, 
                                                    verbose=FALSE))
    doSingleEndCheck <- TRUE    
  } else {
    dnaSet <- readDNAStringSet(seqFilePath)
    doSingleEndCheck <- TRUE    
  }
  
  if(doSingleEndCheck) {
    if(any(duplicated(names(dnaSet)))) {
      stop("Duplicate definition lines found in file(s): ", 
           paste(seqFilePath, collapse=" * "))
    }
    
    if(length(dnaSet)==0) {
      stop("No sequences found")
    }
  }
  
  return(dnaSet)
}

#' Write a fasta file per sample in parallel
#'
#' Given a listed DNAStringSet object return from \code{\link{extractSeqs}}, the
#' function writes a fasta file for each sample as defined in filePath parameter.
#'
#' @param dnaSet listed DNAStringSet object containing sequences to be written.
#' @param filePath a path write the fasta files per sample. Default is current
#' working directory.
#' @param filePrefix prefix the filenames with a string. Default is 'processed' 
#' followed by samplename.
#' @param prependSamplenames Prepend definition lines with samplenames. 
#' Default is TRUE. Make sure the dnaSet parameter is a named list where 
#' names are used as samplenames.
#' @param format either fasta (the default) or fastq.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'
#' @seealso \code{\link{findBarcodes}}, \code{\link{read.SeqFolder}}, 
#' \code{\link{extractSeqs}}, \code{\link{addListNameToReads}}
#'
#' @note
#' \itemize{
#'   \item Writing of the files is done using \code{\link{writeXStringSet}} 
#'   with parameter append=TRUE. This is to aggregate reads from a sample 
#'   which might be present in more than one sector. 
#'   \item If data is paired end, then each pair will be written separately 
#'   with designations in the filename as well as in the definition line as 
#'   (at)pairX(at) appended at the end.
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered 
#'   before running the function. One can use any of the following 
#'   \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#' }
#'
#' @export
#' 
#' @examples
#'  
#' \donttest{
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' samples <- c('Roth-MLV3p-CD4TMLVWell6-Tsp509I', 
#' 'Roth-MLV3p-CD4TMLVWell6-MseI', 'Roth-MLV3p-CD4TMLVwell5-MuA')
#' seqs <- extractSeqs(seqProps, sector='2', samplename=samples, feature="primed")
#' write.listedDNAStringSet(seqs)
#' }
write.listedDNAStringSet <- function(dnaSet, filePath=".", 
                                     filePrefix="processed", 
                                     prependSamplenames=TRUE, format="fasta", 
                                     parallel=FALSE) {
  stopifnot(class(dnaSet)=="list")
  
  if(filePrefix=="") { filePrefix <- NA }
  filePath <- normalizePath(filePath, mustWork=TRUE)
  
  dp <- if(parallel) { bpparam() } else { SerialParam() }
  for(x in seq_along(dnaSet)) {
    bpmapply(function(outputSeqs, samplename) {
      if(is.list(outputSeqs)) {
        list.test <- any(sapply(outputSeqs, length)>0)
      } else {
        list.test <- length(outputSeqs)>0
        outputSeqs <- list("tempy"=outputSeqs)
      }
      
      if(list.test) {
        for(p in names(outputSeqs)) {
          pairname <- ifelse(p=="tempy", NA, p)
          
          if(is.null(names(outputSeqs[[p]]))) {
            message("No names attribute found for ", samplename,
                    " ... using artifically generated names")
            names(outputSeqs[[p]]) <- 
              paste("read",1:length(outputSeqs[[p]]),sep="-")
          }
          
          ## remove '.' at the beginning of the filename incase filePrefix is empty
          filename <- paste(na.omit(c(filePrefix, samplename, pairname, 
                                      ifelse(format=="fastq", "fastq", "fa"))), 
                            collapse=".")
          filename <- paste(filePath, filename, sep="/")
          
          if(prependSamplenames) {
            names(outputSeqs[[p]]) <- paste(samplename, names(outputSeqs[[p]]), 
                                            sep="-")
          }
          
          if(p!="tempy") {
            names(outputSeqs[[p]]) <- paste0(names(outputSeqs[[p]]), 
                                             "@",pairname,"@")
          }
          
          writeXStringSet(outputSeqs[[p]], filepath=filename, format=format, 
                          append=TRUE) 
        } 
      } else {
        message("No reads written for ", samplename)
      }       
    }, as(dnaSet[[x]],"list"), names(dnaSet[[x]]), BPPARAM=dp)    
  }
}

# Check args and set defaults for functions dealing with reads. This function
# checks all the arguments passed to a function related to aligning or trimming 
# a read and then sets default values for internal use. Evaluation of this 
# function happens in the parent function.
.checkArgs_SEQed <- function() {
  
  checks <- expression( 
    if("subjectSeqs" %in% names(formals())) { 
      if(is.null(subjectSeqs) | length(subjectSeqs)==0) {
        stop("subjectSeqs paramter is empty. Please supply reads to be aligned")
      }      
      
      ## give names if not there for troubleshooting purpose in later steps
      removeSubjectNamesAfter <- FALSE
      if(is.null(names(subjectSeqs))) {
        removeSubjectNamesAfter <- TRUE
        names(subjectSeqs) <- paste("read", 1:length(subjectSeqs))
      } 
    },
    
    if("patternSeq" %in% names(formals())) { 
      if(is.null(patternSeq) | length(patternSeq)==0) {
        stop("patternSeq paramter is empty. Please supply reads to be aligned")
      } else if (length(patternSeq)>1) {
        stop("More than 1 patternSeq defined. Please only supply one pattern.")
      }
    },
    
    if("reads" %in% names(formals())) { 
      if(is.null(reads) | length(reads)==0) {
        stop("reads paramter is empty. Please supply reads to be aligned")
      }
    },
    
    if("Vector" %in% names(formals())) { 
      if(is.null(Vector) | length(Vector)==0) {
        stop("Vector paramter is empty. Please supply Vector to be aligned")
      }
      
      if(!grepl("DNAString",class(Vector))) {
        stop("Vector paramter is not of DNAString class")
      }
    },
    
    if("parallel" %in% names(formals())) { 
      dp <- if(parallel & .Platform$OS.type != "windows") { bpparam() } else { SerialParam() }
    },
    
    if("sampleInfo" %in% names(formals())) { 
      stopifnot(is(sampleInfo,"SimpleList"))
    },
    
    if("feature" %in% names(formals())) {
      if(is.null(feature)) {
        stop("Please define a feature to extract.")
      }
    }      
  )
  
  eval.parent(checks)
}

# Check args and set defaults for functions dealing with hits of aligned reads.
# This function checks all the arguments passed to a function related to 
# checking hits of an aligned read and then sets default values for internal 
# use. Evaluation of this function happens in the parent function.
.checkArgsSetDefaults_ALIGNed <- function() {
  
  checks <- expression(
    if("pslFile" %in% names(formals()) | "files" %in% names(formals())) { 
      if(is.null(files) | length(files)==0) {
        stop("files parameter empty. Please supply a filename to be read.")
      }
      
      if(any(grepl("\\*|\\$|\\+|\\^",files))) {
        ## vector of filenames
        files <- list.files(path=dirname(files), pattern=basename(files), 
                            full.names=TRUE)      
      }
      
      if(length(files)==0) { 
        stop("No file(s) found with given paramter in files:", files) 
      }
    },

    if("psl.rd" %in% names(formals())) {
      if(!is.null(psl.rd)) {
        if(length(psl.rd)==0 | !is(psl.rd,"GRanges")) {
          stop("psl.rd paramter is empty or not a GRanges object")
        }
      }
    },

    if("parallel" %in% names(formals())) { 
      dp <- if(parallel & .Platform$OS.type != "windows") { bpparam() } else { SerialParam() }
    }
  )
  
  eval.parent(checks)
}

#' Reads a BAM/SAM file and converts it into a PSL like format. 
#'
#' Given filename(s), the function reads the BAM/SAM file, converts into a PSL 
#' like format. Any other file format will yield errors or erroneous results.
#' This is intended to be used independently with other short read aligners.
#'
#' @param bamFile BAM/SAM filename, or vector of filenames, or a pattern of 
#' files to import.
#' @param removeFile remove the file(s) supplied in bamFile paramter after 
#' importing. Default is FALSE.
#' @param asGRanges return a GRanges object. Default is TRUE
#'
#' @return a GRanges or GAlignments object reflecting psl file type.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{blatSeqs}}, 
#' \code{\link{read.blast8}}, \code{\link{read.psl}}, 
#' \code{\link{pslToRangedObject}}, \code{\link{pairUpAlignments}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' read.BAMasPSL(bamFile="processed.*.bam$")
#' read.BAMasPSL(bamFile=c("sample1hits.bam","sample2hits.bam"))
#' }
read.BAMasPSL <- function(bamFile=NULL, removeFile=TRUE, asGRanges=TRUE) {
  if(is.null(bamFile) | length(bamFile)==0) {
    stop("bamFile parameter empty. Please supply a filename to be read.")
  }
  
  if (any(grepl("\\*|\\$|\\+|\\^",bamFile))) {
    ## vector of filenames
    bamFile <- list.files(path=dirname(bamFile), 
                          pattern=basename(bamFile), full.names=TRUE)      
  }
  
  if(length(bamFile)==0) { 
    stop("No file(s) found with given paramter in bamFile:", bamFile) 
  }
  
  if("hits" %in% names(bamFile)) {
    bamFile <- bamFile[["hits"]]
  }
  
  param <- ScanBamParam(what=c("qname"),
                        flag=scanBamFlag(isPaired = NA, isProperPair = NA, 
                                         isUnmappedQuery = FALSE, 
                                         hasUnmappedMate = NA,
                                         isMinusStrand = NA, 
                                         isMateMinusStrand = NA,
                                         isFirstMateRead = NA, 
                                         isSecondMateRead = NA,
                                         isSecondaryAlignment = NA, 
                                         isDuplicate = NA,
                                         isNotPassingQualityControls = NA),
                        tag=c("CC", "CT", "CP", "CG"))
  hits <- lapply(bamFile, readGAlignments, param=param)
  hits <- do.call(c, hits)
  
  if(length(hits)==0) {
    if(removeFile) { file.remove(bamFile) }
    stop("No hits found")
  }
  
  ## remove extra tag columns if not present ##
  for(f in c("CC", "CT", "CP", "CG")) {
    if(all(is.na(mcols(hits)[,f]))) {
      mcols(hits)[[f]] <- NULL
    }
  }
  
  message("1. Ordering by qName")
  mcols(hits)$qName <- mcols(hits)$qname
  mcols(hits)$qname <- NULL
  hits <- hits[order(mcols(hits)$qName)]
  
  ## add few required columns present in PSL format which are crucial to filtering
  message("Adding relevent PSL columns")
  stopifnot(identical(cigarWidthAlongQuerySpace(cigar(hits)), qwidth(hits)))
  mcols(hits)$qSize <- qwidth(hits)
  
  bore <- cigarRangesAlongQuerySpace(cigar(hits), ops="M")
  mcols(hits)$matches <- sapply(width(bore), sum)
  mcols(hits)$qStart <- min(start(bore))
  mcols(hits)$qEnd <- max(end(bore))
  
  bore <- cigarRangesAlongQuerySpace(cigar(hits), ops="X")
  mcols(hits)$misMatches <- sapply(width(bore), sum)
  
  bore <- cigarRangesAlongQuerySpace(cigar(hits), ops="I")
  mcols(hits)$qNumInsert <- sapply(bore, length)
  mcols(hits)$qBaseInsert <- sapply(width(bore), sum)
  
  bore <- cigarRangesAlongReferenceSpace(cigar(hits), ops="I")
  mcols(hits)$tNumInsert <- sapply(bore, length)
  mcols(hits)$tBaseInsert <- sapply(width(bore), sum)
  
  if(asGRanges) {
    message("Converting to GRanges object...")
    hits.gr <- as(hits, "GRanges")
    mcols(hits.gr) <- cbind(DataFrame(cigar=cigar(hits), ngap=njunc(hits)),
                            mcols(hits))
    hits <- hits.gr
    rm(hits.gr)    
  }
  
  if(removeFile) { file.remove(bamFile) }
  
  return(hits)
}

#' Pair up alignments in a GRanges object
#'
#' Given a GRanges object, the function uses specified gaplength parameter to 
#' pair up reads where the qName column ends with "atpersand pairname atpersand" 
#' which is outputted by \code{\link{extractSeqs}}.
#'
#' @param psl.rd a GRanges object with qNames ending in "atpersand pairname atpersand".
#' @param maxGapLength maximum gap allowed between end of pair1 and start of 
#' pair2. Default is 2500 bp.
#' @param sameStrand should pairs be aligned to the same strand or in same 
#' orientationin the reference genome? Default is TRUE. This is 'TRUE' because 
#' pair2 reads are reverseComplemented when reading in data in 
#' \code{\link{findBarcodes}}
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'
#' @return a GRanges object with reads paired up denoted by "paired" column. 
#' Improper pairs or unpaired reads are returned with "paired" column as FALSE.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{blatSeqs}}, 
#' \code{\link{read.blast8}}, \code{\link{read.psl}}, 
#' \code{\link{getIntegrationSites}}, \code{\link{read.BAMasPSL}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' psl.rd <- read.BAMasPSL(bamFile=c("sample1hits.bam","sample2hits.bam"))
#' pairUpAlignments(psl.rd)
#' }
pairUpAlignments <- function(psl.rd=NULL, maxGapLength=2500, 
                             sameStrand=TRUE, parallel=TRUE) {
  dp <- NULL
  
  .checkArgsSetDefaults_ALIGNed()
  stopifnot("qName" %in% names(mcols(psl.rd)))
  qName <- NULL
  
  ### identify pairs ###
  mcols(psl.rd)$pair <- sub(".+@(.+)@","\\1", mcols(psl.rd)$qName)
  mcols(psl.rd)$qName <- sub("(.+)@.+@","\\1", mcols(psl.rd)$qName)
  stopifnot(!any(is.na(mcols(psl.rd)$pair)))
  
  ### check pairs ###
  p1 <- mcols(psl.rd)$qName[mcols(psl.rd)$pair=="pair1"]
  p2 <- mcols(psl.rd)$qName[mcols(psl.rd)$pair=="pair2"]
  loners <- unique(c(setdiff(p1, p2), setdiff(p2, p1)))
  loners <- subset(psl.rd, mcols(psl.rd)$qName %in% loners)
  mcols(loners)$paired <- FALSE
  psl.rd <- subset(psl.rd, !mcols(psl.rd)$qName %in% mcols(loners)$qName)  
  rm("p1","p2")
  
  ### pair up pairs ###
  # collapse pairs first, then tag successfully merged pairs by checking counts
  # add respective metadata back to successfully collapsed reads
  # for unsuccessful ones, match which pairs overlap with the collapsed area then...
  # subset those reads and see they're multihits by any chance...
  # if not tag then as paired=FALSE
  
  bore <- GRanges(seqnames=paste(as.character(seqnames(psl.rd)), 
                                 psl.rd$qName, sep="@tempy@"), 
                  ranges=ranges(psl.rd), strand=strand(psl.rd))
  reduced <- reduce(bore, min.gapwidth=maxGapLength, 
                    ignore.strand=(!sameStrand)) 
  rm(bore)
  
  reduced <- GRanges(seqnames=sub("(.+)@tempy@.+","\\1",
                                  as.character(seqnames(reduced))), 
                     ranges=ranges(reduced), strand=strand(reduced),
                     qName=sub(".+@tempy@(.+)","\\1",
                               as.character(seqnames(reduced))))
  counts <- table(mcols(reduced)$qName)
  
  proper <- reduced[mcols(reduced)$qName %in% names(counts[counts==1]),]
  bore <- subset(psl.rd, 
                 mcols(psl.rd)$qName %in% names(counts[counts==1]) &
                   mcols(psl.rd)$pair=='pair1')
  rows <- match(mcols(proper)$qName, mcols(bore)$qName)
  stopifnot(identical(mcols(proper)$qName, mcols(bore)$qName[rows]))
  mcols(proper) <- mcols(bore[rows])
  mcols(proper)$paired <- TRUE
  
  reduced <- reduced[mcols(reduced)$qName %in% names(counts[counts!=1]),]
  psl.rd <- subset(psl.rd, 
                   mcols(psl.rd)$qName %in% names(counts[counts!=1]))    
  rm(counts)
  
  reduced <- split(reduced, mcols(reduced)$qName)
  psl.rd <- split(psl.rd, mcols(psl.rd)$qName)
  psl.rd <- psl.rd[names(reduced)]
  stopifnot(identical(names(reduced),names(psl.rd)))
  
  pair <- bpmapply(function(x,y) {
    res <- as.data.frame(findOverlaps(x,y))
    res$pair <- mcols(y)$pair[res$subjectHits]
    res$qName <- mcols(y)$qName[res$subjectHits]
    test <- aggregate(pair~qName+queryHits, data=res, 
                      FUN=function(x) paste(sort(x), collapse=""))
    
    # remove cases where a pair gives rise to many reduced 
    # regions but one may ... 
    # contain one pair and other contains both
    hasBoth <- subset(test, pair=="pair1pair2")
    test <- subset(test, pair!="pair1pair2" & !qName %in% hasBoth$qName)
    
    # despite the collapse if each pair is yielding it's own hit 
    # then chances are...
    # you're on different chromosome, too far apart, or 
    # different strands!
    if(nrow(test)>0) {
      unpaired <- y[subset(res, queryHits %in% test$queryHits)$subjectHits]
      mcols(unpaired)$paired <- FALSE
    }
    
    # reduced is >1 but each collapse region has both reads...
    # hence a multihit #            
    res <- subset(res, queryHits %in% hasBoth$queryHits)
    if(nrow(res)>0) {
      paired <- x[unique(res$queryHits)]
      mcols(paired) <- mcols(y[res[res$pair=="pair1","subjectHits"]])
      mcols(paired)$paired <- TRUE
    }
    
    rm("res","test","hasBoth")
    
    toReturn <- c()
    if(exists("unpaired")) {
      toReturn <- c(toReturn, unpaired)
    } 
    if(exists("paired")) {
      toReturn <- c(toReturn, paired)
    }  
    do.call(c, toReturn)
  }, reduced, psl.rd, SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)
  
  pair <- suppressWarnings(do.call(c, pair))
  
  c(proper,pair,loners)
}

#' Start/Stop a gfServer instance
#'
#' Start or Stop a gfServer with indexed reference genome to align batch of 
#' sequences using BLAT gfServer/gfClient protocol.
#'
#' @param seqDir absolute or relative path to the genome index (nib/2bit files).
#' @param host name of the machine to run gfServer on. Default: localhost
#' @param port a port number to host the gfServer with. Default is 5560.
#' @param gfServerOpts a character vector of options to be passed to gfServer 
#' command on top of server defaults. Default: c(repMatch=112312, stepSize=5, 
#' tileSize=10, maxDnaHits=10). Set this to NULL to start gfServer with defaults.
#'
#' @return system command status for executing gfServer command.
#'
#' @seealso \code{\link{stopgfServer}}, \code{\link{read.psl}}, 
#' \code{\link{blatSeqs}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' #startgfServer(seqDir="/usr/local/blatSuite34/hg18.2bit",port=5560)
#' #stopgfServer(port=5560)
startgfServer <- function(seqDir=NULL, host="localhost", port=5560, 
                          gfServerOpts=c(repMatch=112312, stepSize=5, 
                                         tileSize=10, maxDnaHits=10)) {
  
  if(length(system("which gfServer",intern = TRUE))==0) { 
    stop("Command gfServer for BLAT not found!")
  }
  
  if(is.null(seqDir)) {
    stop("Please define the path of nib/2bit files containing the ",
         "indexed reference sequence(s)")
  }
  
  cmd <- sprintf("gfServer start %s %i %s %s &", 
                 host, port, 
                 ifelse(!is.null(gfServerOpts),
                        paste(paste("-",names(gfServerOpts),sep=""), 
                              gfServerOpts, collapse=" ", sep="="),
                        ""), 
                 normalizePath(seqDir, mustWork = TRUE))
  message(cmd)
  system(cmd)        
  
  ## wait for server to load & be ready
  message("Loading BLAT server...please wait.")    
  searchCMD <- sprintf("gfServer status %s %s", host, port)
  while(system(searchCMD,ignore.stderr=TRUE)!=0) {
    cat(".")
    Sys.sleep(10)
  }
}

#' @rdname startgfServer
stopgfServer <- function(host="localhost", port=NULL) {
  if(is.null(port)) {
    stop("Please define the port gfServer is running on.")
  }
  
  cmd <- sprintf("kill `ps ax | grep '%s' | grep -v 'grep' | awk '{print $1}'`", 
                 paste("gfServer start", host, port))
  system(cmd)
}

#' Align a listed DNAStringSet to a reference using gfClient or standalone BLAT.
#'
#' Align sequences from a listed DNAStringSet object returned from 
#' \code{\link{extractSeqs}} to an indexed reference genome using 
#' gfServer/gfClient protocol or using standalone BLAT and return the psl file 
#' as a GRanges object. This function heavily relies on defaults of 
#' \code{\link{blatSeqs}}.
#'
#' @param dnaSetList DNAStringSet object containing sequences to be aligned against the reference.
#' @param ... parameters to be passed to \code{\link{blatSeqs}}.
#'
#' @return a list of GRanges object reflecting psl file type per set of sequences.
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#'  \code{\link{startgfServer}}, \code{\link{stopgfServer}}, 
#'  \code{\link{blatSeqs}}, \code{\link{read.psl}}, 
#'  \code{\link{pslToRangedObject}}, \code{\link{read.blast8}}
#'
#' @export
#' 
blatListedSet <- function(dnaSetList=NULL, ...) {
  if(is.null(dnaSetList)) {
    stop("dnaSetList is empty. Please supply a listed DNAStringSet ",
         "returned from extractSeqs() to be aligned against a reference")
  }
  
  sapply(names(dnaSetList), function(x) {
    do.call(GRangesList, sapply(names(dnaSetList[[x]]), function(y) {
      if(length(dnaSetList[[x]][[y]])>0) {
        message("BLATing ",y)
        outFiles <- blatSeqs(query=dnaSetList[[x]][[y]], ...)
        read.psl(outFiles, bestScoring=TRUE, asGRanges=TRUE, removeFile=TRUE, 
                 parallel=FALSE)
      }
    }))
  })
}

#' Convert psl dataframe to GRanges
#'
#' Convert psl dataframe to GRanges object using either the query or target as 
#' the reference data column. 
#'
#' @param x dataframe reflecting psl format
#' @param useTargetAsRef use target(tName) or query(qName) as the chromosome or 
#' the reference data. Default is TRUE.
#' @param isblast8 the input dataframe blast8 format output from BLAT. 
#' Default is FALSE.
#'
#' @return a GRanges object reflecting psl file type.
#'
#' @seealso \code{\link{read.psl}}, \code{\link{read.blast8}}, 
#' \code{\link{blatListedSet}}
#'
#' @export
#'
#' @examples
#' data(psl)
#' psl <- head(psl)
#' pslToRangedObject(psl)
#' pslToRangedObject(psl, useTargetAsRef=FALSE)
pslToRangedObject <- function(x, useTargetAsRef=TRUE, isblast8=FALSE) {
  if(useTargetAsRef) {
    metadataCols <- c(setdiff(names(x), c("tName","tStart","tEnd","strand")),
                      ifelse(isblast8, NA, "tStarts"))
    out <- GRanges(seqnames=x$tName, IRanges(start=x$tStart, end=x$tEnd),
                   strand=x$strand)     
  } else {
    metadataCols <- c(setdiff(names(x), c("qName","qStart","qEnd","strand")),
                      ifelse(isblast8, NA, "qStarts"))
    out <- GRanges(seqnames=x$qName, IRanges(start=x$qStart, end=x$qEnd),
                   strand=x$strand)
  }
  
  for(f in na.omit(metadataCols)) {
    mcols(out)[[f]] <- x[,f]
  } 
  
  out
}

#' Split DNA sequences into smaller files.
#'
#' Given a vector of sequences or DNAStringSet or a FASTA filename, the function 
#' splits it into smaller pieces as denoted by totalFiles parameter.
#'
#' @param x a DNAStringSet object, or a FASTA filename.
#' @param totalFiles an integer indicating how many files to create. Default is 4.
#' @param suffix a word to add to each file created. Default is "tempy".
#' @param filename name of the file if x is a DNAStringSet object. 
#' Default is "queryFile.fa".
#'
#' @return a vector of filename names created.
#'
#' @seealso \code{\link{blatSeqs}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' seqs <- DNAStringSet(sapply(sample(c(100:1000), 500), 
#' function(size) paste(sample(DNA_BASES, size, replace=TRUE), collapse=""))) 
#' splitSeqsToFiles(seqs,5,"tempyQ","myDNAseqs.fa")
#' }
splitSeqsToFiles <- function(x, totalFiles=4, suffix="tempy", 
                             filename="queryFile.fa") {
  if(is.atomic(x)) {
    message("Splitting file ",x)
    totalSeqs <- length(fasta.info(x, use.names=FALSE))
    chunks <- round(totalSeqs/totalFiles)
    ## incase totalSeqs is lower than number of files to be created!
    chunks <- ifelse(chunks>0, chunks, totalSeqs) 
    
    starts <- seq(0, totalSeqs, by=chunks) ## create chunks of starts    
    for(skippy in starts[starts!=totalSeqs]) {
      filename.out <- paste(x, skippy, runif(1), suffix, sep=".")
      ## no need to read the entire file...save memory by reading in N lines
      query.tmp <- readBStringSet(x,nrec=chunks, skip=skippy) 
      writeXStringSet(query.tmp, filepath=filename.out, format="fasta")            
    }
    return(list.files(path=dirname(x), 
                      pattern=paste0(basename(x),".*", suffix, "$"), 
                      full.names=TRUE))
  } else if (class(x)=="DNAStringSet") {
    message("Splitting Reads.")
    totalSeqs <- length(x)
    chunks <- round(totalSeqs/totalFiles)
    starts <- seq(1, totalSeqs, by=chunks)
    stops <- unique(c(seq(chunks, totalSeqs, by=chunks), totalSeqs))
    stopifnot(length(starts)==length(stops))        
    for(skippy in 1:length(starts)) {
      filename.out <- paste(filename, skippy, runif(1), suffix, sep=".")            
      writeXStringSet(x[starts[skippy]:stops[skippy]], filepath=filename.out,
                      format="fasta")            
    }            
    return(list.files(path=".", 
                      pattern=paste(filename,".*",suffix,"$",sep=""), 
                      full.names=TRUE))
  } else {
    stop("Dont know what is supplied in parameter x.")
  }
}

#' Align sequences using BLAT.
#'
#' Align batch of sequences using standalone BLAT or gfServer/gfClient protocol 
#' against an indexed reference genome. Depending on parameters provided, the 
#' function either aligns batch of files to a reference genome using gfClient or
#' takes sequences from query & subject parameters and aligns them using 
#' standalone BLAT. If standaloneBlat=FALSE and gfServer is not launched 
#' apriori, this function will start one using \code{\link{startgfServer}} 
#' and kill it using \code{\link{stopgfServer}} upon successful execution. 
#'
#' @param query an object of DNAStringSet, a character vector of filename(s), 
#' or a path/pattern of fasta files to BLAT. Default is NULL.
#' @param subject an object of DNAStringSet, a character vector, or a path to 
#' an indexed genome (nibs,2bits) to serve as a reference or target to the query.
#' Default is NULL. If the subject is a path to a nib or 2bit file, then 
#' standaloneBlat will not work!
#' @param standaloneBlat use standalone BLAT as suppose to gfServer/gfClient 
#' protocol. Default is TRUE.
#' @param port the same number you started the gfServer with. Required if 
#' standaloneBlat=FALSE. Default is 5560.
#' @param host name of the machine running gfServer. Default is 'localhost' and 
#' only used when standaloneBlat=FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is registered, 
#' then a serial version is ran using \code{\link{SerialParam}}.
#' @param numServers launch >1 gfServer and load balance jobs? This only
#' applies when parallel=TRUE and standaloneBlat=FALSE. Enable this option only
#' if the machine has a lot of RAM! Option ignored if launched gfServer is found
#' at specified host and port. Default is 1. 
#' @param gzipResults gzip the output files? Default is TRUE.
#' @param blatParameters a character vector of options to be passed to 
#' gfClient/BLAT command except for 'nohead' option. Default: c(minIdentity=90,
#' minScore=10, stepSize=5, tileSize=10, repMatch=112312, dots=50, maxDnaHits=10, 
#' q="dna", t="dna", out="psl"). Be sure to only pass parameters accepted by 
#' either BLAT or gfClient. For example, if repMatch or stepSize parameters are 
#' specified when using gfClient, then the function will simply ignore them! 
#' The defaults are configured to align a 19bp sequence with 90\% identity.
#'
#' @return a character vector of psl filenames. Each file provided is split by 
#' number of parallel workers and with read number denoting the cut. Files are 
#' cut in smaller pieces to for the ease of read & write into a single R session. 
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, 
#' \code{\link{startgfServer}}, \code{\link{stopgfServer}}, \code{\link{read.psl}},
#' \code{\link{splitSeqsToFiles}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples
#'  
#' \donttest{
#' blatSeqs(dnaSeqs, subjectSeqs, blatParameters=c(minIdentity=90, minScore=10, 
#' tileSize=10, dots=10, q="dna", t="dna", out="blast8"))
#' blatSeqs(dnaSeqs, "/usr/local/genomeIndex/hg18.2bit", standaloneBlat=FALSE)
#' blatSeqs("mySeqs.fa", "/usr/local/genomeIndex/hg18.2bit", standaloneBlat=FALSE)
#' blatSeqs("my.*.fa", "/usr/local/genomeIndex/hg18.2bit", standaloneBlat=FALSE)
#' }
blatSeqs <- function(query=NULL, subject=NULL, standaloneBlat=TRUE, port=5560, 
                     host="localhost", parallel=TRUE, numServers=1L,
                     gzipResults=TRUE,
                     blatParameters=c(minIdentity=90, minScore=10, stepSize=5, 
                                      tileSize=10, repMatch=112312, dots=50, 
                                      maxDnaHits=10, q="dna", t="dna", 
                                      out="psl")) {
  
  if(length(system("which blat",intern = TRUE))==0) { 
    stop("Command blat not found!")
  }
  
  ## get all BLAT options from the system for comparison to blatParameters later
  suppressWarnings(blatOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
                                          grep("\\s+-.+=", 
                                               system("blat",intern=TRUE), 
                                               value=TRUE))))
  
  ##suppressWarnings(gfClientOpts <- unique(sub("\\s+-(\\S+)=\\S+.+", "\\1", 
  ##                                            grep("\\s+-.+=", 
  ##                                                 system("gfClient", 
  ##                                                        intern=TRUE), 
  ##                                                 value=TRUE))))

  gfClientOpts <- c("dots", "minScore", "minIdentity", "out", "maxIntron")

  gfServerOpts <- c("tileSize","stepSize","minMatch","maxGap", "trans","log",
                    "seqLog","syslog","logFacility","mask","repMatch",
                    "maxDnaHits", "maxTransHits","maxNtSize","maxAsSize",
                    "canStop")                                               
  
  if(!standaloneBlat) { 
    message("Using gfClient protocol to perform BLAT.")
    if(is.null(port)) {
      stop("The port paramter is empty. ",
           "Please define the port used to start gfServer with")
    }
  }
  
  ## check the subject parameter
  if(is.null(subject) | length(subject)==0) {
    stop("The subject parameter is empty. ", 
         "Please supply subject sequences or a path to 2bit or nib files to ",
         "serve as reference/target")
  } else {
    subjectFile <- NULL
    if(is.atomic(subject)) {
      if (any(grepl("\\.2bit$|\\.nib$", subject, ignore.case=TRUE))) {
        if(standaloneBlat) { 
          stop("Standalone BLAT cannot be used when subject is an indexed ",
               "nib or 2bit file.") 
        }
        indexFileDir <- dirname(subject)
        subjectFile <- list.files(path=indexFileDir, 
                                  pattern=basename(subject), full.names=TRUE)
        if(length(subjectFile)==0) { 
          stop("The file(s) supplied in subject parameter doesn't exist.") }
      } else {
        ## change object type if necessary for troubleshooting purpose in later 
        ## steps
        subject <- DNAStringSet(subject)
      }
    }
    
    if(is.null(subjectFile)) {
      ## subjectFile is still null so it means that subject is a DNAStringSet
      if(is.null(names(subject))) { ## add names of subject if not present
        names(subject) <- paste("subject", 1:length(subject))
      }
      
      ## write out the subject sequences into a fasta file
      filename.seq <- paste("subjectFile.fa",runif(1),"tempyS",sep=".")      
      writeXStringSet(subject, filepath=filename.seq, format="fasta")                                  
      subjectFile <- filename.seq
    }
  }
  
  ## check the query parameter
  if(is.null(query) | length(query)==0) {
    stop("The query parameter is empty. Please supply reads to be aligned")
  } else {
    queryFiles <- NULL
    if(is.atomic(query)) {
      if (any(grepl("\\.fna$|\\.fa$|\\.fastq$|\\.fasta$|\\*", query, 
                    ignore.case=TRUE))) {
        ## detect whether query paramter is a regex or list of files
        if(any(grepl("\\*|\\$|\\+|\\^",query))) {
          queryFiles <- list.files(path=dirname(query), pattern=basename(query), 
                                   full.names=TRUE)            
        } else {
          queryFiles <- query
        }
        
        if(parallel) {
          ## split the fasta files into smaller chunks for parallel BLATing
          queryFiles <- unlist(sapply(queryFiles,
                                      function(f) 
                                        splitSeqsToFiles(f, bpworkers(),
                                                         "tempyQ")), 
                               use.names=FALSE)                    
        }
      } else {
        ## change object type if necessary for troubleshooting purpose in later steps
        query <- DNAStringSet(query)
      }
    }
    
    if(is.null(queryFiles)) {
      ## queryFiles is still null so it means that query is a DNAStringSet           
      if(is.null(names(query))) {  ## fix names of query if not present
        names(query) <- paste("read", 1:length(query),sep="-")
      }  
      
      ## write out the query sequences into fasta files
      if(parallel) {
        queryFiles <- splitSeqsToFiles(query, bpworkers(), "tempyQ")
      } else {
        queryFiles <- paste("queryFile.fa",runif(1),"tempyQ",sep=".")          
        writeXStringSet(query, filepath=queryFiles, format="fasta")                
      }
    }
  }
  
  ## perform the Blatting of queryFiles vs subjectFile/indexFiles using 
  ## gfClient/standalone BLAT  
  
  ## do some formatting ##
  queryFiles <- as.character(queryFiles)
  subjectFile <- as.character(subjectFile)
    
  dp <- if(parallel){ bpparam() } else { SerialParam() }
  
  ## BLAT it ##
  if(standaloneBlat) {        
    blatOpts <- blatParameters[names(blatParameters) %in% blatOpts]
    stopifnot(length(subjectFile)==1)
    filenames <- bplapply(queryFiles, function(x) {
      filename.out <- paste(x, blatOpts["out"], sep=".")
      cmd <- paste("blat", paste(paste0("-",names(blatOpts)), blatOpts, 
                                 collapse=" ", sep="="), "-noHead", 
                   subjectFile, x, filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$",x)) { system(sprintf("rm %s",x)) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out
    }, BPPARAM=dp)
    
    if(grepl("\\.tempyS$",subjectFile)) { system(sprintf("rm %s",subjectFile)) }
    
  } else {
    # start the gfServer if not started already! #
    killFlag <- FALSE
    port <- port + 0:(numServers-1)
    for(n in 1:numServers) {
      searchCMD <- sprintf("gfServer status %s %s", host, port[n])
      if(system(searchCMD,ignore.stderr=TRUE)!=0) {
        message(sprintf("Starting gfServer # %s.", n))
        startgfServer(seqDir=subjectFile, host=host, port=port[n], 
                      gfServerOpts=blatParameters[names(blatParameters) 
                                                  %in% gfServerOpts])
        killFlag <- TRUE
      }
    }
    
    gfClientOpts <- blatParameters[names(blatParameters) %in% gfClientOpts]
    stopifnot(length(subjectFile)>0)
    filenames <- bpmapply(function(port, x) {
      filename.out <- paste(x, gfClientOpts["out"], sep=".")
      cmd <- paste("gfClient",
                   paste(paste0("-", names(gfClientOpts)), gfClientOpts, 
                         collapse=" ", sep="="), "-nohead", host, port, "/", x, 
                   filename.out)
      message(cmd)
      system(cmd)
      
      ## no need to save splitted files!
      if(grepl("\\.tempyQ$", x)) { file.remove(x) } 
      
      if(gzipResults) { 
        system(paste("gzip", filename.out))
        filename.out <- paste(filename.out, "gz", sep=".") 
      }
      filename.out      
    }, rep(port, length=length(queryFiles)), queryFiles, 
    SIMPLIFY=FALSE, USE.NAMES=FALSE, BPPARAM=dp)  
    
    ## only kill if the gfServer was started from within this function ## 
    if(killFlag) {
      # stop to conserve memory #
      message("Kill gfServer.")        
      sapply(port, function(x) stopgfServer(port=x))
    }
  }
  return(unlist(filenames))
}

#' Return PSL file columns with classes
#'
#' Print out required fields & classes of PSL file format
#'
#' @param withClass return classes for each column.
#' @return vector of PSL column names
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}}, 
#' \code{\link{startgfServer}}, \code{\link{blatSeqs}}, \code{\link{read.blast8}}, 
#' \code{\link{read.BAMasPSL}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' pslCols()
#' 
pslCols <- function(withClass=TRUE) {
  cols <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", 
            "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
            "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", 
            "blockCount", "blockSizes", "qStarts", "tStarts")
  
  cols.class <- c(rep("numeric",8), rep("character",2), rep("numeric",3),
                  "character", rep("numeric",4), rep("character",3))
  
  if(withClass) {
    structure(cols.class, names=cols)
  } else {
    cols
  }
}

#' Read PSL file(s) outputted by BLAT
#'
#' Given filename(s), the function reads the PSL file format from BLAT as a 
#' data frame and performs basic score filtering if indicated. Any other file 
#' format will yield errors or erroneous results. Make sure there is no 
#' header row! See required columns in \code{\link{pslCols}}.
#'
#' @param pslFile PSL filename, or vector of filenames, or a pattern of files 
#' to import.
#' @param bestScoring report only best scoring hits instead of all hits. 
#' Default is TRUE. Score is calculated by 
#' matches-misMatches-qBaseInsert-tBaseInsert.
#' @param asGRanges return a GRanges object instead of a dataframe.
#'  Default is FALSE
#' @param removeFile remove the PSL file(s) after importing. 
#' Default is FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'
#' @return a dataframe reflecting psl file type. If asGRanges=TRUE, 
#' then a GRanges object.
#'
#' @note If parallel=TRUE, then be sure to have a parallel backend registered 
#' before running the function. One can use any of the following 
#' \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#' \code{\link{startgfServer}}, \code{\link{blatSeqs}}, 
#' \code{\link{read.blast8}}, \code{\link{read.BAMasPSL}},
#' \code{\link{pslToRangedObject}}, \code{\link{write.psl}}
#'
#' @export
#'
#' @examples 
#' data(psl)
#' pslFile <- tempfile()
#' write.psl(psl, filename = pslFile)
#' head(read.psl(pslFile=pslFile))
#' \dontrun{
#' # read many PSL files matching the regex #
#' psl <- read.psl(pslFile="processed.*.psl$")
#' }
read.psl <- function(pslFile=NULL, bestScoring=TRUE, asGRanges=FALSE, 
                     removeFile=TRUE, parallel=FALSE) {
  qName <- dp <- NULL
  files <- pslFile
  .checkArgsSetDefaults_ALIGNed()
    
  ## setup psl columns + classes
  cols <- pslCols()
  
  hits <- bplapply(files, function(x) {
    message(x)
    ## add extra fields incase pslx format ##
    ncol <- max(count.fields(x, sep = "\t"))
    if(ncol > length(cols)) {
      for(f in 1:(ncol-length(cols))) {
        cols[paste0("V",f)] <- "character"
      }
    }
    hits.temp <- read.delim(x, header=FALSE, col.names=names(cols), 
                            stringsAsFactors=FALSE, colClasses=cols)    
    if(bestScoring) {  
      ## do round one of bestScore here to reduce file size          
      hits.temp$score <- with(hits.temp, 
                              matches-misMatches-qBaseInsert-tBaseInsert)
      isBest <- with(hits.temp, ave(score, qName, FUN=function(x) x==max(x)))
      hits.temp <- hits.temp[as.logical(isBest),]
      rm("isBest")
    }
    hits.temp    
  }, BPPARAM=dp)  
  hits <- unique(rbind.fill(hits))
  
  if(nrow(hits)==0) {
    if(removeFile) { file.remove(pslFile) }
    stop("No hits found")
  }
  
  ## do round two of bestScore incase any got missed in round one
  if(bestScoring) {
    message("\t cherry picking!")
    hits$score <- with(hits, matches-misMatches-qBaseInsert-tBaseInsert)    
    isBest <- with(hits, ave(score, qName, FUN=function(x) x==max(x)))
    hits <- hits[as.logical(isBest),]
    rm("isBest")
  }
  
  if(asGRanges) {
    hits <- pslToRangedObject(hits, useTargetAsRef=TRUE)
  }
  
  message("2. Ordering by qName")
  if(is(hits,"GRanges")) {
    hits <- sort(hits, by=~qName)
  } else {
    hits <- arrange(hits, qName)
  }  
  
  if(removeFile) { file.remove(pslFile) }
  
  return(hits)
}

#' Write PSL file from dataframe or GRanges
#'
#' Given a data frame or GRanges object, the function write a tab deliminated PSL file
#'
#' @param x data frame or GRanges object with required columns for psl file format.
#' @param filename name for the output PSL file. Default is "out.psl"
#' @param header include PSL header line. Default is FALSE.
#' @param includeOtherCols nclude other non PSL specific columns from x in 
#' the output. Default is FALSE.
#'
#' @return name of the output PSL file
#'
#' @seealso \code{\link{read.psl}}, \code{\link{blatSeqs}}, 
#' \code{\link{read.blast8}}, \code{\link{read.BAMasPSL}}, 
#' \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' data(psl)
#' pslFile <- tempfile()
#' write.psl(psl, filename = pslFile)
write.psl <- function(x, filename="out.psl", header=FALSE, 
                      includeOtherCols=FALSE) {
  if(is.null(x) | length(x)==0) {
    stop("x parameter is empty or of length 0.")
  }
  
  if(!is(x, "data.frame")) {
    isGRanges <- is(x, "GRanges")
    x <- BiocGenerics::as.data.frame(x)
    if(isGRanges) {
      names(x)[1:3] <- c("tName","tStart","tEnd")
    }
    x$width <- NULL
  }
  
  ## PSL columns ##
  cols <- pslCols(withClass=FALSE)
  
  ## missing columns ##
  missing <- setdiff(cols, names(x))   
  if(length(missing)>0) {
    stop("Following columns missing from the input: ", 
         paste(missing,collapse=", "))
  }
  
  if(includeOtherCols) {
    cols <- union(cols, names(x))  
  }
  
  write.table(x[,cols], file=filename, sep="\t", row.names=FALSE, quote=FALSE, 
              col.names=header)
  
  return(filename)
}

#' Read blast8 file(s) outputted by BLAT
#'
#' Given filename(s), the function reads the blast8 file format from BLAT as a 
#' data frame and performs basic score filtering if indicated. Any other file 
#' format will yield errors or erroneous results.
#'
#' @param files blast8 filename, or vector of filenames, or a pattern of files 
#' to import.
#' @param asGRanges return a GRanges object instead of a dataframe.
#'  Default is TRUE Saves memory!
#' @param removeFile remove the blast8 file(s) after importing. 
#' Default is FALSE.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#'
#' @return a dataframe or GRanges object reflecting blast8 file type.
#'
#' @note If parallel=TRUE, then be sure to have a parallel backend registered 
#' before running the function. One can use any of the following 
#' \code{\link{MulticoreParam}} \code{\link{SnowParam}}
#'
#' @seealso \code{\link{pairwiseAlignSeqs}}, \code{\link{vpairwiseAlignSeqs}},
#' \code{\link{startgfServer}}, \code{\link{blatSeqs}}, \code{\link{read.psl}}
#'
#' @export
#'
#' @examples 
#' # this function works similar to read.psl #
#' #read.blast8(files="processed.*.blast8$")
#' #read.blast8(files=c("sample1hits.blast8","sample2hits.blast8"))
#'
read.blast8 <- function(files=NULL, asGRanges=FALSE,
                        removeFile=TRUE, parallel=FALSE) {
  qName <- dp <- NULL

  .checkArgsSetDefaults_ALIGNed()

  ## setup blast8 columns + classes
  cols <- c("qName", "tName", "identity", "span", "misMatches", "gaps", 
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")
  cols.class <- c(rep("character",2), rep("numeric",10))
  cols <- structure(cols.class, names=cols)
  
  hits <- bplapply(files, function(x) {
    message(x)
    ## add extra fields incase pslx format ##
    ncol <- max(count.fields(x, sep = "\t"))
    if(ncol > length(cols)) {
      for(f in 1:(ncol-length(cols))) {
        cols[paste0("V",f)] <- "character"
      }
    }
    hits.temp <- read.delim(x, header=FALSE, col.names=names(cols), 
                            stringsAsFactors=FALSE, colClasses=cols)
    hits.temp$strand <- with(hits.temp, ifelse(tStart>tEnd, "-","+"))
    
    # switch tStart & tEnd for cases where strand=='-' 
    # since it's reversed in blast8 format.
    rows <- hits.temp$strand=='-'
    tstarts <- hits.temp$tEnd[rows]
    tends <- hits.temp$tStart[rows]
    hits.temp$tStart[rows] <- tstarts
    hits.temp$tEnd[rows] <- tends
    rm("tstarts","tends","rows")
    hits.temp
  }, BPPARAM=dp)
  hits <- unique(rbind.fill(hits))
  
  if(nrow(hits)==0) {
    if(removeFile) { file.remove(files) }
    stop("No hits found")
  }
  
  if(asGRanges) {
    hits <- pslToRangedObject(hits, useTargetAsRef=TRUE, isblast8=TRUE)
  }
  
  message("3. Ordering by qName")
  if(is(hits,"GRanges")) {
    hits <- sort(hits, by=~qName)
  } else {
    hits <- arrange(hits, qName)
  }
  
  if(removeFile) { file.remove(files) }
  
  return(hits)
}

#' Obtain integration sites from BLAT output
#'
#' Given a GRanges object from \code{\link{read.psl}}, the function uses 
#' specified filtering parameters to obtain integration sites and maintain 
#' sequence attrition. The function will remove any non-best scoring alignments
#' from the object if not already filtered apriori.
#'
#' @param psl.rd a GRanges object reflecting psl format where tName is 
#' the seqnames.
#' @param startWithin upper bound limit of where the alignment should start 
#' within the query. Default is 3.
#' @param alignRatioThreshold cuttoff for (alignment span/read length). 
#' Default is 0.7.
#' @param genomicPercentIdentity cuttoff for (1-(misMatches/matches)). 
#' Default is 0.98.
#' @param correctByqStart use qStart to correct genomic position. This would
#' account for sequencing/trimming errors. 
#' Position=ifelse(strand=="+",tStart-qStart,tEnd+qStart). Default is TRUE.
#' @param oneBased the coordinates in psl files are "zero based half open". 
#' The first base in a sequence is numbered zero rather than one. Enabling 
#' this would add +1 to the start and leave the end as is. Default is FALSE.
#'
#' @return a GRanges object with integration sites which passed all filtering
#' criteria. Each filtering parameter creates a new column to flag if a 
#' sequence/read passed that filter which follows the scheme: 
#' 'pass.FilterName'. Integration Site is marked by new column named 'Position'.
#'
#' @seealso \code{\link{startgfServer}}, \code{\link{read.psl}}, 
#' \code{\link{blatSeqs}}, \code{\link{blatListedSet}}, 
#' \code{\link{findIntegrations}}, \code{\link{pslToRangedObject}}, 
#' \code{\link{clusterSites}}, \code{\link{isuSites}}, 
#' \code{\link{crossOverCheck}}, \code{\link{read.blast8}}
#'
#' @export
#'
#' @examples 
#' data(psl)
#' psl.rd <- pslToRangedObject(psl)
#' getIntegrationSites(psl.rd)
getIntegrationSites <- function(psl.rd=NULL, startWithin=3, 
                                alignRatioThreshold=0.7, 
                                genomicPercentIdentity=0.98, 
                                correctByqStart=TRUE, oneBased=FALSE) {
  stopifnot((is(psl.rd,"GRanges") | is(psl.rd,"GAlignments")) & 
              !is.null(psl.rd) & !is.null(startWithin) & length(psl.rd)!=0 &
              !is.null(alignRatioThreshold) & !is.null(genomicPercentIdentity))
  
  ## check if required columns exist ##
  absentCols <- setdiff(c("qName", "qStart", "qSize","matches", "misMatches", 
                          "qBaseInsert", "tBaseInsert"), 
                        colnames(mcols(psl.rd)))
  
  if(length(absentCols)>0) {
    stop("Following columns are absent from psl.rd object: ",
         paste(absentCols,collapse=","))
  }
  
  ## get the integration position by correcting for any insertions due to 
  # sequencing errors
  if(correctByqStart) {
    mcols(psl.rd)$Position <- ifelse(as.character(strand(psl.rd))=="+",
                                     start(psl.rd)-mcols(psl.rd)$qStart,
                                     end(psl.rd)+mcols(psl.rd)$qStart)
  } else {
    mcols(psl.rd)$Position <- start(flank(psl.rd, width=-1))
  }
  
  ## get +1 based coordinate ##    
  if(oneBased) {
    mcols(psl.rd)$Position <- ifelse(as.character(strand(psl.rd))=="+",
                                     mcols(psl.rd)$Position+1,
                                     mcols(psl.rd)$Position)
  }
  
  # get scores for picking best hits and identify multihits later
  # check if scoring filtering hasn't already been applied by blat functions
  if(!"score" %in% colnames(mcols(psl.rd))) {
    message("Adding score column.")
    mcols(psl.rd)$score <- with(as.data.frame(mcols(psl.rd)), 
                                matches-misMatches-qBaseInsert-tBaseInsert)     
    isBest <- ave(mcols(psl.rd)$score, as.character(mcols(psl.rd)$qName), 
                  FUN=function(x) x==max(x))    
    psl.rd <- psl.rd[as.logical(isBest),]
    rm("isBest")
    cleanit <- gc()
  }    
  
  message("Performing QC checks.")
  # remove rows where the best hit dont start within first X bp
  mcols(psl.rd)$pass.startWithin <- mcols(psl.rd)$qStart<=startWithin
  
  # check if aligned ratio matches the threshold
  mcols(psl.rd)$alignRatio <- mcols(psl.rd)$score/mcols(psl.rd)$qSize
  mcols(psl.rd)$pass.alignRatio <- mcols(psl.rd)$alignRatio >= alignRatioThreshold
  
  # check for %identity    
  mcols(psl.rd)$percIdentity <- 1-(mcols(psl.rd)$misMatches/mcols(psl.rd)$matches)
  mcols(psl.rd)$pass.percIdentity <- 
    mcols(psl.rd)$percIdentity >= genomicPercentIdentity
  
  ## find which query aligned to multiple places with equally good score aka...
  ## multihits
  cloneHits <- table(mcols(psl.rd)$qName)
  cloneHits <- cloneHits[as.character(mcols(psl.rd)$qName)]>1
  mcols(psl.rd)$isMultiHit <- as.logical(cloneHits)
  rm(cloneHits)    
  cleanit <- gc()
  
  mcols(psl.rd)$pass.allQC <- mcols(psl.rd)$pass.percIdentity & 
    mcols(psl.rd)$pass.alignRatio & mcols(psl.rd)$pass.startWithin
  
  return(psl.rd)
}

#' Cluster/Correct values within a window based on their frequency given 
#' discrete factors
#'
#' Given a group of discrete factors (i.e. position ids) and integer values, 
#' the function tries to correct/cluster the integer values based on their 
#' frequency in a defined windowsize.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined. 
#' @param value a vector of integer with values that needs to corrected or 
#' clustered (i.e. Positions). Required if psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping of length posID or psl.rd by 
#' which to pool the rows (i.e. samplenames). Default is NULL. 
#' @param psl.rd a GRanges object returned from \code{\link{getIntegrationSites}}.
#' Default is NULL. 
#' @param weight a numeric vector of weights to use when calculating frequency 
#' of value by posID and grouping if specified. Default is NULL.
#' @param windowSize size of window within which values should be corrected or 
#' clustered. Default is 5.
#' @param byQuartile flag denoting whether quartile based technique should be 
#' employed. See notes for details. Default is TRUE.
#' @param quartile if byQuartile=TRUE, then the quartile which serves as the 
#' threshold. Default is 0.70.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}. 
#' Process is split by the grouping the column.
#' @param sonicAbund calculate breakpoint abundance using 
#' \code{\link{getSonicAbund}}. Default is FALSE.
#'
#' @note The algorithm for clustering when byQuartile=TRUE is as follows: for 
#' all values in each grouping, get a distribution and test if their frequency 
#' is >= quartile threshold. For values below the quartile threshold, test if 
#' any values overlap with the ones that passed the threshold and is within the 
#' defined windowSize. If there is a match, then merge with higher value, else 
#' leave it as is. This is only useful if the distribution is wide and polynodal. 
#' When byQuartile=FALSE, for each group the values within the defined window 
#' are merged with the next highest frequently occuring value, if freuquencies 
#' are tied then lowest value is used to represent the cluster. When psl.rd is 
#' passed, then multihits are ignored and only unique sites are clustered. All 
#' multihits will be tagged as a good 'clusterTopHit'.
#'
#' @return a data frame with clusteredValues and frequency shown alongside with 
#' the original input. If psl.rd parameter is defined then a GRanges object is 
#' returned with three new columns appended at the end: clusteredPosition, 
#' clonecount, and clusterTopHit (a representative for a given cluster chosen 
#' by best scoring hit!). 
#'
#' @seealso \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, 
#' \code{\link{otuSites}}, \code{\link{isuSites}}, \code{\link{crossOverCheck}}, 
#' \code{\link{pslToRangedObject}}, \code{\link{getSonicAbund}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' clusterSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-',
#' 'chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042), 
#' grouping=c('a','a','a','b','b','b','c'))
#' data(psl)
#' psl <- psl[sample(nrow(psl),100),]
#' psl.rd <- getIntegrationSites(pslToRangedObject(psl))
#' psl.rd$grouping <- sub("(.+)-.+","\\1",psl.rd$qName)
#' clusterSites(grouping=psl.rd$grouping, psl.rd=psl.rd)
#' }
clusterSites <- function(posID=NULL, value=NULL, grouping=NULL, psl.rd=NULL, 
                         weight=NULL, windowSize=5L, byQuartile=FALSE, 
                         quartile=0.70, parallel=TRUE, sonicAbund=FALSE) {

  # to avoid 'no visible binding for global variable' NOTE during R check #
  posID2 <- freq <- belowQuartile <- isMax <- isClosest <- val <- NULL
  ismaxFreq <- dp <- NULL
  
  .checkArgsSetDefaults_ALIGNed()
  
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
  } else {
    
    ## dereplicate & converge sites! ##
    if(!"Position" %in% colnames(mcols(psl.rd))) {
      stop("The object supplied in psl.rd parameter does not have Position ",
           "column in it. Did you run getIntegrationSites() on it?")
    }
    
    ## make sure column(s) for picking best hits are there! ##
    isthere <- grepl("score", colnames(mcols(psl.rd)), ignore.case=TRUE)  	
    if(!any(isthere)) {
      message("No 'score' column found in the data. Using 'qEnd' as an ",
              "alternative to pick the best hit!")
      isthere <- grepl("qEnd", colnames(mcols(psl.rd)), ignore.case=TRUE)		
      if(!any(isthere)) {
        stop("No 'qEnd' column found in the data either...can't pick the ",
             "cluster best hit without 'qEnd' or 'score' column present in ",
             "the object supplied in psl.rd :(")
      }            
    }
    
    if(length(which(isthere))>1) {
      message("multiple score based columns found: ", 
              paste(colnames(mcols(psl.rd))[which(isthere)], 
                    collapse=","), " choosing the first one...")
    }
    isthere <- which(isthere)[1]
    
    good.row <- rep(TRUE, length(psl.rd))
    multi.there <- grepl("isMultiHit", colnames(mcols(psl.rd)), 
                         ignore.case=TRUE)		
    if(any(multi.there)) { ## see if multihit column exists
      message("Found 'isMultiHit' column in the data. ",
              "These rows will be ignored for the calculation.")
      good.row <- good.row & !mcols(psl.rd)[[which(multi.there)]]
    }
    
    posIDs <- paste0(as.character(seqnames(psl.rd)), 
                     as.character(strand(psl.rd)))
    values <- mcols(psl.rd)$Position
    if(is.null(weight)) { 
      ## see if sequences were dereplicated before in the pipeline which 
      ## adds counts=x identifier to the deflines
      weight <- suppressWarnings(as.numeric(sub(".+counts=(\\d+)", "\\1",
                                                mcols(psl.rd)$qName)))
      if(all(is.na(weight))) { 
        weight <- NULL } 
      else { 
        weight[is.na(weight)] <- 1 
      }
    }
    
    grouping <- if(is.null(grouping)) { 
      rep("group1",length(values)) 
    } else { 
      grouping 
    }

    ## for sonic abundance ##
    mcols(psl.rd)$groups <- grouping
    
    clusters <- clusterSites(posID=posIDs[good.row], 
                             value=values[good.row], 
                             grouping=grouping[good.row], 
                             weight=weight[good.row], 
                             windowSize=windowSize, 
                             byQuartile=byQuartile, quartile=quartile)
    
    message("Adding clustered data back to psl.rd.")        
    clusteredValues <- with(clusters,
                            split(clusteredValue, 
                                  paste0(posID,value,grouping)))
    groupingVals <- paste0(posIDs, values, grouping)[good.row]
    mcols(psl.rd)$clusteredPosition <- mcols(psl.rd)$Position
    mcols(psl.rd)$clusteredPosition[good.row] <- 
      as.numeric(clusteredValues[groupingVals])
    
    ## add frequency of new clusteredPosition ##
    clusteredValueFreq <- with(clusters, 
                               split(clusteredValue.freq, 
                                     paste0(posID,value,grouping)))
    mcols(psl.rd)$clonecount <- 0
    mcols(psl.rd)$clonecount[good.row] <- 
      as.numeric(clusteredValueFreq[groupingVals])
    rm("clusteredValueFreq","clusteredValues","clusters")
    cleanit <- gc()
    
    ## pick best scoring hit to represent a cluster ##
    ## make sure to avoid multihit rows! ##
    message("Picking best scoring hit to represent a cluster.")   
    groupingVals <- paste0(as.character(seqnames(psl.rd)), 
                           as.character(strand(psl.rd)), 
                           mcols(psl.rd)$clusteredPosition, grouping)
    bestScore <- tapply(mcols(psl.rd)[[isthere]][good.row], 
                        groupingVals[good.row], max)
    isBest <- mcols(psl.rd)[[isthere]][good.row] == 
      bestScore[groupingVals[good.row]]
    
    ## pick the first match for cases where >1 reads with the same 
    ## coordinate had the same best scores ##
    tocheck <- which(isBest)
    res <- tapply(tocheck, names(tocheck), "[[", 1) 
    
    mcols(psl.rd)$clusterTopHit <- FALSE
    mcols(psl.rd)$clusterTopHit[good.row][res] <- TRUE
    mcols(psl.rd)$clusterTopHit[!good.row] <- TRUE
    
    message("Cleaning up!")
    rm("isBest","bestScore","posIDs","values","groupingVals")
    cleanit <- gc()
    
    if(sonicAbund) {
      message("Calculating sonic abundance.")   
      psl.rd <- getSonicAbund(psl.rd=psl.rd, grouping=mcols(psl.rd)$groups,
                              parallel=parallel)
    }
    mcols(psl.rd)$groups <- NULL
    
    return(psl.rd)
  }

  # get frequencies of each posID & value combination by grouping #
  groups <- if(is.null(grouping)) { "" } else { grouping }
  weight2 <- if(is.null(weight)) { 1 } else { weight }
  sites <- arrange(data.frame(posID, value, grouping=groups, 
                              weight=weight2, posID2=paste0(groups, posID), 
                              stringsAsFactors=FALSE), posID2, value)
  sites <- count(sites, c("posID","value","grouping","posID2"), wt_var="weight")    
  rm("groups","weight2")
  
  if(byQuartile) {
    message("Clustering by quartile: ", quartile)
    # obtain the defined quartile of frequency per posID & grouping #
    sites <- arrange(sites, posID2, value, plyr::desc(freq))
    quartiles <- with(sites,
                      tapply(freq, posID2, quantile, probs=quartile, names=FALSE))
    sites$belowQuartile <- with(sites,freq < quartiles[posID2])
    rm(quartiles)
    
    if(any(sites$belowQuartile)) {
      # for values belowQuartile, see if any within defined windowSize of 
      # aboveQuartile
      pos.be <- with(subset(sites,belowQuartile,drop=TRUE),
                     GRanges(IRanges(start=value,width=1), 
                             seqnames=posID2, freq=freq))
      pos.ab <- with(subset(sites,!belowQuartile,drop=TRUE),
                     GRanges(IRanges(start=value,width=1), 
                             seqnames=posID2, freq=freq))
      pos.overlap <- as.data.frame(as.matrix(findOverlaps(pos.be, pos.ab,
                                                          maxgap=windowSize,
                                                          ignore.strand=TRUE)))
      
      # for overlapping values, merge them with the biggest respective 
      # aboveQuartile site
      pos.overlap$freq <- values(pos.ab[pos.overlap[,"subjectHits"]])$freq
      pos.overlap$isMax <- with(pos.overlap, 
                                ave(freq, as.character(queryHits), 
                                    FUN=function(x) x==max(x)))
      pos.overlap$isMax <- as.logical(pos.overlap$isMax)
      pos.overlap <- subset(pos.overlap,isMax, drop=TRUE)
      
      # if there are >1 biggest respective aboveQuartile site, then choose the 
      # closest one ... if tied, then use the latter to represent the site
      counts <- xtabs(isMax~queryHits,pos.overlap)
      if(length(table(counts))>1) {
        toFix <- as.numeric(names(which(counts>1)))
        rows <- pos.overlap$queryHits %in% toFix
        pos.overlap$aboveQuartileValue <- 
          pos.overlap$belowQuartileValue <- pos.overlap$valueDiff <- 0
        pos.overlap$aboveQuartileValue[rows] <- 
          start(pos.ab[pos.overlap[rows,"subjectHits"]])
        pos.overlap$belowQuartileValue[rows] <- 
          start(pos.be[pos.overlap[rows,"queryHits"]])
        pos.overlap$valueDiff[rows] <- with(pos.overlap[rows,],
                                            abs(aboveQuartileValue-
                                                  belowQuartileValue))
        mins <- with(pos.overlap[rows,], 
                     tapply(valueDiff, as.character(queryHits), min))
        pos.overlap$isClosest <- TRUE
        pos.overlap$isClosest[rows] <- with(pos.overlap[rows,], 
                                            valueDiff == 
                                              mins[as.character(queryHits)])
        pos.overlap <- subset(pos.overlap, isMax & isClosest,drop=TRUE)
        rm("counts","mins")
      }
      
      # trickle the changes back to the original dataframe#     
      pos.overlap$clusteredValue <- 
        start(pos.ab[pos.overlap[,"subjectHits"]])    
      pos.overlap$posID2 <- 
        as.character(seqnames(pos.be[pos.overlap[,"queryHits"]]))
      
      # for cases where no overlap was found, try clustering to themselves #
      rows <- which(!1:length(pos.be) %in% pos.overlap$query)
      loners <- pos.be[rows]
      if(length(loners)>0) {
        times.rep <- values(loners)[["freq"]]
        res <- clusterSites(rep(as.character(seqnames(loners)), 
                                times=times.rep),
                            rep(start(loners),times=times.rep),
                            byQuartile=FALSE)
      }
      pos.overlap <- rbind(pos.overlap[,c("queryHits","clusteredValue")], 
                           data.frame(queryHits=rows, 
                                      clusteredValue=
                                        as.numeric(res$clusteredValue)))
      
      sites$clusteredValue <- sites$value
      sites$clusteredValue[sites$belowQuartile][pos.overlap[,"queryHits"]] <- 
        pos.overlap$clusteredValue
      stopifnot(any(!is.na(sites$clusteredValue)))            
    } else {
      message("No sites found below defined quartile. Try to increase ",
              "the quartile or use standard clustering, byQuartile=FALSE.")
    }
  } else {
    message("Clustering by minimum overlap.")
    
    sites <- split(sites, sites$grouping)
    
    sites <- bplapply(sites, function(x) {
      
      ## find overlapping positions using findOverlaps() using 
      ## maxgap adjusted by windowSize!
      sites.gr <- with(x, GRanges(seqnames=posID2, IRanges(start=value, width=1), 
                                  strand="*", freq))
      
      # the key part is drop.self=TRUE,drop.redundant=FALSE..
      # helps overwrite values at later step
      res <- as.data.frame(as.matrix(findOverlaps(sites.gr, drop.self=TRUE, 
                                                  drop.redundant=FALSE,
                                                  select="all", 
                                                  maxgap=windowSize))) 
      if(nrow(res)>0) {
        # add accessory columns to dictate decision making!
        # q = query, s = subject, val = value, freq = frequency of query/subject
        res$q.val <- start(sites.gr)[res$queryHits]
        res$s.val <- start(sites.gr)[res$subjectHits]
        res$q.freq <- sites.gr$freq[res$queryHits]
        res$s.freq <- sites.gr$freq[res$subjectHits]
        res$dist <- with(res,abs(q.val-s.val))
        
        ## do safety checking!
        stopifnot(!any(res$dist>windowSize)) 
        
        # favor a lower value where frequence/cloneCount is tied, 
        # else use the value of the highest frequency!
        res$val <- with(res, ifelse(q.freq==s.freq, 
                                    ifelse(q.val < s.val, q.val, s.val), 
                                    ifelse(q.freq >= s.freq, q.val, s.val))) 
        
        # For cases where there are >1 matches between query & subject...
        # find the one with the highest frequency and merge with that.
        # If all frequencies are the same, then use the lowest 
        # value to represent the cluster!
        res$maxFreq <- with(res, pmax(q.freq, s.freq))    
        res$ismaxFreq <- as.logical(with(res, ave(maxFreq, queryHits, 
                                                  FUN=function(x) x==max(x))))
        res$ismaxFreq <- as.logical(res$ismaxFreq)
        
        ## VIP step...this is what merges high value to low 
        ## value for ties in the hash structure below!!!
        res <- arrange(res, plyr::desc(queryHits), val)
        clustered <- unique(subset(res,ismaxFreq)[,c("queryHits","val")])
        clustered <- with(clustered, split(val, queryHits))
        
        ## make sure there is only one entry per hit this is useful in 
        ## situations when multiple query & subject are off by 1bp
        ## i.e. queryHits: 1,2,3,4; subjectHits: 1,2,3,4; 
        ## vals: 31895692 31895693 31895694 31895695
        clustered <- unlist(sapply(clustered, "[[", 1))        
        
        # trickle results back to sites
        x$clusteredValue <- x$value
        x$clusteredValue[as.numeric(names(clustered))] <- as.numeric(clustered)
        rm("clustered","res")
        cleanit <- gc()
      } else {
        message("No locations found within ", windowSize, "bps for ",
                x$grouping[1], "...no clustering performed!")
        x$clusteredValue <- x$value
      }    
      x
    }, BPPARAM=dp)
    sites <- rbind.fill(sites)
  }
  
  message("\t - Adding clustered value frequencies.")
  # get frequency of clusteredValue
  counts <- count(sites[,-grep("value",names(sites),fixed=TRUE)],
                  c("posID2","clusteredValue"), wt_var="freq")
  names(counts)[grep("freq",names(counts),fixed=TRUE)] <- "clusteredValue.freq"
  sites <- merge(sites,counts)
  
  if(byQuartile) {
    sites <- sites[,c("posID","value","freq","clusteredValue",
                      "clusteredValue.freq","grouping")]
  }
  
  sites$posID2<-NULL
  if(is.null(grouping)) { sites$grouping<-NULL }
  if(is.null(weight)) { sites$weight<-NULL }
  
  return(sites)
}

#' Bin values or make OTUs by assigning a unique ID to them within discrete factors.
#'
#' Given a group of values or genomic positions per read/clone, the function 
#' tries to yield a unique OTU (operation taxinomical unit) ID for the 
#' collection based on overlap of locations to other reads/clones by grouping. 
#' This is mainly useful when each read has many locations which needs to be 
#' considered as one single group of sites.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand).
#' Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, 
#' i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param readID a vector of read/clone names which is unique to each row, 
#' i.e. deflines.
#' @param grouping additional vector of grouping of length posID or psl.rd by 
#' which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a GRanges object returned from \code{\link{clusterSites}}. 
#' Default is NULL. 
#' @param maxgap max distance allowed between two non-overlapping position to 
#' trigger the merging. Default is 5.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}. 
#' Process is split by the grouping the column.
#'
#' @note The algorithm for making OTUs of sites is as follows: 
#' \itemize{
#'  \item for each grouping & posID, fix values off by maxgap parameter
#'  \item create bins of fixed values per readID
#'  \item assign arbitrary numeric ID to each distinct bins above & obtain its frequency
#'  \item perform overlap b/w readIDs with only one value (singletons) to readIDs with >1 value (non-singletons)
#'  \item   - for any overlapping values, tag non-singleton readID with the ID of singleton readID
#'  \item   - if non-singleton readID matched with more than one singleton readID, then pick on at random
#'  \item for any non-tagged & non-singleton readIDs, perform an overlap of values within themselves using the maxgap parameter
#'  \item   - tag any overlapping positions across any readID with the ID of most frequently occuring bin
#'  \item positions with no overlap are left as is with the original arbitrary ID
#' }
#' 
#' @return a data frame with binned values and otuID shown alongside the 
#' original input. If psl.rd parameter is defined, then a GRanges object.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{isuSites}},
#' \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, 
#' \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' otuSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), 
#' value=c(1000,1003,5832,1000,12324,65738,928042), 
#' readID=paste('read',sample(letters,7),sep='-'), 
#' grouping=c('a','a','a','b','b','b','c'))
otuSites <- function(posID=NULL, value=NULL, readID=NULL, grouping=NULL, 
                      psl.rd=NULL, maxgap=5, parallel=TRUE) {
  clusteredValue <- dp <- NULL
  .checkArgsSetDefaults_ALIGNed()
  
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
    stopifnot(!is.null(readID))
  } else {
    
    ## find the otuID by Positions ##
    isthere <- grepl("clusteredPosition", colnames(mcols(psl.rd)), 
                     ignore.case=TRUE)
    if(!any(isthere)) {
      message("No 'clusteredPosition' column found in the data. 
                    Using 'Position' as an alternative to use as values.")
      isthere <- grepl("Position", colnames(mcols(psl.rd)), 
                       ignore.case=TRUE)  
      if(!any(isthere)) {
        stop("The object supplied in psl.rd parameter does not have ",
             "'clusteredPosition' or 'Position' column in it. ",
             "Did you run getIntegrationSites() or clusterSites() on it?")
      }
    }
    
    good.rows <- TRUE
    value <- mcols(psl.rd)[[which(isthere)]]
    posID <- paste0(as.character(seqnames(psl.rd)), 
                    as.character(strand(psl.rd)))
    if("qName" %in% colnames(mcols(psl.rd))) {
      readID <- mcols(psl.rd)$qName
    } else if ("Sequence" %in% colnames(mcols(psl.rd))) {
      readID <- mcols(psl.rd)$Sequence
    } else {
      stop("No readID type column found in psl.rd object.")
    }
    
    grouping <- if(is.null(grouping)) { rep("A",length(psl.rd)) } else { grouping }
    
    otus <- otuSites(posID=posID[good.rows], 
                      value=value[good.rows], 
                      readID=readID[good.rows], 
                      grouping=grouping[good.rows], 
                      parallel=parallel)
    
    message("Adding otuIDs back to psl.rd.")        
    otuIDs <- with(otus,
                   split(otuID,
                         paste0(posID,value,readID,grouping)))
    mcols(psl.rd)$otuID <- NA
    mcols(psl.rd)$otuID[good.rows] <- 
      as.numeric(otuIDs[paste0(posID,value,readID,grouping)[good.rows]])
    
    message("Cleaning up!")
    rm("otus","otuIDs","value","posID","readID","grouping","good.rows")
    cleanit <- gc()
    
    return(psl.rd)
  }
  
  groups <- if(is.null(grouping)) { "A" } else { grouping }
  
  ## fix values off by maxgap parameter for sanity! ##
  sites.clustered <- clusterSites(posID, value, groups, windowSize=maxgap, 
                                  parallel=parallel)
  
  sites <- data.frame(posID, value, readID,                      
                      grouping=groups, stringsAsFactors=FALSE)
  sites <- merge(sites, sites.clustered)
  sites$posID2 <- with(sites, paste0(posID,clusteredValue))
  rm("groups","sites.clustered")
  
  ## get unique positions per readID by grouping 
  ## use tapply instead of ddply() or by() because it's a lot faster on 
  ## larger datasets
  counts <- with(sites, tapply(posID2, paste0(grouping,readID), 
                               function(x) {
                                 uniques <- sort(unique(x))
                                 list(paste(uniques,collapse=","), 
                                      length(uniques))
                               }))
  reads <- unique(sites[,c("grouping","readID")])
  reads$posIDs <- sapply(counts[with(reads,paste0(grouping,readID))],"[[", 1)
  reads$counts <- sapply(counts[with(reads,paste0(grouping,readID))],"[[", 2)
  
  # create initial otuID by assigning a numeric ID to each collection of 
  # posIDs per grouping
  reads$otuID <- unlist(
    lapply(lapply(with(reads,split(posIDs,grouping)), as.factor), as.numeric)
  ) 
  sites <- merge(arrange(sites,grouping,readID), 
                 arrange(reads[,c("grouping","readID","counts","otuID")],
                         grouping,readID), 
                 by=c("grouping","readID"), all.x=TRUE)
  sites$posID2 <- NULL
  rm(reads)
  sites <- arrange(sites, grouping, posID, clusteredValue)
  sites$readID <- as.character(sites$readID)
  sites$grouping <- as.character(sites$grouping)
  
  sites.gr <- with(sites, GRanges(seqnames=posID, IRanges(start=clusteredValue,
                                                          width=1), 
                                  strand="*", readID, grouping, counts, 
                                  otuID, newotuID=otuID))
  mcols(sites.gr)$grouping <- as.character(mcols(sites.gr)$grouping)
  mcols(sites.gr)$readID <- as.character(mcols(sites.gr)$readID)
  mcols(sites.gr)$check <- TRUE
  sites.gr <- sort(sites.gr)
  cleanit <- gc()
  
  ## see if readID with a unique/single location matches up to a readID with >1
  ## location, if yes then merge
  mcols(sites.gr)$singles <- mcols(sites.gr)$counts==1
  if(any(mcols(sites.gr)$singles)) {
    message('Merging non-singletons with singletons if any...')    
    sites.gr.list <- split(sites.gr, mcols(sites.gr)$grouping)
    sites.gr <- bplapply(sites.gr.list, function(x) {
      sigs <- subset(x, mcols(x)$singles)
      nonsigs <- subset(x, !mcols(x)$singles)
      res <- findOverlaps(nonsigs, sigs, maxgap=maxgap, select="first")
      rows <- !is.na(res)
      if(any(rows)) {
        res <- data.frame(queryHits=which(rows), subjectHits=res[rows])
        res$sigsOTU <- mcols(sigs)$otuID[res$subjectHits]
        
        res$sigsReadID <- mcols(sigs)$readID[res$subjectHits]
        res$nonsigsReadID <- mcols(nonsigs)$readID[res$queryHits]
        
        res$sigPosID <- paste0(as.character(seqnames(sigs)),
                               start(sigs))[res$subjectHits]
        res$nonsigPosID <- 
          paste0(as.character(seqnames(nonsigs)), start(nonsigs))[res$queryHits]
        
        ## if >1 OTU found per nonsigsReadID...choose lowest ID ## 
        bore <- with(res, split(sigsOTU, nonsigsReadID))
        bore <- sapply(sapply(sapply(bore, unique, simplify=FALSE), 
                              sort, simplify=FALSE), 
                       "[[", 1)
        res$OTU <- bore[res$nonsigsReadID]                                                        
        
        ## if >1 OTU found per nonsigPosID...choose lowest ID ##                          
        res$OTU2 <- res$OTU
        counts <- with(res, 
                       tapply(OTU2, nonsigPosID, function(x) length(unique(x)))
        )                              
        totest <- names(which(counts>1))
        while(length(totest)>0) { 
          #print(length(totest))
          for(i in totest) {
            rows <- res$nonsigPosID == i
            rows <- res$nonsigsReadID %in% res$nonsigsReadID[rows]
            rows <- rows | res$nonsigPosID %in% res$nonsigPosID[rows]
            res$OTU2[rows] <- min(res[rows,"OTU"])
          }
          counts <- 
            with(res, tapply(OTU2, nonsigPosID, function(x) length(unique(x))))                              
          totest <- names(which(counts>1))
        }
        
        bore <- sapply(with(res, split(OTU2, nonsigsReadID)), unique)
        
        if(!is.numeric(bore)) { 
          ## safety check incase >1 OTU found per readID
          bore <- sapply(bore, min)
        }
        rows <- mcols(nonsigs)$readID %in% names(bore)
        mcols(nonsigs)[rows,"newotuID"] <- bore[mcols(nonsigs)[rows,"readID"]]
        mcols(nonsigs)[rows,"check"] <- FALSE
        mcols(sigs)[mcols(sigs)$readID %in% res$sigsReadID, "check"] <- FALSE
      }
      c(sigs,nonsigs)
    }, BPPARAM=dp)
    sites.gr <- unlist(GRangesList(sites.gr), use.names=FALSE)
    rm(sites.gr.list)
  }
  mcols(sites.gr)$singles <- NULL
  
  ## perform non-singletons overlap of values within maxgap ##
  ## merge OTUs with overlapping positions within same grouping ##
  message('Performing non-singletons overlap...')
  goods <- subset(sites.gr, !mcols(sites.gr)$check)
  sites.gr <- subset(sites.gr, mcols(sites.gr)$check)
  sites.gr.list <- split(sites.gr, mcols(sites.gr)$grouping)
  sites.gr <- bplapply(sites.gr.list, function(x) {		    
    res <- findOverlaps(x, maxgap=maxgap, drop.self=TRUE, 
                        drop.redundant=TRUE, select="all")
    if(length(res)>0) {
      res <- as.data.frame(res)
      res$queryOTU <- mcols(x)$otuID[res$queryHits]
      res$subjectOTU <- mcols(x)$otuID[res$subjectHits]
      res$subjectReadID <- mcols(x)$readID[res$subjectHits]
      res$subjectPosID <- paste0(as.character(seqnames(x)), 
                                 start(x))[res$subjectHits]
      res$queryPosID <- paste0(as.character(seqnames(x)), 
                               start(x))[res$queryHits]
      
      ## if >1 OTU found per subjectReadID...choose lowest ID ## 
      bore <- with(res, split(queryOTU, subjectReadID))
      bore <- sapply(sapply(sapply(bore, unique, simplify=FALSE), 
                            sort, simplify=FALSE), "[[", 1)
      res$OTU <- bore[as.character(res$subjectReadID)]                                                    
      
      ## if >1 OTU found per subjectPosID...choose lowest ID ##                          
      res$OTU2 <- res$OTU                        
      counts <- with(res, 
                     tapply(OTU2, subjectPosID, function(x) length(unique(x))))                          
      totest <- names(which(counts>1))
      
      while(length(totest)>0) {
        #print(length(totest))
        for(i in totest) {
          rows <- res$subjectPosID == i 
          rows <- rows | res$subjectReadID %in% res$subjectReadID[rows]                               
          rows <- rows | res$subjectOTU %in% res$subjectOTU[rows]
          rows <- rows | res$subjectPosID %in% res$subjectPosID[rows]
          rows <- rows | res$queryOTU %in% res$subjectOTU[rows]
          
          res$OTU2[rows] <- min(res[rows,"OTU2"])
        }                            
        counts <- with(res, tapply(OTU2, subjectPosID, 
                                   function(x) length(unique(x)))
        )                          
        totest <- names(which(counts>1))
      }
      
      bore <- sapply(with(res, split(OTU2, subjectHits)), unique)
      ## safety check incase >1 OTU found per subjectHits
      stopifnot(is.numeric(bore)) 
      
      rows <- as.numeric(names(bore))
      mcols(x)[rows,"newotuID"] <- as.numeric(bore)
      mcols(x)[rows, "check"] <- FALSE
      
      bore <- sapply(with(res, split(OTU2, subjectReadID)), unique)
      if(!is.numeric(bore)) { 
        ## safety check incase >1 OTU found per readID
        bore <- sapply(bore, min)
      }
      rows <- mcols(x)$readID %in% names(bore)
      mcols(x)[rows,"newotuID"] <- bore[mcols(x)[rows,"readID"]]
      mcols(x)[rows,"check"] <- FALSE
    }
    x
  }, BPPARAM=dp)
  sites.gr <- unlist(GRangesList(sites.gr), use.names=FALSE)  
  sites.gr <- c(sites.gr, goods)
  rm("sites.gr.list","goods")
  cleanit <- gc()
  
  ## trickle the OTU ids back to sites frame ##    
  ots.ids <- sapply(split(mcols(sites.gr)$newotuID,
                          paste0(mcols(sites.gr)$readID, 
                                 mcols(sites.gr)$grouping)),
                    unique)
  
  if(!is.numeric(ots.ids) | any(sapply(ots.ids, length)>1)) {
    stop("Something went wrong merging non-singletons. ",
         "Multiple OTUs assigned to one readID most likely!")
  }
  sites$otuID <- as.numeric(unlist(ots.ids[with(sites,paste0(readID,grouping))],
                                   use.names=FALSE))
  
  stopifnot(any(!is.na(sites$otuID)))
  cleanit <- gc()
  
  if(is.null(grouping)) { sites$grouping <- NULL }
  
  sites$clusteredValue <- NULL; sites$clusteredValue.freq <- NULL
  sites$freq <- NULL; sites$counts <- NULL
  
  return(sites)
}

#' Bin values or make ISUs by assigning a unique ID to them within discrete factors.
#'
#' Given a group of values or genomic positions per read/clone, the function 
#' tries to yield a unique ISU (Integration Site Unit) ID for the collection 
#' based on overlap of locations to other reads/clones by grouping. This is 
#' mainly useful when each read has many locations which needs to be considered 
#' as one single group of sites.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand).
#' Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, 
#' i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param readID a vector of read/clone names which is unique to each row, 
#' i.e. deflines.
#' @param grouping additional vector of grouping of length posID or psl.rd by 
#' which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a GRanges object returned from \code{\link{clusterSites}}. 
#' Default is NULL. 
#' @param maxgap max distance allowed between two non-overlapping position to
#' trigger the merging. Default is 5.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}. 
#' Process is split by the grouping the column.
#'
#' @note The algorithm for making isus of sites is as follows: for each readID
#' check how many positions are there. Separate readIDs with only position from 
#' the rest. Check if any readIDs with >1 position match to any readIDs with 
#' only one position. If there is a match, then assign both readIDs with the 
#' same ISU ID. Check if any positions from readIDs with >1 position match any
#' other readIDs with >1 position. If yes, then assign same ISU ID to all 
#' readIDs sharing 1 or more positions.
#'
#' @return a data frame with binned values and isuID shown alongside the 
#' original input. If psl.rd parameter is defined, then a GRanges object 
#' where object is first filtered by clusterTopHit column and the isuID 
#' column appended at the end.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{isuSites}}, 
#' \code{\link{crossOverCheck}}, \code{\link{findIntegrations}}, 
#' \code{\link{getIntegrationSites}}, \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' isuSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-','chr16-','chr11-'), 
#' value=c(rep(1000,2),5832,1000,12324,65738,928042), 
#' readID=paste('read',sample(letters,7),sep='-'), 
#' grouping=c('a','a','a','b','b','b','c'))
isuSites <- function(posID=NULL, value=NULL, readID=NULL, grouping=NULL, 
                     psl.rd=NULL, maxgap=5, parallel=TRUE) {

  res <- otuSites(posID=posID, value=value, readID=readID, grouping=grouping, 
                   psl.rd=psl.rd, maxgap=maxgap, parallel=parallel)
  
  if(is(res,"GRanges") | is(res,"GAlignment")) {
    cols <- grep('otu',colnames(mcols(res)))
    colnames(mcols(res))[cols] <- gsub("otu","isu",colnames(mcols(res))[cols])
  } else {
    cols <- grep('otu',colnames(res))
    colnames(res)[cols] <- gsub("otu","isu",colnames(res)[cols])
  }
  
  res
}

#' Remove values/positions which are overlapping between discrete groups based 
#' on their frequency.
#'
#' Given a group of discrete factors (i.e. position ids) and integer values, 
#' the function tests if they overlap between groups. If overlap is found, 
#' then the group having highest frequency of a given position wins, else the 
#' position is filtered out from all the groups. The main use of this function 
#' is to remove crossover sites from different samples in the data.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). 
#' Required if psl.rd parameter is not defined.
#' @param value a vector of integer locations/positions that needs to be binned, 
#' i.e. genomic location. Required if psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping of length posID or psl.rd by 
#' which to pool the rows (i.e. samplenames). Default is NULL.
#' @param weight a numeric vector of weights to use when calculating frequency 
#' of value by posID and grouping if specified. Default is NULL.
#' @param windowSize size of window within which values should be checked. 
#' Default is 1.
#' @param psl.rd a GRanges object. Default is NULL. 
#'
#' @return a data frame of the original input with columns denoting whether a 
#' given row was a Candidate and isCrossover. If psl.rd parameter is defined, 
#' then a GRanges object with 'isCrossover', 'Candidate', and 'FoundIn' columns
#'  appended at the end.
#'
#' @seealso  \code{\link{clusterSites}}, \code{\link{otuSites}}, 
#' \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, 
#' \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' crossOverCheck(posID=c('chr1-','chr1-','chr1-','chr1-','chr2+','chr15-','
#' chr16-','chr11-'), value=c(rep(1000,3),5832,1000,12324,65738,928042), 
#' grouping=c('a','a','b','b','b','b','c','c'))
crossOverCheck <- function(posID=NULL, value=NULL, grouping=NULL, 
                           weight=NULL, windowSize=1, psl.rd=NULL) {
  
  .checkArgsSetDefaults_ALIGNed()
  
  # to avoid 'no visible binding for global variable' NOTE during R check #
  qgroup <- sgroup <- qfreq <- sfreq <- isBest <- NULL
  
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
  } else {     
    
    posID <- paste0(as.character(seqnames(psl.rd)), 
                    as.character(strand(psl.rd)))
    if("clusteredPosition" %in% colnames(mcols(psl.rd))) {
      message("Using clusteredPosition column from psl.rd as the value parameter.")
      value <- mcols(psl.rd)$clusteredPosition
      good.row <- mcols(psl.rd)$clusterTopHit
    } else if("Position" %in% colnames(mcols(psl.rd))) {
      message("Using Position column from psl.rd as the value parameter.")
      value <- mcols(psl.rd)$Position
      good.row <- rep(TRUE, length(value))
    } else {
      message("Using start(psl.rd) as the value parameter.")
      value <- start(psl.rd)
      good.row <- rep(TRUE, length(value))
    }
    
    isthere <- grepl("isMultiHit", colnames(mcols(psl.rd)), ignore.case=TRUE)
    if(any(isthere)) { ## see if multihit column exists
      message("Found 'isMultiHit' column in the data. ",
              "These rows will be ignored for the calculation.")
      good.row <- good.row & !mcols(psl.rd)[[which(isthere)]]
    }
    
    if(is.null(weight)) { ## see if clonecount column exists
      isthere <- grepl("clonecount", colnames(mcols(psl.rd)))
      if(any(isthere)) { 
        weight <- mcols(psl.rd)[[which(isthere)]] 
      } else { 
        weight <- weight 
      }
    }
    grouping <- if(is.null(grouping)) { "" } else { grouping }
    crossed <- crossOverCheck(posID[good.row], value[good.row], 
                              grouping=grouping[good.row], 
                              weight=weight[good.row],
                              windowSize=windowSize)
    
    message("Adding cross over data back to psl.rd.")
    crossed$pvg <- with(crossed, paste0(posID,value,grouping))
    totest <- pmin(xtabs(~pvg+isCrossover, crossed),1)
    if(any(rowSums(totest)>1)) {
      stop("Error in crossOverCheck: sampling culprits... ", 
           paste(rownames(totest[rowSums(totest)>1,]), collapse=", "))
    }
    
    crossed <- as(crossed, "DataFrame")
    rows <- match(paste0(posID,value,grouping)[good.row], crossed$pvg)
    
    newCols <- c("Candidate","isCrossover","FoundIn")
    mcols(psl.rd)[newCols] <- NA
    for(newCol in newCols) {
      mcols(psl.rd)[[newCol]] <- as(mcols(psl.rd)[[newCol]], 
                                    class(crossed[[newCol]]))
    }
    
    mcols(psl.rd)[good.row,][!is.na(rows),newCols] <- 
      crossed[rows[!is.na(rows)], newCols]
    
    message("Cleaning up!")
    rm("posID","value","grouping","crossed")
    cleanit <- gc()
    
    return(psl.rd)
  }
  
  # get frequencies of each posID & value combination by grouping #
  groups <- if(is.null(grouping)) { "" } else { grouping }
  weight2 <- if(is.null(weight)) { 1 } else { weight }
  sites <- data.frame(posID, value, grouping=groups, weight=weight2, 
                      stringsAsFactors=FALSE)
  sites <- arrange(sites, grouping, posID, value)
  sites <- count(sites, c("posID","value","grouping"), wt_var="weight")
  rm("groups","weight2")
  
  sites$isCrossover <- sites$Candidate <- FALSE
  sites$FoundIn <- sites$grouping
  
  # find overlapping positions & pick the winner based on frequencies #
  sites.gr <- with(sites, GRanges(seqnames=posID, IRanges(start=value,width=1), 
                                  strand="*", grouping, freq))
  res <- findOverlaps(sites.gr, maxgap=windowSize, drop.self=TRUE, 
                      drop.redundant=FALSE, select="all")
  if(length(res)>0) {
    res <- as.data.frame(res)
    res$qgroup <- mcols(sites.gr)$grouping[res$queryHits]
    res$sgroup <- mcols(sites.gr)$grouping[res$subjectHits]
    res$qfreq <- mcols(sites.gr)$freq[res$queryHits]
    res$sfreq <- mcols(sites.gr)$freq[res$subjectHits]
    res <- ddply(res, .(queryHits), summarize, 
                 FoundIn=paste(sort(unique(c(qgroup,sgroup))),collapse=","),
                 isBest=all(qfreq>sfreq))
    sites$isCrossover[subset(res,!isBest)$queryHits] <- TRUE  
    sites$FoundIn[res$queryHits] <- res$FoundIn
    sites$Candidate[res$queryHits] <- TRUE
  }
  sites
}

## a helper used in summary functions to obtain lengths by object class ##
.uniqueLength <- function(x) {
  if(is(x,"GRanges")) {
    length(unique(x$qName))
  } else if(is(x,"DNAStringSet")) {
    length(unique(names(x)))
  } else {
    length(x)
  }
}

#' Simple summary of a sampleInfo object.
#'
#' Give a simple summary of major attributes in sampleInfo/SimpleList object.
#'
#' @param object sample information SimpleList object, which samples per 
#' sector/quadrant information along with other metadata.
#' @param ... ignored for now.
#'
#' @return a dataframe summarizing counts of major attributes per sample 
#' and sector. 
#'
#' @export
#'
#' @examples 
#' data(FLX_seqProps)
#' sampleSummary(seqProps)
#'
sampleSummary <- function(object, ...) {
  stopifnot(is(object,"SimpleList"))
  message("Total sectors:", paste(names(object$sectors),collapse=","), "\n")

  res <- lapply(names(object$sectors), function(sector) {
    bore <- extractFeature(object, sector=sector,
                           feature="samplename")[[sector]]
    res.df <- data.frame(Sector=sector, SampleName=as.character(bore))
    res.df$SampleName <- as.character(res.df$SampleName)
    for (metaD in c("decoded","primed","LTRed","vectored","linkered",
                    "psl","sites")) {
      res <- extractFeature(object, sector=sector, feature=metaD)[[sector]]
      
      if(is(res,"DataFrame")) {
        res <- t(sapply(as.list(res), sapply, .uniqueLength))
        if(length(res)>0) {
          colnames(res) <- paste(metaD, colnames(res), sep=".")
          res.df <- cbind(res.df, as.data.frame(res)[res.df$SampleName,])
        } else {
          res.df[,metaD] <- NA
        }
      } else {
        res <- sapply(res, .uniqueLength)
        if(length(res)>0) {
          res.df[,metaD] <- res[res.df$SampleName]
        } else {
          res.df[,metaD] <- NA
        }
      }
    }
    res.df        
  })
  rbind.fill(res)
}

#' Calculate breakpoint/sonic abundance of integration sites in a population
#'
#' Given distinct fragment lengths per integration, the function calculates
#' sonic abundance as described in \code{\link{sonicLength}}. This function is 
#' called by \code{\link{clusterSites}} and needs all individual fragments 
#' lengths per position to properly estimate the clonal abundance of an 
#' integration sites in a given population.
#'
#' @param posID a vector of discrete positions, i.e. Chr,strand,Position.
#' Required if psl.rd parameter is not defined.
#' @param fragLen a vector of fragment length per posID. Required if 
#' psl.rd parameter is not defined. 
#' @param grouping additional vector of grouping of length posID or psl.rd by 
#' which to pool the rows (i.e. samplenames). Default is NULL.
#' @param replicateNum an optional vector of the replicate number per grouping 
#' and posID. Default is NULL.
#' @param psl.rd a GRanges object returned from \code{\link{getIntegrationSites}} 
#' Default is NULL.
#' @param parallel use parallel backend to perform calculation with 
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is 
#' registered, then a serial version is ran using \code{\link{SerialParam}}. 
#' Process is split by the grouping the column.
#'
#' @return a data frame with estimated sonic abundance shown alongside with the
#' original input. If psl.rd parameter is defined then a GRanges object is 
#' returned with a new column 'estAbund'.
#'
#' @note For samples isolated using traditional restriction digest method, 
#' the abundance will be inaccurate as it is designed for sonicated or sheared
#' sample preparation method.
#'
#' @seealso \code{\link{clusterSites}}, \code{\link{otuSites}}, 
#' \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}}, 
#' \code{\link{pslToRangedObject}}
#'
#' @export
#'
#' @examples 
#' data("A1",package='sonicLength')
#' A1 <- droplevels(A1[1:1000,])
#' bore <- with(A1, getSonicAbund(locations, lengths, "A", replicates))
#' head(bore)
getSonicAbund <- function(posID=NULL, fragLen=NULL, grouping=NULL, 
                          replicateNum=NULL, psl.rd=NULL, parallel=TRUE) {
  dp <- NULL
  .checkArgsSetDefaults_ALIGNed()
  
  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(fragLen))
  } else {
    
    ## find the abundance by Positions/ClusteredPosition ##
    isthere <- grepl("clusteredPosition", colnames(mcols(psl.rd)), 
                     ignore.case=TRUE)    
    if(!any(isthere)) {
      message("No 'clusteredPosition' column found in psl.rd.
              Using 'Position' as an alternative...")
      isthere <- grepl("Position", colnames(mcols(psl.rd)), ignore.case=TRUE)		
      if(!any(isthere)) {
        stop("No 'Position' column found in psl.rd. either...can't generate
             posID attribute :(")
      }            
    }
    isthere <- which(isthere)[1]
    
    if(!"qEnd" %in% colnames(mcols(psl.rd))) {
      stop("Supplied psl.rd object does not have qEnd column")
    }
    
    posIDs <- paste0(as.character(seqnames(psl.rd)), 
                     as.character(strand(psl.rd)), 
                     mcols(psl.rd)[[isthere]])
    fragLens <- mcols(psl.rd)$qEnd
    
    grouping <- if(is.null(grouping)) { 
      rep("A", length(psl.rd)) 
    } else { 
      grouping 
    }
    
    replicateNum <- if(is.null(replicateNum)) { 
      rep(1, length(psl.rd)) 
    } else { 
      replicateNum 
    }

    res <- getSonicAbund(posIDs, fragLens, grouping, replicateNum, 
                         parallel=parallel)
    
    message("Adding sonic abundance back to psl.rd.")        
    estAbund <- with(res, split(estAbund, paste(posID, grouping)))
    mcols(psl.rd)$estAbund <- 0
    mcols(psl.rd)$estAbund <- as.numeric(estAbund[paste(posIDs, grouping)])
    
    rm("res","estAbund","posIDs","fragLens","grouping","replicateNum")
    cleanit <- gc()
    return(psl.rd)
  }
  
  dfr <- data.frame(posID, fragLen, grouping, row.names=NULL, 
                    stringsAsFactors=FALSE)
  if(!is.null(replicateNum)) {
    dfr$replicateNum <- replicateNum
  } else {
    dfr$replicateNum <- 1
  }
    
  counts.fragLen <- count(count(dfr, c("grouping","posID","fragLen"))[,-4],
                          c("grouping","posID"))
  names(counts.fragLen)[3] <- "fragLenCounts"
  
  dfr <- unique(dfr)
  dfr <- split(dfr, dfr$grouping)
  
  res <- bplapply(dfr, function(x) {
    dummy.theta <- structure(rep(1,length(x$posID)), names=x$posID)
    if(length(unique(x$replicateNum))>1) {
      siteAbund <- tryCatch(with(x, estAbund(factor(posID), fragLen, 
                                             factor(replicateNum))),  
                            error = function(z) list("theta"=dummy.theta))     
    } else {
      siteAbund <- tryCatch(with(x, estAbund(factor(posID), fragLen)),  
                            error = function(z) list("theta"=dummy.theta))
    }
    x$estAbund <- round(siteAbund$theta)[x$posID]
    x
  }, BPPARAM=dp)
  
  res <- rbind.fill(res)
  rownames(res) <- NULL
  res <- merge(unique(res[,c("grouping","posID","estAbund")]), 
               counts.fragLen, all.x=TRUE)
  return(res)
}

#' Prepend name attribute of a list to DNAStringSet
#'
#' Given a named listed DNAStringSet object returned from \code{\link{extractSeqs}}, 
#' the function prepends the sample name to read names.
#'
#' @param dnaSet output from \code{\link{extractSeqs}}
#' @param flatten should the output be unlisted? Default is FALSE.
#'
#' @return listed DNAStringSet with the names attribute prepended with the 
#' name of the list. If flatten is TRUE, then a DNAStringSet object
#'
#' @seealso \code{\link{extractFeature}}, \code{\link{extractSeqs}}, 
#' \code{\link{getSectorsForSamples}}, \code{\link{write.listedDNAStringSet}}
#'
#' @export
#'
#' @examples 
#' load(file.path(system.file("data", package = "hiReadsProcessor"),
#' "FLX_seqProps.RData"))
#' samples <- c('Roth-MLV3p-CD4TMLVWell6-Tsp509I', 
#' 'Roth-MLV3p-CD4TMLVWell6-MseI', 'Roth-MLV3p-CD4TMLVwell5-MuA')
#' seqs <- extractSeqs(seqProps, sector='2', samplename=samples, 
#' feature="genomic")
#' addListNameToReads(seqs, TRUE)
addListNameToReads <- function(dnaSet, flatten=FALSE) {
  stopifnot(class(dnaSet) == "list")
  
  if(all(sapply(dnaSet, class)=="DNAStringSet")) {
    ## this would be list of DNAStringSets ##    
    out <- mapply(function(x, y) {
      names(x) <- paste(y, names(x), sep="-")
      x
    }, dnaSet, names(dnaSet), SIMPLIFY=FALSE)
  } else {
    ## this would be list of lists ##
    out <- sapply(dnaSet, function(i) {
      mapply(function(x, y) {
        names(x) <- paste(y, names(x), sep="-")
        x
      }, i, names(i), SIMPLIFY=FALSE)
    }, simplify=FALSE)
  }
  
  if(flatten) {
    names(out) <- NULL
    out <- do.call(c, unlist(out, use.names=FALSE))
  }
  
  out
}




#--------------------- intSiteRetriever -----------------

.get_unique_sites <- function(sample_ref, conn) {
    sample_ref_in_db <- .get_sample_ref_in_db(sample_ref, conn)
    sites <- tbl(conn, "sites") 
    inner_join(sites, sample_ref_in_db)
}

#' for a given samples get sites positions.
#' @param sample_ref dataframe with 2 cols: sampleName, refGenome
#' @param conn connection to database
#' @return sites dataframe with cols: siteID, chr, strand, position, sampleName, refGenome
#' @export
getUniqueSites <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sites <- .get_unique_sites(sample_ref, conn)
    collect( select(sites, 
        siteID, chr, strand, position, sampleName, refGenome),
        n = Inf
    )
}

.get_multihitpositions <- function(sample_ref, conn) {
    sample_ref_in_db <- .get_sample_ref_in_db(sample_ref, conn)
    multihitpositions <- tbl(conn, "multihitpositions") 
    inner_join(multihitpositions, sample_ref_in_db, by="sampleID")
}

#' lengths distributions for multihits
#' @inheritParams getUniqueSites
#' @return df with cols: sampleName, refGenome, multihitID, length
#' @export
getMultihitLengths <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    samples_multihitpositions <- .get_multihitpositions(sample_ref, conn)
    multihit_lengths <- tbl(conn, "multihitlengths")
    collect(distinct(select(inner_join(samples_multihitpositions, multihit_lengths, by="multihitID"),
        sampleName, refGenome, multihitID, length)),
        n = Inf)
}

.get_breakpoints <- function(sample_ref, conn) {
    sample_ref_sites <- .get_unique_sites(sample_ref, conn)
    breakpoints <- tbl(conn, "pcrbreakpoints") 
    inner_join(sample_ref_sites, breakpoints) 
}

#' breakpoints
#' @inheritParams getUniqueSites
#' @export
getUniquePCRbreaks <- function(sample_ref, conn) {
    breakpoints <- .get_breakpoints(sample_ref, conn)
    collect(select(breakpoints,
        breakpoint, count, position, siteID, chr, strand, sampleName, refGenome),
        n = Inf
    )
# column named kept as in DB ...sites.position AS integration...
}

.check_has_sample_ref_cols <- function(sample_ref) {
    return (all(c("sampleName", "refGenome") %in% names(sample_ref)))
}

.get_sample_table <- function(conn) {
    samples_in_db <- tbl(conn, "samples") 
    select(samples_in_db, sampleID, sampleName, refGenome)
}

.get_sample_ref_in_db <- function(sample_ref, conn) {
    samples_in_db <- tbl(conn, "samples") 
    samples_in_db <- select(samples_in_db, sampleID, sampleName, refGenome, gender)
    inner_join(samples_in_db, sample_ref, by=c('sampleName', 'refGenome'), copy=TRUE)
}

#' do we have samples in database
#' @inheritParams getUniqueSites
#' @return vector of TRUE/FALSE for each row in sample_ref df
#' @export
setNameExists <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    
    sample_ref_in_db <- collect(.get_sample_table(conn), n = Inf)
    if (nrow(sample_ref_in_db) == 0) { # nothing is in db
        return(rep(FALSE, nrow(sample_ref))) 
    }
    ids <- paste0(sample_ref$sampleName, sample_ref$refGenome)
    ids_DB <- paste0(sample_ref_in_db$sampleName, sample_ref_in_db$refGenome)
    ids %in% ids_DB
}

#' counts
#' @inheritParams getUniqueSites
#' @export
getUniqueSiteReadCounts <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites_breakpoints <- .get_breakpoints(sample_ref, conn) 
    sample_ref_sites_breakpoints_grouped <- group_by(
        sample_ref_sites_breakpoints, sampleName, refGenome)
    collect(summarize(sample_ref_sites_breakpoints_grouped, readCount=sum(count)), n = Inf)
}

#' unique counts for integration sites for a given sample(with fixed genome)
#' @inheritParams getUniqueSites
#' @export
getUniqueSiteCounts <- function(sample_ref, conn) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sample_ref_sites <- .get_unique_sites(sample_ref, conn)
    sample_ref_sites_grouped <- group_by(sample_ref_sites, sampleName, refGenome)
    collect(summarize(sample_ref_sites_grouped, uniqueSites=n()), n = Inf)
}


#' creates match random controls.
#' @inheritParams getUniqueSites
#' @param numberOfMRCs how many controls for each site
#' @return df with cols: siteID, position, strand, chr, sampleName, refGenome
#' @export
getMRCs <- function(sample_ref, conn, numberOfMRCs=3) {
    stopifnot(.check_has_sample_ref_cols(sample_ref))
    sites <- .get_unique_sites(sample_ref, conn) 
    sites.metadata <- collect(select(sites, 
        siteID, gender, sampleName, refGenome), n = Inf)

    sites_meta <- data.frame("siteID"=sites.metadata$siteID,
                           "gender"=tolower(sites.metadata$gender))

    stopifnot(length(unique(sites.metadata$refGenome)) == 1)
    ref_genome <- sites.metadata$refGenome[1] # all the same
  
    mrcs <- get_N_MRCs(sites_meta, get_reference_genome(ref_genome), numberOfMRCs)
  
    merge(mrcs, sites.metadata[c("siteID", "sampleName", "refGenome")])
}


#' return BS genome OBJECT(not name) for human readable UCSC format
#' 
#' format examples are: hg18, hg19, hg38 for human
#'                      mm8, mm9 for mouse
#' @note stop if cannot find unique genome from installed BSgenome
#' @seealso getRefGenome
#' @export
get_reference_genome <- function(reference_genome) {
    if ( ! requireNamespace("BSgenome", quietly = TRUE)) {
        stop("BSgenome needed for get_reference_genome function to work. Please install it.", 
            call. = FALSE)
    }
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

#' for a given reference genome and gender generate random positions 
#'
#' @param siteIDs vector of unique siteIDs for use as random seed
#' @param reference_genome BS object reference genome(@seealso get_reference_genome)
#' @param gender 'm' or 'f'
#' @param number_of_positions total number of random positions to generate
#' @param male_chr list of male-specific chromosomes prefixes(only 1 prefix is allowed at present)
#' @return dataframe with columns: 
#'      siteID(numeric), chr(character), strand(character), position(numeric)
#' @export
get_random_positions <- function(siteIDs, reference_genome, gender='m',
                                 number_of_positions=3, male_chr=c("chrY")){
  stopifnot(length(male_chr) == 1)
  stopifnot(length(gender) == 1)
  stopifnot(.check_gender(gender))

  chr_len <- seqlengths(reference_genome)
  stopifnot(any(grepl(male_chr, names(chr_len)))) # male chomosome is in genome
  chr_len <- .get_gender_specific_chr(chr_len, gender, male_chr)
  chr_len <- chr_len[names(chr_len) != "chrM"] #remove mitochondria

  cs <- c(0,cumsum(as.numeric(chr_len)))
  genomeLength <- max(cs)

  seed <- .Random.seed #don't want to screw up global randomness

  mrcs <- lapply(siteIDs, function(x){
    set.seed(x)
    rands <- round(runif(number_of_positions, 1, genomeLength*2)-genomeLength)
    cuts <- cut(abs(rands), breaks=cs, labels=names(chr_len))

    #outputs in format of "siteID", "chr", "strand", "position"
    data_frame(
        "siteID" = rep(x,number_of_positions),
        "chr" = as.character(cuts),
        "strand" = as.character(cut(sign(rands), breaks=c(-1,0,1), labels=c("-", "+"), include.lowest=T)),
        "position" = abs(rands) - cs[match(cuts, names(chr_len))])
  })

  .Random.seed <- seed #resetting the seed
    
  do.call(rbind, mrcs)
}

#' generate random controls for sites
#'
#' @param sites_meta dataframe with columns: siteID, gender
#' @param reference_genome BS genome object. All sites have the same genome.
#' @param number_mrcs_per_site number of MRCs to generate for each site 
#' @return dataframe with columns: siteID, chr, strand, position
#'
#' @note siteID are the same as given by sites_meta df
#' @export
get_N_MRCs <- function(sites_meta, reference_genome, number_mrcs_per_site=3, male_chr="chrY") {
    stopifnot(setequal(names(sites_meta), c("siteID", "gender")))
    stopifnot(number_mrcs_per_site > 0)
    stopifnot(.check_gender(sites_meta$gender))

    num_sites <- nrow(sites_meta)
    tot_num_mrcs <- num_sites * number_mrcs_per_site

    mrcs <- lapply(split(sites_meta, sites_meta$gender), function(sites){
      get_random_positions(sites$siteID, reference_genome, sites[1,"gender"],
                           number_mrcs_per_site, male_chr)
    })

    plyr::unrowname(do.call(rbind, mrcs))
}

# from vector of chromosome lengths with names creates vector for male or female
# @param all_chromosomes vector with length, names(all_chromosomes) are actual names of chromosome
.get_gender_specific_chr <- function(all_chromosomes, gender, male_chr) {
    stopifnot( ! is.null(names(all_chromosomes)))
    stopifnot(length(male_chr) == 1)
    if (gender == 'm') {
        return(all_chromosomes)
    }
    stopifnot(gender =='f')
    chromosome_names <- names(all_chromosomes)
    female_chromosomes <- all_chromosomes[ ! grepl(male_chr, chromosome_names)]
    female_chromosomes
}

# gender can only be male('m') or female('f')
.check_gender <- function(gender) {
    valid <- c('m', 'f')
    values <- unique(gender)
    all(values %in% valid)
}


#----------------- GCcontent -------------------------


getGCpercentage <- function(
    sites, column_prefix, window_size, reference_genome_sequence
){
    stopifnot(length(window_size) == length(names(window_size)))
    metadata <- mcols(sites)

    rangesToCalc <- .expand_trim_GRanges(sites, window_size)

   #seqs will take a lot of memory
    #could split at severe cpu time penelty
    seqs <- getSeq(reference_genome_sequence, rangesToCalc, as.character=F)

    letterFreqs <- letterFrequency(seqs, c("G", "C", "A", "T"))
    rm(seqs)

    GC <- letterFreqs[, c("G", "C")]
    ATGC <- letterFreqs[, c("A", "T", "G", "C")]

    gcContent <- rowSums(GC)/rowSums(ATGC)

    gcContent[!is.finite(gcContent)] <- NA #handled gracefully by pipeUtils

    gcContent <- DataFrame(matrix(gcContent, nrow=length(sites)))

    names(gcContent) <- paste(column_prefix, names(window_size), sep=".")

    mcols(sites) <- cbind(metadata, gcContent)

    sites
}

.expand_trim_GRanges <- function(sites, window_size) {
    nsites <- length(sites)
    strand(sites) = "+" #unimportant for GC and speeds up later calculations

    sites.seqinfo.original <- seqinfo(sites)
    isCircular(seqinfo(sites)) <- rep(FALSE, length(seqinfo(sites)))

    sites <- rep(sites, length(window_size))
    sites <- trim(suppressWarnings(flank(sites,
                                       rep(window_size/2, each=nsites),
                                       both=T)))

    seqinfo(sites) = sites.seqinfo.original
    sites
}

#------------------------ hotROCs ------------------------------------

##' ROC curve areas and variances estimated as per Delong et al, 1988.
##'
##' This is a utility function that users will probably not want to
##' invoke directly.
##' @title ROC.DDCP
##' @param response logical vector
##' @param indep matrix or \code{data.frame} of
##'     \code{length(response)} rows containing numeric values.
##' @seealso DeLong, Elizabeth R., David M. DeLong, and Daniel
##'     L. Clarke-Pearson. \dQuote{Comparing the areas under two or
##'     more correlated receiver operating characteristic curves: a
##'     nonparametric approach.}  Biometrics (1988): 837-845.
##' @return \code{list} with element \code{theta} giving the areas
##'     under the ROC curve for the variables in \code{indep} and
##'     element \code{var} giving the variance-covariance matrix of
##'     the estimates.
##' @author Charles Berry
ROC.DDCP<-
  function(response, indep)
{
    if(!is.logical(response))
        stop("response must be a boolean variable")
    if(any(is.na(response)))
        stop("NA's not allowed in response")
    indep <- as.matrix(indep)
    x <- indep[!response,  , drop = FALSE ]
    y <- indep[response,  , drop = FALSE ]
    m <- sum(!response)
    n <- sum(response)
    x.ranks <- apply(x,2,rank)
    y.ranks <- apply(y,2,rank)
    is.x <- rep( c(TRUE,FALSE), c( m, n ))
    xy.ranks <- lapply(1:ncol(indep),
                       function(z) split(rank(c(x[,z],y[,z])), is.x ))
    V.10 <- 1 - sapply( 1:ncol(indep),
                       function(z) 2*(xy.ranks[[z]][[2]] - x.ranks[, z ])/ n )/2
    V.01 <- sapply( 1:ncol(indep),
                   function(z) 2*(xy.ranks[[z]][[1]]-y.ranks[, z])/ m )/2
    theta.hat <- colMeans(V.10)
    S.10 <- var(V.10) 
    S.01 <- var(V.01)
    ## ensure S is posdef
    S <- S.10/nrow(x) + S.01/nrow(y) + 1e-10*diag(ncol(S.10))
    list(theta = theta.hat, var = S)
}


##' ROC curve areas, variances, and p-values for datasets in which
##' each case has a collection of matching controls.
##'
##' When data are collected under a scheme in which matched controls
##' are used to allow for biased selection and every case has \code{k}
##' controls, the ROC area is estimated as the average fraction of
##' controls whose variable value is less than that of its reference
##' case. Variances are computed from the usual sample variance
##' divided by the number of cases.  It is the responsiibility of the
##' caller to remove \code{NA} values from the arguments.
##' @title ROC.MRC - ROC areas for matched control data
##' @param response logical vector identifying cases or character
##'     vector or factor vector with \dQuote{insertion} marking the
##'     cases. Alternatively, a factor or character vector with
##'     elements matching \code{"insertion"} denoting the cases.
##' @param stratum a vector of \code{length(response)} elements whose
##'     unique values correspond to cases and their matched
##'     controls. There must be exactly one case and \code{k} controls
##'     in every stratum.
##' @param variables numeric a \code{matrix} or \code{data.frame} with
##'     \code{length(response)} rows and two or more columns.
##' @param origin \code{NULL} or a vector of \code{length(response)}
##'     elements to identify different data sources.
##' @param origin.levels optional character vector of origin levels
##' @param ragged.OK logical - \code{TRUE}, if differing numbers of
##'     MRCs of any one origin are acceptable.
##' @return \code{list} with elements \sQuote{ROC} giving a matrix of
##'     ROC curve areas, \sQuote{var} giving a list of variance
##'     matrices, and \sQuote{pvalues} giving a list of matrices
##'     containing pvalues for various contrasts.
##' @export
##' @examples
##' case <- rep(rep(c(TRUE,FALSE),c(1,3)),1000)
##' group <- rep(1:1000,each=4)
##' var1 <- case + rnorm(4000)
##' var2 <- case + rexp(4000)
##' var3 <- case + runif(4000,-1,1)
##' origin <- rep(factor(letters[1:4]),each=1000)
##' ROC.MRC(case,group,cbind(var1,var2,var3),origin)
##' @author Charles Berry
ROC.MRC <-
    function(response,stratum,variables,origin=NULL,origin.levels=NULL,
             ragged.OK=TRUE)
{
    if (any(is.na(variables))){
        res <- colSums(is.na(variables))>0
        stop("NA values not allowed. \nFound in:",
             paste(names(res)[res],collapse="\n\t"))
        }
    stopifnot(all(!is.na(response)))
    stopifnot(all(!is.na(origin)))
    stopifnot(all(!is.na(stratum)))
   if (!is.logical(response)) response <- response == "insertion"
    if (is.null(origin)) origin <- rep(1,length(response))
    if (is.null(origin.levels))
        origin.levels <- as.character(unique(origin))
    stopifnot(is.logical(response))
    nvars <- ncol(variables)
    origin.levels <-
        if (is.factor(origin))
            levels(origin)
        else
            unique(as.character(origin))
    phi.fun <-
        function(x)
    {
        ok.rows <- origin==x
        stratum <- stratum[ok.rows]
        response <- response[ok.rows]
        variables <- as.matrix(variables[ok.rows,])
        rstrata <- stratum[response]
        mstrata <- stratum[!response]
        stopifnot(all(!duplicated(rstrata)))
        nsites <- length(rstrata)
        cmtab <- table(stratum,response)
        nMRCs <- cmtab[1,1]
        if (any(cmtab[,1]!=nMRCs)){
            if (ragged.OK){
                morder <- order(mstrata)
                mindex <- split(which(!response)[morder],
                                mstrata[morder])
                lenMRCs <- lengths(mindex)
                nMRCs <- max(lenMRCs)
                mindex[lenMRCs<nMRCs] <-
                    lapply(mindex[lenMRCs<nMRCs],
                           function(x) c(x,rep(NA,nMRCs-length(x))))
                mindex <- unlist(mindex,use.names=FALSE)
            } else {
                stop("Differing Numbers of MRCs not allowed.")
            }
        } else {
            mindex <- which(!response)[order(mstrata)]
        }
        
        if (any(cmtab[,2]!=1)) stop("MRCs with no matching Integration Site.")
        phi <- 
            sweep(
                array(variables[mindex,],
                      dim=c(nMRCs, nsites, nvars),
                      dimnames=list(NULL,NULL,colnames(variables))),2:3,
                variables[which(response)[order(rstrata)],],
                function(x,y) (sign(y-x)+1)/2)
        colMeans(phi,na.rm=TRUE)
    }

    phi.list <-
        sapply(origin.levels,phi.fun,simplify=FALSE)
    rocz <- sapply(phi.list,colMeans)
    ## inflate the variance by 1e-10 to avoid diff/0.0 in Stats 
    roczVar <- lapply(phi.list,function(x) 1e-10*diag(ncol(x)) + var(x)/nrow(x))
    nullStats <- (rocz-0.5)^2/sapply(roczVar,diag)
    nullPvals <- pchisq(nullStats,df=1,lower.tail=FALSE)
    variableDiffs <-
        do.call(rbind,combn(nvars,2,function(x) rocz[x[1],]-rocz[x[2],],
                simplify=FALSE))
    variableDVars <-
        sapply(roczVar,
           function(x) combn(nvars,2,
                                 function(y) sum(x[y,y]*c(1,-1,-1,1))))
    variableDStats <- variableDiffs^2/variableDVars
    variablePvals <- pchisq(variableDStats,df=1,lower.tail=FALSE)
    attr(variablePvals,"whichRow") <- combn(nvars,2)
    if (length(origin.levels)>1){
        originDVars <- combn(roczVar,2,function(x) diag(x[[1]])+diag(x[[2]]))
        originDiffs <- combn(origin.levels,2,function(x) rocz[,x[1]]-rocz[,x[2]])
        originStats <- originDiffs^2/originDVars
        originPvals <- pchisq(originStats,df=1,lower.tail=FALSE)
        attr(originPvals,"whichCol") <- combn(length(origin.levels),2)
        rownames(originPvals) <- rownames(rocz)
    } else {
        originPvals <- matrix(NA,nrow=nvars,ncol=1L)
    }
    pvals <- list(op=originPvals,vp=variablePvals,np=nullPvals)
    list(ROC=rocz,var=roczVar,pvalues=pvals)
}


##' ROC curve areas, variances, and p-values for datasets in which the
##' cases have ordinary random controls (i.e. not matched on any
##' characteristic).  The ROC areas and corresponding variances are
##' estimated by the method of Delong et al.  The data may derive from
##' multiple datasets. Comparisons of different variables in a dataset
##' compare the ROC curve areas relative to the controls.  Comparisons
##' of different datasets uses only the responses from each dataset
##' and do not involve the random controls.  It is the responsiibility
##' of the caller to remove \code{NA} values from the arguments.
##'
##' A matrix of ROC curve areas and variances.
##' @title ROC.ORC - ROC area matrix
##' @param response logical vector identifying cases or character
##'     vector or factor vector with \code{"insertion"} marking the
##'     cases.
##' @param variables a \code{matrix} or \code{data.frame} with
##'     \code{length(response)} rows and two or more columns.
##' @param origin \code{NULL} or a vector of \code{length(response)}
##'     elements to identify different data sources.
##' @param origin.levels optional character vector of origin levels
##' @return \code{list} with elements \sQuote{ROC} giving a matrix of
##'     ROC curve areas, \sQuote{var} giving a list of variance
##'     matrices, and \sQuote{pvalues} giving a list of matrices
##'     containing pvalues for various contrasts.
##' @export
##' @examples
##' case <- rep(rep(c(TRUE,FALSE),c(1,3)),1000)
##' var1 <- case + rnorm(4000)
##' var2 <- case + rexp(4000)
##' var3 <- case + runif(4000,-1,1)
##' origin <- rep(factor(letters[1:4]),each=1000)
##' ROC.ORC(case,cbind(var1,var2,var3),origin)
##' @author Charles Berry
ROC.ORC <-
    function(response,variables,origin=NULL,origin.levels=NULL)
{
    if (any(is.na(variables))){
        res <- colSums(is.na(variables))>0
        stop("NA values not allowed. \nFound in:",
             paste(names(res)[res],collapse="\n\t"))
    }
    stopifnot(all(!is.na(response)))
    stopifnot(all(!is.na(origin)))    
    if (!is.logical(response)) response <- response == "insertion"
    if (is.null(origin))
        origin <- rep(1,length(response))
    if (is.null(origin.levels))
        origin.levels <- as.character(unique(origin))
    ## res provides all the ROC areas, variances for comparison to
    ## null, and for variable vs variable
    res <-
        lapply(origin.levels,
           function(lev){
                   resp <- response[origin==lev]
                   indep <- variables[origin==lev,]
                   ROC.DDCP(resp,indep)
           })
    ## res2 provides variances for origin to origin comparisons
    pairs <- combn(origin.levels,2)
    res2 <-
        apply(pairs,2,
          function(x) {
                  ok.rows <- response & (origin %in% x)
                  ROC.DDCP(origin[ok.rows]==x[1],
                           variables[ok.rows,])})
    rocz <- 
        do.call(cbind, lapply(res,"[[","theta"))
    colnames(rocz) <- origin.levels
    rownames(rocz) <- colnames(variables)

    nullVars <- sapply(res,function(x) diag(x$var))
    nullStats <- (rocz-0.50)^2/nullVars
    nullPvals <- pchisq(nullStats,df=1L,lower.tail=FALSE)
    ncv <- nrow(rocz)
    variableDVars <-
        sapply(res,
           function(resElt) combn(ncv,2,
                      function(x) sum(resElt$var[x,x]*c(1,-1,-1,1))))
    variableDiffs <-
        do.call(rbind,
                combn(ncv,2,
              function(x) rocz[x[1],]-rocz[x[2],],simplify=FALSE))
    variableDStats <- variableDiffs^2/variableDVars
    variablePvals <- pchisq(variableDStats,df=1,lower.tail=FALSE)
    attr(variablePvals,"whichRow") <- combn(ncv,2)

    originDiffs <- sapply(res2,"[[","theta")-0.50
    originVars <- lapply(res2,function(x) diag(x$var))
    originStats <- originDiffs^2 / do.call(cbind,originVars)
    originPvals <- pchisq(originStats,df=1,lower.tail=FALSE)
    rownames(originPvals) <- rownames(rocz)
    attr(originPvals,"whichCol") <- combn(length(origin.levels),2)
    list(ROC=rocz,
         var=
             list(within.origin=lapply(res,"[[","var"),
                  between.origin=lapply(res2,"[[","var")),
         pvalues=list(op=originPvals,vp=variablePvals,
              np=nullPvals))
}

##' Create an interactive ROC curve heatmap
##'
##'  A file \code{file.path(svg.file.base, "main.svg" ) } is created
##'  which can be viewed with and SVG viewer like the FireFox
##'  browser. There are many other files linked to that file that
##'  allow the user to inspect significance test results.
##' @title ROC curve heatmaps
##' @param roc.res.list An object such as returned by
##'     \code{\link{ROC.ORC}} or \code{\link{ROC.MRC}}
##' @param svg.file.base The name of the directory for in which to
##'     store the resulting files
##' @return \code{NULL} -- the function is run for its side effects.
##' @importFrom grDevices dev.off
##' @importFrom graphics axis box image layout mtext par text title
##' @importFrom stats pchisq var
##' @importFrom utils combn
##' @importFrom RSVGTipsDevice devSVGTips setSVGShapeURL setSVGShapeToolTip
##' @importFrom colorspace diverge_hcl
##' @export
##' @author Charles Berry
ROCSVG <-
    function(roc.res.list,svg.file.base='roc')
{
    ## Purpose: Produce ROC heatmap with dynamic  p-value display
    ## ----------------------------------------------------------------------
    ## Arguments: roc.res - object of ROC.strata()
    ##            file - where to save results
    ## ----------------------------------------------------------------------
    ## Author: CCB, Date:  7 Oct 2008, 13:14

    roc.res <- roc.res.list$ROC
    
    ## colormap:

    
    dcol <- diverge_hcl(21, c = c(100, 0), l = c(50, 90), power = 1)


    ## main svg file:

    ## mainfile <- paste(svg.file.base,"main.svg",sep='-')

    mainfile <- file.path(svg.file.base,"main.svg")

    ## open up left margin

    bigmarg <- c(3,10,4,2)+.1

      ### using a subdir for most svg files and putting the unadorned
    ### version above it causes some headaches with keeping button
    ### relative URLs pointing to the right place. hence the use of
    ### strip.dirname and mk.image(use.base=...)

    
    strip.dirname <- function(buttons) lapply(buttons,function(x) {if(dirname(x['URL'])==svg.file.base) x['URL'] <- basename(x['URL']);x})
    
    ## button for the right margin

    null.URL <- file.path(svg.file.base,"H50.svg")
    relative.null.URL <- basename(null.URL)
    button.list <- list(
                        c(text="<Show Plain Heatmap>",URL=mainfile,tiptitle="Click to:",tipdesc="Clear Annotations"),
                        c(text="<Compare to Area == 0.50>",URL=relative.null.URL,tiptitle="Test Each Area",tipdesc="vs Chance Discrimination")
                        )
    
    ## transform to make image 'look right'

    troc.res <- t(roc.res)
    roc.rows <- nrow(roc.res)
    roc.cols <- ncol(roc.res)

    ### generic call to image

    mk.image <- function(title.=NULL,use.base=TRUE){ # TRUE for subdir images
        image(1:roc.cols,1:roc.rows,troc.res,zlim=c(0,1),axes=FALSE,col=dcol,
              xlab="",ylab="")
            if (!is.null(title.)) title(title.,line=4)
        box()
        sapply(1:roc.rows,
               function(x) {
                    setSVGShapeToolTip(title="Compare rows to:",desc=rownames(roc.res)[x])
                   ru <- if (use.base) basename(rowURL[x]) else rowURL[x]
                   setSVGShapeURL( ru )
                   mtext(rownames(roc.res)[x],side=2,adj=1,at=x,las=1,line=1)
               })

        sapply(1:roc.cols,
               function(x) {
                    setSVGShapeToolTip(title="Compare columns to:",desc=colnames(roc.res)[x])
                   cu <- if (use.base) basename(colURL[x]) else colURL[x]
                   setSVGShapeURL( cu )
                   mtext(colnames(roc.res)[x],side=3,adj=0.5,at=x,las=1,line=1)
               })
    }

    mk.buttons <- function(blist,descend){

      ## button list is a list with each element a vector with
      ## elements named "text", "URL", "tiptitle","tipdesc"

      n.button <- length(blist)
      at.pos <- if (roc.rows<2*n.button) {
        seq(0,roc.rows, length=n.button+2)[-c(1,n.button+2)]
      } else {
        seq(roc.rows-0.5,by=-2,length=n.button)
      }
      for (i in 1:n.button){
        button <- blist[[i]]
        if (!is.null(button['tiptitle'])&nchar(button['tiptitle'])!=0){
          setSVGShapeToolTip(title=button['tiptitle'],desc=button['tipdesc'])
        }
        setSVGShapeURL( button['URL'] ) 
        mtext(button['text'],side=4,adj=0,at=at.pos[i],las=1,line=2,
              col='blue',cex=0.8)
      }
    }

    
    mk.stripe <- function(stripe.ht=1){
      new.mai <- par()$mai
      new.mai[1] <- stripe.ht * 0.65
      new.mai[3] <- 0.25
      par(mai=new.mai)
      image(seq(0,1,length=length(dcol)),1,matrix(seq(0,1,length=length(dcol))),
            zlim=c(0,1),axes=FALSE,col=dcol,
            xlab="",ylab="",main="Color Key",cex.main=0.8)
      axis(1)

    }
    
    rowURL <- file.path(svg.file.base,paste('row',1:roc.rows,'svg',sep='.'))
    colURL <- file.path(svg.file.base,paste('col',1:roc.cols,'svg',sep='.'))

    ### functions to add text overlays:
    
    do.overlay <- function(ov) do.call(text,mk.overlay(ov)) 
    mk.overlay <- function(ovmat) # ovmat is a character matrix
      c(
        subset( data.frame(x=as.vector(col(ovmat)),y=as.vector(row(ovmat)),labels=as.vector(ovmat)),
               labels!=""),
        adj=0.5)

    mk.stars <- function(pvmat){
      x <- array("",dim(pvmat))
      x[] <- as.character( cut( pvmat, c(0, 0.001, 0.01, 0.05, 1), c("***","**","*",""), include.lowest=TRUE))
      x
    }

### here set up the overlays:

    opvals <- roc.res.list$pvalues$op
    vpvals <- roc.res.list$pvalues$vp
    nullpvals <- roc.res.list$pvalues$np
    
    opstars <- mk.stars( opvals )
    vpstars <- mk.stars( vpvals )
    nullpstars <- mk.stars( nullpvals )

    matchCol <- function(x,indx) x[2:1,][x==indx] 
    isCol <- function(x,indx) colSums(x==indx)==1

    omasks <- lapply(1:roc.cols,function(x){
      res <- array("",dim(roc.res))
      res[,x] <- "--"
      wc <- attr(opvals,'whichCol')
      res[,matchCol( wc, x )] <- opstars[ , isCol( wc,x)  ]
      res
    })
 
     vmasks <- lapply(1:roc.rows,function(x){
      res <- array("",dim(roc.res))
      res[x,] <- "|"
      wr <- attr(vpvals,'whichRow')
      res[matchCol(wr,x),] <- vpstars[  isCol(wr,x), ]
      res
    })


    ## draw the images:
    
    ## base image

    ## allow enough space for row labels to fit
    max.rowchar <- max(nchar(rownames(roc.res)))+10
    char.inch <- c(0.10,0.15) ### this is a guess - cin seems not to work for SVGDevice
    right.padding <- 1.0 + char.inch[1]*max(sapply(button.list,function(x) nchar(x['text'])))
    par.args <- list(xpd=NA,mai=c(4,max.rowchar,8,2)*rep(rev(char.inch),2)+c(0,0,0,right.padding))
    
    ## be sure columns are wide enough

    max.colchar <- max(nchar(colnames(roc.res)))+1
    min.width <- max(6,max.colchar) * char.inch[1] * roc.cols

    ## be sure to make enough room for text or it will not display
    ## properly

    vpad <- 0.025*roc.rows
    dev.width <- ceiling(min.width + sum(par.args$mai[c(2,4)]))
    min.height <- ceiling(  char.inch[2]*roc.rows )
    dev.height <- vpad+min.height + sum(par.args$mai[c(1,3)])

    layout.args <- list( matrix(1:2,ncol=1), heights=c(dev.height,1) )

    dir.create( svg.file.base, showWarnings = FALSE )

    
    devSVGTips(file = mainfile, width=dev.width,height=dev.height+1)
    do.call(layout, layout.args)
    do.call(par, par.args )
    mk.image("ROC Curve Areas")
    mk.buttons( button.list[ -1 ] )
    mk.stripe()
    dev.off()

    ## populate subdir with overlays:
    
    for (i in seq(along=vmasks)){
      devSVGTips(file = rowURL[i],width=dev.width,height=dev.height+1)
      do.call(layout, layout.args)
      do.call( par, par.args )
      mk.image("Rows Compared",use.base=TRUE)
      mk.buttons( strip.dirname(button.list) )
      do.overlay( vmasks[[i]])
      mk.stripe()
      dev.off()
    }

    for (i in seq(along=omasks)){
      devSVGTips(file = colURL[i],width=dev.width,height=dev.height+1)
      do.call(layout, layout.args)
      do.call( par, par.args )
      mk.image("Columns Compared",use.base=TRUE)
      mk.buttons( strip.dirname(button.list) )
      do.overlay( omasks[[i]])
      mk.stripe()
      dev.off()
    }


    devSVGTips(file = null.URL,width=dev.width,height=dev.height+1)
    do.call(layout, layout.args)
    do.call( par, par.args )
    mk.image("Compare to Chance Discrimination",use.base=TRUE)
    mk.buttons( strip.dirname( button.list[ -grep(relative.null.URL,sapply(button.list,'[','URL')) ] ))
    do.overlay( nullpstars )
    mk.stripe()
    dev.off()
  
    
    invisible()

}

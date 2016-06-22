library(stringr, quietly = TRUE)
library(plyr, quietly = TRUE)
options(stringsAsFactors = FALSE)

#### make sure the folder is clean ####
if( length(list.files(".", pattern="*.RData", recursive=TRUE))>0 |
   length(list.files(".", pattern="Data/*.fasta", recursive=TRUE))>0 )
    stop("Previous run trace left, first run\ngit clean -df\n" )


#### run intSiteValidation test data and wait until finish ####
cmd <- 'Rscript ../../intSiteCaller.R -j intSiteValidation'
message(cmd)
if( system(cmd)!=0 ) stop(cmd, " not executed")
message("This test should finish in 10 minutes if the queues are not busy.\n")

#### track running process ####
minutes <- 0
mybjobsid <- function() {
    system("bjobs -w | grep intSiteValidation", intern=TRUE)
}
while( length(mybjobsid())>0 ) {
    message("Running_minutes: ", minutes, "\tjobs: ", length(mybjobsid()))
    Sys.sleep(60)
    minutes <- minutes + 1
}
message("Run stopped after: ", minutes, " minutes\n")

### 1. check md5 for RData objects ####
message("\nChecking md5 digest for RData files")
source("../../check_rdata_md5.R")

### 2. check attriton table ####
message("\nChecking attrition tables")
cmd <- "Rscript ../../check_stats.R > testrun.attr"
##message(cmd)
system(cmd)

attr.old <- read.table("intSiteValidation.attr", header=TRUE)
attr.old$workdir <- NULL
attr.new <- read.table("testrun.attr", header=TRUE)
attr.new$workdir <- NULL
if( !identical(attr.old, attr.new) ) {
    message("FAIL")
    q()
} else {
   message("PASS") 
}



q(save="no")


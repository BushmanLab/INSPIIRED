library('getopt');

opt = getopt(c("noRsetup", 's', 0, "logical"))

if (is.null(opt$noRsetup))
{
   # R packages to be installed
   r_packages <- c('Rcpp', 'dplyr', 'RMySQL', 'RSQLite', 'optparse', 'yaml', 'igraph', 'argparse', 'devtools', 'knitr', 'vegan', 'sonicLength', 'reldist', 'PubMedWordcloud', 'vcd', 'RSVGTipsDevice')

   install.packages(r_packages, repos='http://cran.us.r-project.org')

   source("http://bioconductor.org/biocLite.R")
   biocLite()
   biocLite(c('ShortRead', 'hiAnnotator', 'BSgenome', 'GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg18', 'BSgenome.Mmusculus.UCSC.mm9'), suppressUpdates=TRUE, ask=FALSE)

   library('devtools')
   install_github('BushmanLab/intSiteRetriever')
   install_github('BushmanLab/GCcontent')
}

# Download software repositories

system("git clone -b qsub https://github.com/BushmanLab/intSiteCaller")
system("git clone -b SQLite https://github.com/BushmanLab/intSiteUploader")
system("git clone -b sqlite https://github.com/BushmanLab/geneTherapyPatientReportMaker")
system("git clone -b sqlite https://github.com/BushmanLab/genomicHeatmapMaker")
system("git clone -b sqlite https://github.com/BushmanLab/EpigeneticHeatmapMaker")

setwd('input')
system("wget http://www.bushmanlab.org/assets/doc/INSPIIRED_demoDataSet.tar")
system("tar xvf INSPIIRED_demoDataSet.tar")
system("rm -f INSPIIRED_demoDataSet.tar")

setwd("../genomicHeatmapMaker")
system("wget http://www.bushmanlab.org/assets/doc/pipeUtils_1.3.5.tar.gz")
system("R CMD INSTALL pipeUtils_1.3.5.tar.gz")

setwd("../EpigeneticHeatmapMaker")
system("wget http://www.bushmanlab.org/assets/doc/INSPIIRED_EpigeneticData.tar")
system("tar xvf INSPIIRED_EpigeneticData.tar")
system("rm -f INSPIIRED_EpigeneticData.tar")

setwd("../..")
message('Setup complete')

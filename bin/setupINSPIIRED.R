#    This source code file is a component of the large INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
library('getopt');

opt = getopt(c("noRsetup", 's', 0, "logical"))

if (is.null(opt$noRsetup))
{
   # R packages to be installed
   r_packages <- c('Rcpp', 'dplyr', 'RMySQL', 'RSQLite', 'optparse', 'yaml', 'igraph', 'argparse', 'devtools', 'knitr', 'vegan', 'sonicLength', 'reldist', 'PubMedWordcloud', 'vcd', 'RSVGTipsDevice', 'gridExtra')

   install.packages(r_packages, repos='http://cran.us.r-project.org')

   source("http://bioconductor.org/biocLite.R")
   biocLite()
   biocLite(c('ShortRead', 'hiAnnotator', 'BSgenome', 'GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg18', 'BSgenome.Mmusculus.UCSC.mm9'), suppressUpdates=TRUE, ask=FALSE)

   library('devtools')
   install_github('BushmanLab/intSiteRetriever')
   install_github('BushmanLab/GCcontent')
}

# Download software repositories
setwd('components')
system("git clone -b qsub https://github.com/BushmanLab/intSiteCaller")
system("git clone https://github.com/BushmanLab/intSiteUploader")
system("git clone -b sqlite-sample-managment https://github.com/BushmanLab/geneTherapyPatientReportMaker")
system("git clone https://github.com/BushmanLab/genomicHeatmapMaker")
system("git clone https://github.com/BushmanLab/EpigeneticHeatmapMaker")

setwd('../inputs')
system("wget http://www.bushmanlab.org/assets/doc/INSPIIRED_demoDataSet.tar")
system("tar xvf INSPIIRED_demoDataSet.tar")
system("rm -f INSPIIRED_demoDataSet.tar")
system("ln -s ../../INSPIIRED.yml demoDataSet/INSPIIRED.yml")

setwd("../components/genomicHeatmapMaker")
system("wget http://www.bushmanlab.org/assets/doc/pipeUtils_1.3.5.tar.gz")
system("R CMD INSTALL pipeUtils_1.3.5.tar.gz")

setwd("../EpigeneticHeatmapMaker")
system("wget http://www.bushmanlab.org/assets/doc/INSPIIRED_EpigeneticData.tar")
system("tar xvf INSPIIRED_EpigeneticData.tar")
system("rm -f INSPIIRED_EpigeneticData.tar")

setwd("../..")
message('Setup complete')

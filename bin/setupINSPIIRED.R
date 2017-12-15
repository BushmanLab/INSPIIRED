#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
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

library('getopt');

opt = getopt(c("noRsetup", 's', 0, "logical"))

if (is.null(opt$noRsetup))
{
   library(devtools)
   install_version("reldist", version = "1.6-2", repos = "http://cran.us.r-project.org")
   install_version("sonicLength", version = "1.4.4", repos = "http://cran.us.r-project.org")
   install_version("PubMedWordcloud", version = "0.3.2", repos = "http://cran.us.r-project.org")
   install_version("vcd", version = "1.4-0", repos = "http://cran.us.r-project.org")
   install_version("RSVGTipsDevice", version = "1.0-4", repos = "http://cran.us.r-project.org")
   
   source("http://bioconductor.org/biocLite.R")
   biocLite()
   biocLite(c('hiAnnotator', 'BSgenome.Hsapiens.UCSC.hg18', 'BSgenome.Mmusculus.UCSC.mm9'), suppressUpdates=TRUE, ask=FALSE)
}

# Download software repositories
setwd('components')
system("git clone https://github.com/BushmanLab/intSiteUploader")
system("git clone -b intSiteCaller-deployment https://github.com/BushmanLab/intSiteCaller")
system("git clone -b geneTherapyPatientReportMaker-deployment https://github.com/BushmanLab/geneTherapyPatientReportMaker")
system("git clone -b genomicHeatmapMaker-deployment https://github.com/BushmanLab/genomicHeatmapMaker")
system("git clone -b EpigeneticHeatmapMaker-deployment https://github.com/BushmanLab/EpigeneticHeatmapMaker")

setwd('../inputs')
system("wget http://www.bushmanlab.org/assets/doc/INSPIIRED/INSPIIRED_demoDataSet.tar")
system("tar xvf INSPIIRED_demoDataSet.tar")
system("rm -f INSPIIRED_demoDataSet.tar")
system("ln -s ../../INSPIIRED.yml demoDataSet/INSPIIRED.yml")

setwd("../components/genomicHeatmapMaker")

setwd("../EpigeneticHeatmapMaker")
system("wget http://www.bushmanlab.org/assets/doc/INSPIIRED/INSPIIRED_EpigeneticData.tar")
system("tar xvf INSPIIRED_EpigeneticData.tar")
system("rm -f INSPIIRED_EpigeneticData.tar")

setwd("../..")
message('Setup complete')

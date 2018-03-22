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

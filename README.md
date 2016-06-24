## INSPIIRED

**Set up Conda using the Bash shell.**  
While installing Conda, agree to the license and agree to allow the setup script to update your .bashrc file.
```
%> wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
%> bash Miniconda2-latest-Linux-x86_64.sh
%> source ~/.bashrc
```

**Setup up INSPIIRED.**
````
%> git clone https://github.com/BushmanLab/INSPIIRED
%> cd INSPIIRED
%> conda env create -f bin/INSPIIRED.conda.yml
%> source activate INSPIIRED
%> export INSPIIRED=$(pwd)
%> Rscript bin/setupINSPIIRED.R
```

**Identify integration sites.**
```
%> cd inputs/demoDataSet
%> Rscript $INSPIIRED/components/intSiteCaller/intSiteCaller.R -j demo
%> Rscript $INSPIIRED/components/intSiteCaller/check_stats.R
```

**Upload data to local database.**   
The data will be uploaded to a SQLite database defined in the INSPIIRED.yml configuration file.
```
%> Rscript $INSPIIRED/components/intSiteUploader/intSiteUploader.R
```

**Create HTML patient report.**  
The report will be outputted to the analysis directory with the file name format (project).(patient id).(date).html
```
%> Rscript $INSPIIRED/components/geneTherapyPatientReportMaker/makeGeneTherapyPatientReport.R demo.csv
```

**Convert report from HTML to PDF.**  
* Note that the report name in the example below contains the current time and will be different on your system.
```
%> Rscript $INSPIIRED/components/geneTherapyPatientReportMaker/printReportToPdf.R SCIDn1.pP1.20160622.html
```

**Create an interactive genomic heat map.**  
The heat map files will be outputed to a directory named genomicHeatmap/.  
The heat map and associated files are all SVG (Scalable Vector Graphics) files which can be viewed with most browsers.
```
%> Rscript $INSPIIRED/components/genomicHeatmapMaker/genomic_heatmap_from_db.R -o genomicHeatmap demo.csv
```

**Create a epigenic heat map.**  
The heat map files will be outputted to a directory named epiHeatmap/.
```
%> Rscript $INSPIIRED/components/EpigeneticHeatmapMaker/epi_heatmap_from_db.R -o epiHeatmap  demo.epi.csv
```

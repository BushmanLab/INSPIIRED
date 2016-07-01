## INSPIIRED

INSPIIRED is a software suite designed to study viral integration sites and the longitudinal outcomes of gene therapy 
patients. Each component of the INSPIIRED software suite is available through GitHub. For detailed instructions on how  
to make use of each component, please visit its repository:  

[intSiteCaller](https://github.com/BushmanLab/intSiteCaller)  
[intSiteUploader](https://github.com/BushmanLab/intSiteUploader)  
[geneTherapyPatientReportMaker](https://github.com/BushmanLab/geneTherapyPatientReportMaker)  
[genomicHeatmapMaker](https://github.com/BushmanLab/genomicHeatmapMaker)  
[epigeneticHeatmapMaker](https://github.com/BushmanLab/EpigeneticHeatmapMaker)  

INSPIIRED is designed to run in a Linux environment and is best suited to run on a computational cluster where its  
memory requirements are proportional to the volume of data being processed. When using a computational cluster,  
update INSPIIRED's configuration file to instruct it to use either the bsub or qsub job submission methods.  
A single, multi-core computer may be used provided that qsub is installed in order to distribute INPIIRED's jobs  
across its cores. Alternatively, a single core computer (or multi-core machine without qsub installed) may be used  
by setting the configuration variable 'parallelize'  to 'no'.   


**Below is an example of how to set up INSPIIRED on a 64-bit Linux machine**  

**Set up Conda.**  
While installing Conda, agree to the license and agree to allow the setup script to update your .bashrc file.  
If you are not using a 64-bit Linux machine then please visit http://conda.pydata.org/miniconda.html  
and download the appropriate Python 2.x version of Conda for your machine. 
```
%> wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
%> bash Miniconda2-latest-Linux-x86_64.sh
%> source ~/.bashrc
%> conda config --add channels 'bioconda'  
%> conda config --add channels 'r'  
```

**Setup up INSPIIRED.**  
Running the 'git clone' command below will begin the installation process in the same directory from which it is called.    
INSPIIRED depends upon file paths defined in its configuration file (INPIIRED.yml) shown at the bottom of this page.  
If you install INSPIIRED in a location other than your home directory then you **must update** the paths in this configuration file.  
At the bottom of this page are known installation problems and their solutions.  
<br>
The setup script installs / updates a number of R libraries and may take up to 10 minutes to execute.  
````
%> git clone https://github.com/BushmanLab/INSPIIRED
%> cd INSPIIRED
%> export INSPIIRED=$(pwd)
%> conda env create -f bin/INSPIIRED.conda.yml
%> source activate INSPIIRED
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
(!) Note that the report name in the example below contains the current time and will be different on your system.
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
<br>

**INSPIIRED configuration file:**
```
# Log file path.
# This file will be written to the analysis directory.
logFile : intSiteCaller.log

# Database configuration.
# If the database parameter is set to 'sqlite' then the provided sqliteIntSitesDB and sqliteSampleManagement databases will be used.
# If the parameter is set to 'mysql' the [mysqlConnectionGroup] connection credentials defined in your ~/.my.cnf file will be used.
dataBase : sqlite
sqliteIntSitesDB : ~/INSPIIRED/databases/SQLite/intSites.db
sqliteSampleManagement : ~/INSPIIRED/databases/SQLite/specimenManagement.db
mysqlConnectionGroup :  intsites_miseq.read

# The parallelize parameter instructs INSPIRRED how to distribute jobs over the available processors.
# Allowed values are: no, qsub, bsub
# If 'no' is selected then all jobs will be ran serially.
# If using qsub on a single, multi-core machine, then set the forceQsubPath to 'Yes'.
parallelize : qsub
forceQsubPath : Yes

# Debuging
# If this parameter is set to Yes than a number of intermediate temp files will be retained.
debug : No

# Run id
# This id is used to distinguish the sequencing run from other previous runs.
runId : demoRun

# Maximum size for each sequence file chunk
chunkSize : 30000

# Epigenetic heat map data source
epigeneticDataDirectory : ~/INSPIIRED/components/EpigeneticHeatmapMaker/Epigenetic

# Log in to system hosting vector and sequencing data files
remoteUser : everett@microb120.med.upenn.edu

# Directory that holds the vector information file defiled in the vectorSeq column of the provided sampleInfo.csv file
vectorDataPath : .

# Directory that holds the R1, R2, and I1 sequencing run gzipped FASTQ files.
# Defined paths should be absolute or relative to the analysis directory.
SequencingFiles:
  I1 :  Data/Undetermined_S0_L001_I1_001.fastq.gz
  R1 :  Data/Undetermined_S0_L001_R1_001.fastq.gz
  R2 :  Data/Undetermined_S0_L001_R2_001.fastq.gz

# Processing parameters
ProcessingParameters:
  qualityThreshold     : '?'
  badQualityBases      : 5
  qualitySlidingWindow : 10
  mingDNA              : 20
  minPctIdent          : 95
  maxAlignStart        : 5
  maxFragLength        : 2500
```
<br>

####Known installation issues.
On some systems, the following error may arise with the setupINSPIIRED.R script when installing R libraries:
```
awk: symbol lookup error: libreadline.so.6: undefined symbol: PC  
```
This can only be resolved by replacing the conda library with a system library resource, i.e.
```
mv /home/everett/miniconda2/envs/INSPIIRED/lib/libreadline.so.6   /home/everett/miniconda2/envs/INSPIIRED/lib/__libreadline.so.6
ln -s /lib/x86_64-linux-gnu/libreadline.so.6   /home/everett/miniconda2/envs/INSPIIRED/lib/libreadline.so.6  
```
Please note that the system level libreadline.so.6 library may be in a different locatin on your system.

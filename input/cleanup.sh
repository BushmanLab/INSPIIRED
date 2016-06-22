qselect -u everett | xargs qdel
rm -rf *.RData Data/*.done Data/*.fasta *.log logs/* *.err *.done *.qsub Data/demultiplexedReps GTSP* *.2bit miseqid.txt Uninfected* clone* error*

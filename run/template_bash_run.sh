#!/bin/bash 

path_to_code=/work/lcvmm/cgDNAmc_parallel_rahul/cgDNAmc
runName=run_for_seq_NBR
paramset=$path_to_code/cgDNAparamset4.methMH.ME.txt
seqfile=Seq_NBR.txt 
 
nbrconf=100000 
drop=12 

$path_to_code/run_cgDNAmc_from_seq -e t0 -l $runName -i $seqfile -p $paramset -a $nbrconf -d $drop -j n  

mv "$runName"_t0_intr_* "$runName"_t0-intr 
mv "$runName"_t0_* "$runName"_t0 
rm "$runName"-*.log
#rm $seqfile 

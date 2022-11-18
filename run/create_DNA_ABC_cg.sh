sed  's/Prmset.txt/cgDNAparamset4.methMH.ME.txt/g' template_bash_run.sh_orig > template_bash_run.sh 

for i in {1..100}
#for i in {501..553}
do
sbatch Parallel.run ${i} DNA_ABC_CG
done

echo "cgDNAparamset4.methMH.ME.txt" >> DNA_ABC_cg/README


import shlex, subprocess, os, multiprocessing, sys

from joblib import Parallel, delayed
import numpy as np

sublist =  int(sys.argv[1]) ;

file_id = open(("../../SeqSubLists/Sublist_%d.txt" % (sublist) ),"r")
seq_file = file_id.readlines()
n_seq = len(seq_file )

inputs = range(n_seq) 

def processInput(line):
		
	# Extract the sequence  
	seq = str(seq_file[line])
	seq = seq[0:-1]
    
	GCends = "GC"
	seq = GCends + seq + GCends

	# Define the main Input/Output variable's name for cgDNAmc
	Seq_File_Name = ('Seq_%d.txt' % (line))

	# Write the sequence file text
	fid = open(Seq_File_Name,"w")
	fid.write("%s" % (seq))
	fid.close() 
    
    # Run .sh
	subprocess.call( ('cp ../../template_bash_run.sh bash_run_%d.sh' % (line) ) , shell=True )
	subprocess.call( ( 'sed -i \'s/NBR/%d/g\' bash_run_%d.sh' % (line,line) )  , shell=True )	
	subprocess.call(('./bash_run_%d.sh' % (line)))
    
	Run_File_Name = ("run_for_seq_%s" % (line))
    
	# Extract tangent-tangent correlation and compute l_p and l_d
	bp , ttc = np.loadtxt(( "%s_t0" % (Run_File_Name)), unpack=True)
	X = bp.T.dot(bp)
	y_lp = np.log(np.array(ttc).reshape(-1, 1)) 
	lp =-X/bp.T.dot(y_lp)
    
	ttc_intr = np.loadtxt(( "%s_t0-intr" % (Run_File_Name)), usecols=1)
	y_ld=y_lp-np.log(np.array(ttc_intr).reshape(-1, 1)) 
	ld = -X/bp.T.dot(y_ld)
    
	# Erase all the files for this simulation
	subprocess.call( ("rm bash_run_%d.sh %s %s_t0*"  % ( line, Seq_File_Name, Run_File_Name ) )  , shell=True )     
    
	# Return results 
	return [lp[0], ld[0]]
    
num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

np.savetxt(('results_%d.txt' % (sublist) ), results)


#README

-------------------------------------------------------
 cgNA+mc C++ (2022) --- developed at LCVMM, EPFL (https://lcvmwww.epfl.ch/)
-------------------------------------------------------  

cgNA+mc code is the evolution of the cgDNA Monte Carlo code developed in the LCVMM lab (supervised by Prof. John H. Maddocks) at EPFL.
This C++ package is originally developed by Jaroslaw Glowacki (jarek.glowacki@gmail.com) corresponding to the cgDNA model.
The current version cgNA+mc is extended by Alessandro patelli and Rahul Sharma.
Public version: last updated by Rahul Sharma (rs25.iitr@gmail.com, rahul.sharma@epfl.ch) on Nov 2022

-------------------------------------------------------

If you use this code in any publication, please cite
Mitchell, Jonathan S., Jaroslaw Glowacki, Alexandre E. Grandchamp, Robert S. Manning, and John H. Maddocks. "Sequence-dependent persistence lengths of DNA." Journal of chemical theory and computation 13, no. 4 (2017): 1539-1555.

------------------------------------------------------- 
cgNA+mc code 
This package allows for running Monte Carlo (MC) simulations of dsNA molecules using 
the cgNA+ nearest-neighbour rigid base and rigid phosphate model. 
The code allows for direct MC sampling from a Gaussian distribution. 

The underlying cgNA+ model is described in the nest section. 
------------------------------------------------------- 

cgNA+ model

cgNA+ is a software (Matlab or Octave or Python) package for
predicting the ground-state conformation and stiffness
matrix for double-stranded nucleic acids (dsNAs) fragment of any given sequence.
Free energy differences between two configurations of a
molecule of dsNA in standard environmental conditions, can then be computed.

The ground-state conformation is provided in the CURVES+
definition (for the bases and similar coordinates for phosphates) 
of the dsNA structural coordinates (both a
non-dimensional version and the original unscaled version),
and also as a PDB file of atomic coordinates. The PDB file
can be used in the program 3DNA to obtain 3DNA structural
coordinates if desired. The ground-state stiffness matrix
is provided for non-dimensional version of the CURVES+ helical coordinates.

#----------------------------------------------------------------------------
A user-friendly web version of this program is available at cgDNAweb.epfl.ch

If you use cgDNAweb.epfl.ch in relation to any publication please cite:

cgNA+web: A web based visual interface to the cgNA+ sequence dependent statistical mechanics model of double-stranded nucleic acids
R. Sharma, A. S. Patelli, L. De Bruin, and J.H. Maddocks
Submitted
The cgNA+web site is an evolution from a pre-cursor site (still available at https://cgdnaweb.epfl.ch/view2) and is described in:

cgDNAweb: a web interface to the cgDNA sequence-dependent coarse-grain model of double-stranded DNA.
L. De Bruin, J.H. Maddocks
Nucleic Acids Research 46, issue W1 (2018), p. W5-W10
DOI:10.1093/nar/gky351
#----------------------------------------------------------------------------

The current updated cgNA+web version is an interface to the enhanced coarse-grain model cgNA+ of sequence-dependent statistical mechanics of double-stranded nucleic acids. cgNA+ includes parameter sets for dsDNA in an epigenetic sequence alphabet, dsRNA, and DNA:RNA hybrid as described in detail in

cgNA+: A sequence-dependent coarse-grain model of double-stranded nucleic acids.
R. Sharma, EPFL Thesis #9792, Under the supervision of J. H. Maddocks
Download the PDF here https://lcvmwww.epfl.ch/publications/data/phd/18/PhD_thesis_final.pdf
The extended cgNA+ parameter sets are built on the cgDNA+ model, which itself extends the original cgDNA model by the inclusion of an explicit description of phosphate groups. The cite for the cgDNA+ model itself is:

A sequence-dependent coarse-grain model of B-DNA with explicit description of bases and phosphate groups parametrised from large scale Molecular Dynamics simulations.
A. S. Patelli, EPFL Thesis #9522, Under the supervision of J. H. Maddocks
Download the PDF here https://lcvmwww.epfl.ch/publications/data/phd/15/EPFL_TH9552.pdf
More generally the cgDNA family of models has its own web page, which includes citations to all other related codes and articles.

More information is available at

http://lcvmwww.epfl.ch/cgDNA
#----------------------------------------------------------------------------


%-------------------------------------------------------
% For the impatient...
%-------------------------------------------------------
Run the examples_i.py in Examples directory and also see the basic description of the functions in input.py. 
For more details, please read the codes in functions directory. 

%-------------------------------------------------------
% cgNA+mc package contents
%-------------------------------------------------------
Four basic directories:
1. Functions: contains all the necessary functions
2. Parametersets: Conatains the different parameter sets obtaind from different MD simulations. 
3. Examples: contains all the input files and output file for few example sequences. 
4. work: This is the working directory. You can change its name and just change the necessary 
	 updates in the input.py file. 



#---------------------------------------
#Details on parameter sets
#---------------------------------------

[dna_ps1] DNA PS1 (cgDNA+ model)
- Palindromic sequence library
- 3 microseconds of Amber MD time series
- bsc1 force field and SPC/E water model
- 150mM of K+ counter-ions (Dang parameters)
- Maximum entropy/likelihood truncation
- Fitting functional: Kullback-Leibler divergence with model pdf in first argument.
- Dinucleotide model with specific blocks for the dimers at the ends.

[dna_ps2] DNA PS2 (cgNA+ model, recommended)
- Palindromic sequence library
- 10 microseconds of Amber MD time series
- bsc1 force field, TIP3P water model
- 150mM of K+ counter-ions (Joung and Cheatham parameters)
- Maximum entropy/likelihood truncation
- Fitting functional: Kullback-Leibler divergence with model pdf in first argument.
- Dinucleotide model with specific blocks for the dimers at the ends.
- Also, contains parameters for modified CpG steps
- C is referred to as M and H when methylated and hydroxymethylated, respectively and
- G is referred to as N and K when complementary C is methylated and hydroxymethylated, respectively
- Note only CpG steps can be modified i.e. allowed steps are MN, MG, CN or HK, HG, CK 
- hydroxymethylated and methylated steps are allowed in the same sequence but not adjacent

[rna_ps2] RNA PS2 (cgNA+ model, recommended)
- Palindromic sequence library
- 10 microseconds of Amber MD time series
- OL3 force field, TIP3P water model
- 150mM of K+ counter-ions (Joung and Cheatham parameters)
- Maximum entropy/likelihood truncation
- Fitting functional: Kullback-Leibler divergence with model pdf in first argument.
- Dinucleotide model with specific blocks for the dimers at the ends.
- input accepts both U and T but then internally change T to U

[drh_ps2] DNA:RNA Hybrid (DRH) PS2 (cgNA+ model, recommended)
- Same sequence library but not palindromic
- 10 microseconds of Amber MD time series
- bsc1 and OL3 force field for DNA and RNA strand, respectively, TIP3P water model
- 150mM of K+ counter-ions (Joung and Cheatham parameters)
- Maximum entropy/likelihood truncation
- Fitting functional: Kullback-Leibler divergence with model pdf in first argument.
- Dinucleotide model with specific blocks for the dimers at the ends (only GC ends)
- Only accepts sequence in A, T, C, G and must be with GC ends


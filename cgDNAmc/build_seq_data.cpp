/*
 * Copyright 2016 Jaroslaw Glowacki
 * jarek (dot) glowacki (at) gmail (dot) com
 *
 * This file is part of cgDNAmc.
 *
 * cgDNAmc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cgDNAmc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cgDNAmc.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 * A program to compute <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a>
 * stiffness matrix @f$K@f$ (more exactly its lower triangular factor @f$L@f$,
 * of the Cholesky decomposition @f$K=LL^T@f$) and shape vector
 * @f$\widehat{\mathbf{w}}@f$ for a given DNA sequence using a given
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> parameter set (see below for
 * details). The data is saved in a file (under a
 * provided name) in a format defined by @ref save_cgdna_chol_stiff_and_w_hat()
 * and suitable for @ref load_cgdna_stiff_and_w_hat() in @ref run_cgDNAmc.cpp.
 * Section @ref sec_theory shows how the data is used.
 *
 * This code provides a reference of how to use some of the functions of the
 * module @ref grp_cgDNAutils.
 *
 * <b>Command line arguments</b>
 *
 * To run the binary 3 command line arguments need to be provided. If run
 * without command line arguments the binary prints usage instructions similar
 * to the following. If not all required arguments are provided the binary
 * prints information about the missing ones. The arguments can be provided in
 * the standard POSIX `getopt` format for command line options:
 *     - `-p` <BR>
 *     <b>required</b> <BR>
 *     the argument of this option is a path to a file with a
 *     <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> parameter set
 *     (see below)
 *     - `-s` <BR>
 *     <b>required</b> <BR>
 *     the argument of this option is a path to a file with a single string
 *     consisting of the DNA alphabet letters @f$\{\mathtt{A}, \mathtt{T},
 *     \mathtt{G}, \mathtt{C}\}@f$ -- the DNA sequence to construct the
 *     <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> stiffness matrix and
 *     shape vector for
 *     - `-o` <BR>
 *     <b>required</b> <BR>
 *     the argument of this option is a path to a file where the data should
 *     be stored
 *
 * A <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> parameter set is necessary
 * to build the stiffness matrix @f$K@f$ and shape vector
 * @f$\widehat{\mathbf{w}}@f$ for any sequence. The cgDNAmc package comes with
 * a parameter set
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/repo/cgDNAparamset2.txt">`cgDNAparamset2.txt`</a>
 * which is state of the art as of year 2015.
 * In order to facilitate use of other parameter sets as they become available
 * a Matlab/Octave script
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/repo/save_parameter_set.m">`save_parameter_set.m`</a>
 * is provided with the cgDNAmc package. As input the script takes a
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> parameter set given as a
 * Matlab/Octave structure of the format defined by the
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> package. The output is a file
 * with the requested name with the parameter set in a format defined by the
 * @ref load_cgdna_parameter_set() function inside @ref build_seq_data.cpp.
 *
 * <b>Example</b>
 *
 * The following example illustrates a particular use case.
 * Example input file used below and generated output file are
 * available in the directory
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/examples">`examples`</a>
 * separately or as a ZIP archive
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/examples/examples.zip">`examples.zip`</a>.
 * The names of the files below are links to those examples.
 *
 * If the compiled binary for build_seq_data.cpp and the parameter set file
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/repo/cgDNAparamset2.txt">`cgDNAparamset2.txt`</a>.
 * are located in the directory `/path/to/cgDNAmc` an example run could be:
 * @verbatim
 $ /path/to/cgDNAmc/build_seq_data  -p /path/to/cgDNAmc/cgDNAparamset2.txt  -s my_seq.txt  -o my_seq_data.txt
 @endverbatim
 * The result of the above run is a file called
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/examples/my_seq_data.txt">`my_seq_data.txt`</a>,
 * which contains:
 *     - the DNA sequence read from
 *     <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/examples/my_seq.txt">`my_seq.txt`</a>
 *     - <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> stiffness matrix
 *     @f$K@f$ (more precisely the lower triangular factor @f$L@f$ of the
 *     Cholesky decomposition @f$K=LL^{T}@f$ of the stiffness matrix).
 *     The @f$L@f$ factor is stored in triangular band form, with `UPLO = 'L'`,
 *     as described
 *     <a href="http://www.netlib.org/lapack/lug/node124.html">here</a>.
 *     - the <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> shape coefficients
 *     vector
 *
 * in the format defined by the @ref save_cgdna_chol_stiff_and_w_hat() function.
 * The file is ready for use with the @ref run_cgDNAmc.cpp code.
 *
 * @note An alternative way of generating sequence data for simulations is to use the
 * <a href="http://lcvmwww.epfl.ch/software/cgDNAmc/repo/save_seq_data.m">`save_seq_data.m`</a>
 * Matlab/Octave script to save data generated in a "non-standard" way.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
// The standard POSIX getopt function to handle command line options
#include <unistd.h>

#include <algebra3d/Vector3.h>
#include <algebra3d/Quaternion.h>
#include <cgDNArecon/reconstruct.h>

#include <helpers.h>
#include <cgDNAutils.h>

using namespace std;
using namespace algebra3d;
using namespace cgdna_recon;

/**
 * The main function of the program that builds cgDNA stiffness matrix and shape
 * vector for a given sequence using a given parameter set. The data is saved in
 * a file under the provided name.
 *
 * Run without command line arguments to see usage instructions.
 *
 * @param[in] argc The number of command line arguments
 * @param[in] argv The array of command line arguments
 * @return An error code; 0 - success
 */
int main(int argc, char **argv) {
    
	int result;

	vector<double> dimer_stiff[36];
	vector<double> dimer_shapes[36];
    vector<double> end1_stiff[16];
    vector<double> end1_shapes[16];
    vector<double> end2_stiff[16];
    vector<double> end2_shapes[16];
	
	string seq;

	vector<double> oligo_stiff;
	vector<double> oligo_shape;

	ifstream seq_stream;
	ofstream out_stream;

	string seq_filename;
	string params_filename;
	string out_filename;

	// Information for getopt
	int c;
	char options[] = { 'p', // parameter set file
			's', // sequence file
			'o' // output file
			};
	string options_msg[] =
			{
					"        the argument of this option is a path to a file with a cgDNA parameter\n"
							"        set to use\n", //
					"        the argument of this option is a path to a file with the DNA sequence\n"
							"        to compute cgDNA coefficients for; the sequence file should contain\n"
							"        a single string of the DNA letters A, T G or C describing the requested\n"
							"        sequence\n",
					"        the argument of this option is a path to a file where the computed\n"
							"        coefficients cgDNA data should\n"
							"        be stored\n" };
	string optionts_example[] = { "cgDNAparamset2.txt", //
			"my_seq.txt", //
			"my_seq_data.txt" };
	int num_opts = sizeof(options) / sizeof(options[0]);
	// All arguments are required
	int num_opts_required = num_opts;
	vector<bool> options_found(num_opts, false);
	string optionts_str = "";

	string binary = argv[0];

	// Check the number of input arguments
	if (argc == 1) {
		print_usage_info(binary, options, num_opts, num_opts_required,
				options_msg, optionts_example);
		return 1;
	}

	// Create an options string
	optionts_str = "";
	for (int i = 0; i < num_opts; ++i) {
		optionts_str += options[i];
		optionts_str += ":";
	}

	// Get the command line options with values
	while (true) {
		c = getopt(argc, argv, optionts_str.c_str());
		// Parameter set file
		if (c == options[0]) {
			params_filename = optarg;
			options_found[0] = true;
		}
		// DNA sequence file
		else if (c == options[1]) {
			seq_filename = optarg;
			options_found[1] = true;
		}
		// Output file
		else if (c == options[2]) {
			out_filename = optarg;
			options_found[2] = true;
		}
		// Incorrect option argument
		else if (c == -1) {
			break;
		}
		// Unrecognized argument or the end
		else {
			return c;
		}
	}
	// Unrecognized argument
	if (optind < argc) {
		cerr << "Error!\nRedundant command line argument: " << endl;
		cerr << argv[optind++] << endl;
		return num_opts + 1;
	}

	// Have all arguments been provided?
	if (!check_missing_opts(options, num_opts_required, options_found,
			options_msg)) {
		return num_opts + 2;
	}

	// Verify the command line arguments and initialize everything
	seq_stream.open(seq_filename.c_str());
	if (!seq_stream.good()) {
		cerr << "Couldn't open sequence file: '" << seq_filename << "'" << endl;
		seq_stream.close();
		return 1;
	}

	// Make sure the output file can be opened before any computations are done
	out_stream.open(out_filename.c_str());
	if (!out_stream.good()) {
		cerr << "Couldn't open output file for writing: '" << out_filename
				<< "'" << endl;
		seq_stream.close();
		out_stream.close();
		return 1;
	}
	out_stream.close();

	// Read the sequence
	seq_stream >> seq;
	seq_stream.close();
    

	// Read the provided parameter set into sparse symmetric matrix storage
	
    result = load_cgdna_parameter_set(params_filename, dimer_stiff,
                                      dimer_shapes, end1_stiff ,end1_shapes,end2_stiff,end2_shapes);

	if (result != 0) {
		return result;
	}
    

	// Build the stiffness matrix K and weighted shape vector c = K * \hat{w}
	result = build_cgdna_stiff_and_sigma(seq, dimer_stiff, dimer_shapes,
			end1_stiff,end1_shapes,end2_stiff,end2_shapes, oligo_stiff, oligo_shape);

	if (result != 0) {
		return result;
	}
    
   
    

	// Compute the Cholesky decomposition of K and use it to compute \hat{w}
	result = compute_cgdna_chol_stiff_and_w_hat(oligo_stiff, oligo_shape);

	if (result != 0) {
		return result;
	}

	// Save the data in a file with the provided name
	save_cgdna_chol_stiff_and_w_hat(out_filename, seq, oligo_stiff, oligo_shape,
			10);
    
    
}

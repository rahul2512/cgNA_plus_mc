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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <helpers.h>
#include <cgDNAutils.h>

using namespace std;


/*****************************************************************************/
int get_dimer_param_index(string dimer) {
	static const string BASES[4] = { "A", "C", "G", "T"};

	// Certainly incorrect
	if (dimer.size() != 2) {
		return -1;
	}

	dimer[0] = toupper(dimer[0]);
	dimer[1] = toupper(dimer[1]);

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			if (dimer == BASES[i] + BASES[j]) {
				return i * 4 + j;
			}
		}
	}
    
    if (dimer == "MN") {
        return 16;
    }else if (dimer == "NM") {
        return 17;
    }else if (dimer == "AM") {
        return 18;
    }else if (dimer == "TM") {
        return 19;
    }else if (dimer == "GM") {
        return 20;
    }else if (dimer == "CM") {
        return 21;
    }else if (dimer == "NT") {
        return 22;
    }else if (dimer == "NA") {
        return 23;
    }else if (dimer == "NC") {
        return 24;
    }else if (dimer == "NG") {
        return 25;
    }else if (dimer == "HI") {
        return 26;
    }else if (dimer == "IH") {
        return 27;
    }else if (dimer == "AH") {
        return 28;
    }else if (dimer == "TH") {
        return 29;
    }else if (dimer == "GH") {
        return 30;
    }else if (dimer == "CH") {
        return 31;
    }else if (dimer == "IT") {
        return 32;
    }else if (dimer == "IA") {
        return 33;
    }else if (dimer == "IC") {
        return 34;
    }else if (dimer == "IG") {
        return 35;
    }
    
    

	return -1;
}

/*****************************************************************************/
int get_dimer_end_param_index(string end_dimer, string endindex) {
    static const string BASES[4] = { "A", "C", "G", "T"};
    
    // Certainly incorrect
    if (end_dimer.size() <= 1) {
        return -1;
    } else if ( end_dimer.size() > 2){
      // Check format of the end label : "XX_endn", XX the dimer, n the endindex (5 for front (5'end), 3 for back (3'end))
        if (!((end_dimer.size()==7)&(end_dimer.substr(2,4) == "_end") &(end_dimer.substr(6,1) == endindex))){
            return -1;
        }
    }
    
    
    end_dimer[0] = toupper(end_dimer[0]);
    end_dimer[1] = toupper(end_dimer[1]);
    
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // the endimer label ends with "_endn" in paramset input file, but we only want the two bases
            if ((end_dimer.substr(0,1) == BASES[i]) & (end_dimer.substr(1,1) == BASES[j])) {
                return i * 4 + j;
            }
        }
    }
    
    return -1;
}



/*****************************************************************************/
int read_block_and_c(istream &param_stream, int size, vector<double> &stiff,
		vector<double> &shape) {
	double tmp;

	// The stiffness matrix
	// Here we read it directly into the LAPACK band storage
	// so that it is easy to build banded stiffness for a given sequence
	// The matrices are stored column major
	stiff.resize(size * size);
	for (int c = 0; c < size; ++c) {
		// Skip the elements until the main diagonal
		for (int r = 0; r < c; ++r) {
			param_stream >> tmp;
		}
		// Read everything from diagonal down
		for (int r = 0; r < size - c; ++r) {
			param_stream >> stiff[c * size + r];
		}
	}

	// The shape vector
	shape.resize(size);
	for (int c = 0; c < size; ++c) {
		param_stream >> shape[c];
	}
    

	return 0;
}

/*****************************************************************************/
int load_cgdna_parameter_set(string filename,
                             vector<double> dimer_stiff[36],
                             vector<double> dimer_shape[36],
                             vector<double> end1_stiff[16],
                             vector<double> end1_shape[16],
                             vector<double> end2_stiff[16],
                             vector<double> end2_shape[16],
                             std::ostream &log_stream) {

	string label;
	int ind;
	bool indRead[36];

	// Open the file
	ifstream param_file(filename.c_str());
	if (!param_file.good()) {
		log_stream << "Error!\nCouldn't open cgDNA parameter file" << endl;
		param_file.close();
		return 1;
	}

	for (int i = 0; i < 36; ++i) {
		indRead[i] = false;
	}

	// First comes the dimer data
	for (int i = 0; i < 36; ++i) {
		// First comes the label
		param_file >> label;
        
		ind = get_dimer_param_index(label);
    
		if (ind < 0) {
			log_stream << "Error!\nIncorrect dimer label: '" << label << "'"
					<< endl;
			param_file.close();
			return 2;
		}
		read_block_and_c(param_file, 18, dimer_stiff[ind], dimer_shape[ind]);

		indRead[ind] = true;
	}

	for (int i = 0; i < 36; ++i) {
		if (!indRead[i]) {
			log_stream << "Error!\nNot all dimer data given in '" << filename
					<< "'" << endl;
			param_file.close();
			return 3;
		}
		indRead[i] = false;
	}

	// Then comes the front end data
	for (int i = 0; i < 16; ++i) {
		// First comes the label
		param_file >> label;
		ind = get_dimer_end_param_index(label,"5");
		if (ind < 0) {
			log_stream << "Error!\nIncorrect dimer label: '" << label << "'"
					<< endl;
			param_file.close();
			return 4;
		}
		read_block_and_c(param_file, 18, end1_stiff[ind], end1_shape[ind]);
        
        
		indRead[i] = true;
	}

	for (int i = 0; i < 16; ++i) {
		if (!indRead[i]) {
			log_stream << "Error!\nNot all dimer data given in '" << filename
					<< "'" << endl;
			param_file.close();
			return 5;
		}
	}
    
    
    // Then comes the back end data
    for (int i = 0; i < 16; ++i) {
        // First comes the label
        param_file >> label;
        ind = get_dimer_end_param_index(label,"3");
        if (ind < 0) {
            log_stream << "Error!\nIncorrect dimer label: '" << label << "'"
            << endl;
            param_file.close();
            return 6;
        }
        read_block_and_c(param_file, 18, end2_stiff[ind], end2_shape[ind]);
        indRead[i] = true;
    }
    
    for (int i = 0; i < 16; ++i) {
        if (!indRead[i]) {
            log_stream << "Error!\nNot all dimer data given in '" << filename
            << "'" << endl;
            param_file.close();
            return 7;
        }
    }


	param_file.close();

	return 0;
}

/*****************************************************************************/
int build_cgdna_stiff_and_sigma(string seq,
                                vector<double> dimer_stiff[36],
                                vector<double> dimer_shapes[36],
                                vector<double> end1_stiff[16],
                                vector<double> end1_shapes[16],
                                vector<double> end2_stiff[16],
                                vector<double> end2_shapes[16],
                                vector<double> &oligo_stiff,
                                vector<double> &oligo_sigma, std::ostream &log_stream) {

	int dimer_param_ind;
	int num_bp = seq.size();
	int num_pars = 12 * num_bp - 6;

	// Check the sequence
	if (num_bp < 2) {
		log_stream << "Error!\nSequence too short: '" << seq
				<< "' (at least lenth 2 needed)" << endl;
		return 1;
	}
    
    
    
    if(get_dimer_end_param_index(seq.substr(0,2),"5") == -1) {
        log_stream << "Error!\nIncorrect dimer : '" << seq[0] << seq[1]
        << "' at front end "<<endl;
        return 2;
    }
    
	for (int i = 1; i < num_bp-2; ++i) {
		if (get_dimer_param_index(seq.substr(i,2)) == -1) {
			log_stream << "Error!\nIncorrect dimer : '" << seq[i] << seq[i+1]
					<< "' at position " << i + 1 << endl;
			return 3;
        }
	}
    
    if(get_dimer_end_param_index(seq.substr(num_bp-2,2),"3") == -1) {
        log_stream << "Error!\nIncorrect dimer : '" << seq[num_bp-2] << seq[num_bp-1]
        << "' at back end " << endl;
        return 4;
    }

    

	oligo_stiff.resize(num_pars * 18, 0.0);
	oligo_sigma.resize(num_pars, 0.0);
	// Set all to 0
	for (int i = 0; i < (int) oligo_stiff.size(); ++i) {
		oligo_stiff[i] = 0.0;
	}

	for (int i = 0; i < (int) oligo_sigma.size(); ++i) {
		oligo_sigma[i] = 0.0;
	}
    
    

	// The parameters are assumed to be already in the banded storage
	// so the banded oligomer stiffness is just given by simple summation
    
    // First the front end
        dimer_param_ind = get_dimer_end_param_index(seq.substr(0, 2),"5");
    
    
        //stiffness
        for (int j = 0; j < 18 * 18; ++j) {
            oligo_stiff[j] += end1_stiff[dimer_param_ind][j];
                   }
    
        // weighted shape vector c = K*\hat{w}
        for (int j = 0; j < 18; ++j) {
        oligo_sigma[j] += end1_shapes[dimer_param_ind][j];
            
        }
    
    // Then the middle
	for (int i = 1; i < num_bp - 2; ++i) {
        dimer_param_ind = get_dimer_param_index(seq.substr(i, 2));

		//stiffness
		for (int j = 0; j < 18 * 18; ++j) {
			oligo_stiff[18 * 12 * i + j] += dimer_stiff[dimer_param_ind][j];
		}


		//weighted shape vector c = K*\hat{w}
		for (int j = 0; j < 18; ++j) {
			oligo_sigma[12 * i + j] += dimer_shapes[dimer_param_ind][j];
		}
		
    }
    
    // Then the back end
        dimer_param_ind = get_dimer_end_param_index(seq.substr(num_bp-2, 2),"3");
    
        //stiffness
        for (int j = 0; j < 18 * 18; ++j) {
            oligo_stiff[18 * 12 * (num_bp-2) + j] += end2_stiff[dimer_param_ind][j];
        }
    
        //weighted shape vector c = K*\hat{w}
        for (int j = 0; j < 18; ++j) {
        oligo_sigma[12 * (num_bp-2) + j] += end2_shapes[dimer_param_ind][j];
        
        }
    
	return 0;
}

/*****************************************************************************/
int compute_cgdna_chol_stiff_and_w_hat(std::vector<double> &oligo_stiff,
		std::vector<double> &oligo_shape) {
	unsigned int num_pars = oligo_shape.size();
	int result;

	// Check if the sizes of the banded matrix storage and the weighted shape
	// vector are consistent
	if (num_pars % 12 != 6 && oligo_stiff.size() == num_pars * 18) {
		return 1;
	}

	// Compute Cholesky decomposition
	result = compute_cholesky_band(17, oligo_stiff);
	if (result != 0) {
		return result;
	}

	// Compute \hat{w} by solving K * \hat{w} = c
	result = solve_band_trian('N', 17, oligo_stiff, oligo_shape);
	if (result != 0) {
		return 2 * result;
	}
	result = solve_band_trian('T', 17, oligo_stiff, oligo_shape);

	return 3 * result;
}

/*****************************************************************************/
int save_cgdna_chol_stiff_and_w_hat(std::string filename,
		const std::string &description, const std::vector<double> &oligo_stiff,
		const std::vector<double> &oligo_w_hat, int num_digits,
		std::ostream &log_stream) {
	int num_pars = oligo_w_hat.size();
	int num_diags = oligo_stiff.size() / num_pars;
	int width;

	ofstream out_stream(filename.c_str());

	// Check that the file was opened correctly
	if (!out_stream.good()) {
		log_stream << "Error!\nCouldn't open the file'" << filename
				<< "' for writing" << endl;
		out_stream.close();
		return 1;
	}

	// There cannot be a new line in the description
	if (description.find('\n') != std::string::npos) {
		log_stream << "Error!\nDescrition with a new line character provided"
				<< endl;
		out_stream.close();
		return 2;
	}

	// The first line is just the provided description
	out_stream << description << endl;

	// Indicate size of data
	out_stream << num_pars << " " << num_diags << endl;

	out_stream << std::scientific;
	// We want num_digits including the one before . hence the -1
	out_stream.precision(num_digits - 1);
	width = num_digits + 6;

	// Print the matrix in column-major order
	// Subsequent rows will be columns of the band storage
	for (int i = 0; i < (int) oligo_stiff.size(); ++i) {
		out_stream << setw(width) << oligo_stiff[i];
		if (i % num_diags == num_diags - 1) {
			out_stream << "\n";
		} else {
			out_stream << " ";
		}
	}
	out_stream << endl;

	// Print the shape vector
	for (int i = 0; i < (int) oligo_w_hat.size(); ++i) {
		out_stream << setw(width) << oligo_w_hat[i] << "\n";
	}

	out_stream.close();

	return 0;
}

/*****************************************************************************/
int load_cgdna_stiff_and_w_hat(std::string filename, std::string &seq,
		std::vector<double> &oligo_stiff, std::vector<double> &oligo_w_hat,
		std::ostream &log_stream) {
	int num_cols;
	int num_rows;

	ifstream in_stream(filename.c_str());

	// Check that the file was opened correctly
	if (!in_stream.good()) {
		log_stream << "Error!\nCouldn't open sequence data file: '" << filename
				<< "'" << endl;
		in_stream.close();
		return 1;
	}

	// Fist comes the sequence line
	getline(in_stream, seq);

	// Read the size of data
	// In case of band storage num_cols is the dimension of the matrix
	// while num_rows is (band size + 1)
	in_stream >> num_cols >> num_rows;

	// Prepare storage for data
	oligo_stiff.resize(num_cols * num_rows);
	oligo_w_hat.resize(num_cols);

	// Read the matrix (column-major)
	for (unsigned int i = 0; i < oligo_stiff.size(); ++i) {
		in_stream >> oligo_stiff[i];
	}

	// Read the shape vector
	for (unsigned int i = 0; i < oligo_w_hat.size(); ++i) {
		in_stream >> oligo_w_hat[i];
	}

	in_stream.close();

	return 0;
}

/*****************************************************************************/
double eval_cgdna_en(const std::vector<double> &oligo_stiff,
		const std::vector<double> &oligo_w_hat, std::vector<double> config,
		std::ostream &log_stream) {
	unsigned int num_coords = config.size();

	// Check the sizes of the cgDNA marteix and vector
	if (num_coords % 12 != 6) {
		log_stream << "Error!\nIncorrect size of the configuration vector"
				<< endl;
		return 1;
	}

	// Check the size of the configuration vector
	if (oligo_w_hat.size() != num_coords) {
		log_stream
				<< "Error!\nDimensions of the ground state shape vector inconsistent\n"
				<< "with the configuration vector"
				<< endl;
		return 2;
	}

	// Check the size of the stiffness matrix
	if (oligo_stiff.size() != num_coords * 18) {
		log_stream
				<< "Error!\nDimensions of the stiffness matrix inconsistent\n"
				<< "with the configuration vector"
				<< endl;
		return 3;
	}

	// config <- config_orig - oligo_w_hat
	add_vecs(-1.0, oligo_w_hat, config);

	// config <- L^T * (config_orig - oligo_w_hat)
	//           = L^T * config
	mult_band_trian_mat_vec('T', 17, oligo_stiff, config);

	// Evaluate the energy as
	// (config_orig - oligo_w_hat)^T   *   K     * (config_orig - oligo_w_hat)
	// = (config_orig - oligo_w_hat)^T * L * L^T * (config_orig - oligo_w_hat)
	// =                          config^T * config
	return 0.5 * get_dot_prod(config, config);
}

/*****************************************************************************/
double compute_jacobian(const vector<double> &config) {
	double jacobian = 1.0;
	double tmp;
	int num_rot_pars = config.size() / 6 - 1;

	for (int i = 0; i < num_rot_pars; i++) {
		tmp = config[i * 6] * config[i * 6]
				+ config[i * 6 + 1] * config[i * 6 + 1]
				+ config[i * 6 + 2] * config[i * 6 + 2];

		tmp = (1.0 + 0.01 * tmp);

		jacobian /= tmp * tmp;
	}

	return jacobian;
}

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
 * Helper function for computations for the cgDNA model.
 */

#ifndef CGDNA_UTILS_H_
#define CGDNA_UTILS_H_

#include <string>
#include <vector>

/**
 * @defgroup grp_cgDNAutils Utility functions for the cgDNA model
 */

/**
 * @defgroup grp_paramset cgDNA parameter set loading and helper functions
 *
 * @ingroup grp_cgDNAutils
 */


/**
 * @ingroup grp_paramset
 *
 * Returns the index in the data structures of data referring to the dimer
 * described by the given string (2 characters). This function was added to
 * centralize the indexing to avoid errors.
 *
 * @param[in] dimer The string (2 characters, A, G, T or C) representing
 * the dimer (lower or upper case)
 * @return the index in the data structure of data referring to the base
 * described by the given character or -1 if the dimer is not recognized.
 */
int get_dimer_param_index(std::string dimer);

/**
 * @ingroup grp_paramset
 *
 * Returns the index in the data structures of data referring to the end dimer
 * described by the given string (2 characters). This function was added to
 * centralize the indexing to avoid errors.
 *
 * @param[in] end_dimer The string (2 characters, A, G, T or C) representing
 * the dimer (lower or upper case)
 * @param[in] string endindex: (5 for front (5'end), 3 for back (3'end))
 * @return the index in the data structure of data referring to the base
 * described by the given character or -1 if the dimer is not recognized.
 */
int get_dimer_end_param_index(std::string dimer, std::string endindex);



/**
 * @ingroup grp_paramset
 *
 * Reads a single pair of stiffness block and weighted shape vector of
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a>
 * parameter set data. The stiffness block is read directly into symmetric band
 * storage (`UPLO = 'L'`), as described
 * <a href="http://www.netlib.org/lapack/lug/node124.html">here</a>.
 * As a consequence half of the storage memory are 0s, but this simplifies
 * construction of stiffness matrix in band storage
 * (see @ref build_cgdna_stiff_and_sigma())
 *
 * @note This is just a helper function for @ref load_cgdna_parameter_set()
 *
 * @param[in] param_stream The stream to read from
 * @param[in] size The size of the shape vector and a stiffness block
 * @param[out] stiff The storage to read the block into. The block is read into
 * symmetric band storage directly
 * @param[out] shape The storage to read the weighted shape vector into.
 * @return 0 on successful execution, a non-zero value on error.
 */
int read_block_and_c(std::istream &param_stream, int size,
		std::vector<double> &stiff, std::vector<double> &shape);

/**
 * @ingroup grp_paramset
 *
 * Reads an entire <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> parameter set
 * from a text file formatted as:
 * @code
 // 16 times dimer data:
 label // (e.g AT)
 K // (an 18 x 18 block in standard format)
 c // (an 18-vector)
 // 4 times base data:
 label // (e.g G)
 K // (a 6 x 6 block in standard format)
 c // (a 6-vector)
 @endcode
 * Logs errors in the given log stream.
 *
 * @param[in] filename The name of the file to read from.
 * @param[out] dimer_stiff An array to store dimer stiffness blocks data.
 * @param[out] dimer_shape An array to store dimer weighted shape vectors data.
 * @param[out] end1_stiff  An array to store front end dimer stiffness blocks data.
 * @param[out] end2_stiff  An array to store back end dimer stiffness blocks data.
 * @param[out] end1_base  An array to store front end dimer weighted shape vectors data.
 * @param[out] end2_base  An array to store back end dimer weighted shape vectors data.
 * @param[out] log_stream The stream to log errors into.
 * @return 0 on successful execution, a non-zero value on error.
 *
 * @note the function uses @ref read_block_and_c() to read each pair of
 * stiffness block, weighted shape vector
 */
int load_cgdna_parameter_set(std::string filename,
		std::vector<double> dimer_stiff[36],
		std::vector<double> dimer_shape[36],
        std::vector<double> end1_stiff[16],
        std::vector<double> end1_shape[16],
        std::vector<double> end2_stiff[16],
        std::vector<double> end2_shape[16],
		std::ostream &log_stream = std::cerr);

/**
 * @defgroup grp_cgDNA_coeffs Building, saving and loading cgDNA stiffness and shape vector
 *
 * @ingroup grp_cgDNAutils
 */

/**
 * @ingroup grp_cgDNA_coeffs
 *
 * Builds the <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> stiffness matrix
 * (in symmetric band storage) and a weighted shape vector for the given DNA
 * sequence using the provided <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a>
 * parameter set data.
 *
 * @param[in] seq The DNA sequence to build the data for.
 * @param[in] dimer_stiff An array with
 * <a href="http:lcvmwww.epfl.ch/cgDNA">cgDNA</a> dimer stiffness blocks data.
 * @param[in] dimer_shape An array with
 * <a href="http:lcvmwww.epfl.ch/cgDNA">cgDNA</a> dimer weighted shape vectors data.
 * @param[in] base_stiff An array with
 * <a href="http:lcvmwww.epfl.ch/cgDNA">cgDNA</a> base stiffness blocks data.
 * @param[in] base_shape  An array with
 * <a href="http:lcvmwww.epfl.ch/cgDNA">cgDNA</a> base weighted shape vectors data.
 * @param[out] oligo_stiff Symmetric band storage for the computed stiffness
 * matrix (see <a href="http:www.netlib.org/lapack/lug/node124.html">LAPACK documentation</a>).
 * @param[out] oligo_sigma Storage for the computed weighted shape vector.
 * @param[out] log_stream The stream to log errors into (`std#cerr` by default).
 * @return 0 on successful execution, a non-zero value on error.
 */
int build_cgdna_stiff_and_sigma(std::string seq,
                                std::vector<double> dimer_stiff[36],
                                std::vector<double> dimer_shape[36],
                                std::vector<double> end1_stiff[16],
                                std::vector<double> end1_shape[16],
                                std::vector<double> end2_stiff[16],
                                std::vector<double> end2_shape[16],
                                std::vector<double> &oligo_stiff,
                                std::vector<double> &oligo_sigma,
                                std::ostream &log_stream = std::cerr);

/**
 * @ingroup grp_cgDNA_coeffs
 *
 * Turns a <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> stiffness matrix
 * @f$K@f$ into its lower triangular factor @f$L@f$ of Cholesky decomposition
 * @f$K=LL^{T}@f$ and a weighted shape vector @f$\mathbf{\sigma}@f$ into the
 * shape vector by solving the linear system with the computed factor.
 *
 * @param[inout] oligo_stiff On entry stiffness matrix in symmetric band storage
 * (see <a href="http://www.netlib.org/lapack/lug/node124.html">LAPACK documentation</a>);
 * on exit its @f$L@f$ Cholesky factor in the same form
 * @param[inout] oligo_shape On entry a weighted shape vector; on exit the shape
 * vector that is a result of a linear solve using the computed Cholesky
 * decomposition
 * @return 0 on successful execution, a non-zero value on error.
 */
int compute_cgdna_chol_stiff_and_w_hat(std::vector<double> &oligo_stiff,
		std::vector<double> &oligo_shape);

/**
 * @ingroup grp_cgDNA_coeffs
 *
 * Stores the representation of a <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a>
 * stiffness matrix (column major) and a shape vector into a text file in the
 * following format
 * @code
 sequence // the fist line is the sequence
 N M // The dimensions of the representation of stiffness matrix
     // N - number of columns, M -  number of rows (column major storage)
 K   // N lines of M numbers - a representation of the stiffness matrix
     // (column major storage)
 w   // the shape vector
 @endcode
 * @param[in] filename The name of the file to save the data into.
 * @param[in] seq The sequence the data refers to
 * @param[in] oligo_stiff The storage of the stiffness matrix to save.
 * @param[in] oligo_w_hat The shape vector to save.
 * @param num_digits The number of digits of precision to save the numbers with.
 * @param[out] log_stream The stream to log errors into.
 * @return 0 on successful execution, a non-zero value on error.
 */
int save_cgdna_chol_stiff_and_w_hat(std::string filename,
		const std::string &seq, const std::vector<double> &oligo_stiff,
		const std::vector<double> &oligo_w_hat, int num_digits,
		std::ostream &log_stream = std::cerr);

/**
 * @ingroup grp_cgDNA_coeffs
 *
 * Loads a DNA sequence and a representation of the stiffness matrix and a shape
 * vector (the inverse covariance and average of the Gaussian to draw from)
 * from a file in the format:
 * @code
 sequence // the fist line is the sequence
 N M // The dimensions of the representation of the stiffness
     // N - number of columns, M -  number of rows (column major storage)
 K   // N lines of M numbers - a representation of the stiffness matrix
     // (column major)
 w   // the shape vector of size N
 @endcode
 *
 * @note The matrix is assumed to be in column major order so `N` columns of
 * size `M` are stored in oligo_stiff, independent of their interpretation (for
 * interpretation see @ref run_cgDNAmc.cpp).
 *
 * @see save_cgdna_chol_stiff_and_w_hat()
 *
 * @param[in] filename The name of the file to read the data from.
 * @param[out] seq The sequence the data refers to.
 * @param[out] oligo_stiff The storage for the representation of the stiffness
 * matrix.
 * @param[out] oligo_w_hat The storage for shape vector to read.
 * @param[out] log_stream The stream to log errors into.
 * @return 0 on successful execution, a non-zero value on error.
 */
int load_cgdna_stiff_and_w_hat(std::string filename, std::string &seq,
		std::vector<double> &oligo_stiff, std::vector<double> &oligo_w_hat,
		std::ostream &log_stream = std::cerr);

/**
 * @defgroup grp_cgDNA_en Evaluating the cgDNA energy for a given configuration
 *
 * @ingroup grp_cgDNAutils
 */

/**
 * @ingroup grp_cgDNA_en
 *
 * Evaluates the cgDNA energy for a given configuration of an oligomer using the
 * ground state shape vector and stiffness matrix as loaded by
 * @ref load_cgdna_stiff_and_w_hat.
 *
 * @param[in] oligo_stiff The  representation of the stiffness matrix.
 * @param[in] oligo_w_hat The ground state shape vector.
 * * @param[in] config The config vector to evaluate the energy for.
 * @param[out] log_stream The stream to log errors into.
 * @return The computed value of the cgDNA energy.
 *
 * @see load_cgdna_stiff_and_w_hat()
 */
double eval_cgdna_en(const std::vector<double> &oligo_stiff,
		const std::vector<double> &oligo_w_hat,
		std::vector<double> config,
		std::ostream &log_stream = std::cerr);

/**
 * @defgroup grp_conf_func Functions of cgDNA configuration
 *
 * @ingroup grp_cgDNAutils
 */

/**
 * @ingroup grp_conf_func
 *
 * Computes the Jacobian of the provided
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> configuration.
 *
 * @param[in] config The configuration to compute the Jacobin for.
 * @return The Jacobian of the provided
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> configuration.
 */
double compute_jacobian(const std::vector<double> &config);

#endif /* CGDNA_UTILS_H_ */

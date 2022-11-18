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
 * General helper function for the Monte Carlo simulations
 */

#ifndef HELPERS_H_
#define HELPERS_H_

#include <stdint.h>
#include <vector>

#include <rand.h>

// Forward declarations
/**
 * @namespace algebra3d
 * @brief Dependency: Code for algebra of 3D transformations.
 *
 * This namespace is defined by the algebra3d project at
 * http://lcvmwww.epfl.ch/software/algebra3d.
 */
namespace algebra3d {
/**
 * @brief A class to store vectors of 3 elements.
 *
 * This class is defined by the algebra3d project at
 * http://lcvmwww.epfl.ch/software/algebra3d.
 */
template<class fpType> class Vector3;
typedef Vector3<double> Vector3D;
/**
 * @brief A class to store quaternions.
 *
 * This class is defined by the algebra3d project at
 * http://lcvmwww.epfl.ch/software/algebra3d.
 */
template<class fpType> class Quaternion;
typedef Quaternion<double> QuaternionD;
} // namespace algebra3d

// Excluded from Doxygen documentation
/** @cond */

// Forward declarations of the BLAS and LAPACK Fortran functions
extern "C" {
void daxpy_(int *n, const double *alpha, const double *x, int *incx, double *y,
		int *incy);

double ddot_(int *n, const double *x, int *incx, const double *y, int *incy);

void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
		const double *A, const int *lead_dim_A, const double *x,
		const int *increm_x, const double *beta, const double *y,
		const int *increm_y);

void dpbtrf_(const char *up_lo, const int *n, const int *kd, double *AB,
		const int *lead_dim_AB, int *info);

void dtbtrs_(const char *up_lo, const char *trans, const char *diag,
		const int *n, const int *kd, const int *num_rhs, const double *AB,
		const int *lead_dim_AB, double *B, const int *lead_dim_B, int *info);

void dtbmv_(const char *up_lo, const char *trans, const char *diag,
		const int *n, const int *kd, const double *A,
		const int *lead_dim_A, double *x, int *incx);
}

/** @endcond */

/**
 * @defgroup grp_helpers General helpers
 */

/**
 * @defgroup grp_timing Time measurement
 *
 * @ingroup grp_helpers
 */

/**
 * @ingroup grp_timing
 *
 * Returns the number of microseconds from the Epoch
 *
 * @return The number of microseconds from the Epoch
 */
uint_fast64_t get_time();

/**
 * @defgroup grp_random_wrap Wrappers for pseudo-random number generation
 *
 * @ingroup grp_helpers
 *
 * These wrappers make it easy to replace the pseudo-random number generator
 * (PRNG) globally -- just modify these wrappers to use the PRNG of your choice.
 * Originally functions from @ref grp_random are used.
 *
 */

/**
 * @ingroup grp_random_wrap Wrappers for pseudo-random number generation
 *
 * Returns a random number from /dev/urandom, that can be used to seed the
 * random number generator
 *
 * return A random number from /dev/urandom
 */
uint_fast64_t get_random_seed();

/**
 * @ingroup grp_random_wrap
 *
 * Initializes the pseudo-random number generator. Seeds the uniform number
 * generator with the provided value and makes sure the necessary data for the
 * ZIGNOR version of the Ziggurat algorithm is precomputed.
 *
 * @see seed_xorshift1024
 * @see precompute_zig_nor_data
 */
void init_rand_generator(uint_fast64_t seed);

/**
 * @ingroup grp_random_wrap
 *
 * Generates a random 64-bit unsigned integer number with uniform distribution.
 * The definition is in the header to make the function inline for performance.
 *
 * @return A random 64-bit unsigned integer number with uniform distribution.
 */
inline uint_fast64_t get_uni_rand_i() {
	return get_uni_rand_xorshift1024star_i();
}

/**
 * @ingroup grp_random_wrap
 *
 * Generates a random number with uniform distribution in (0, 1).
 * The definition is in the header to make the function inline for performance.
 *
 * @return A random number with uniform distribution in (0, 1).
 */
inline double get_uni_rand_d() {
	return get_uni_rand_xorshift1024star_d();
}

/**
 * @ingroup grp_random_wrap
 *
 * Generates a random number from normal distribution with 0 mean and
 * 1 std. dev.
 * The definition is in the header to make the function inline for performance.
 *
 * @return A random number from normal distribution with 0 mean and
 * 1 std. dev.
 */
inline double get_normal_rand_d() {
	return get_normal_rand_zignor_d();
}

/**
 * @defgroup grp_lin_alg Linear algebra helpers of BLAS/LAPACK routines
 *
 * @ingroup grp_helpers
 */

/**
 * @ingroup grp_lin_alg
 *
 * Performs the standard DAXPY operation, i.e. @f${\bf v_2} \leftarrow
 * \alpha{\bf v_1} + {\bf v_2}@f$.
 *
 * @param[in] alpha The value of the scalar weigth @f$\alpha@f$
 * @param[in] vec1 The vector @f${\bf v_1}@f$.
 * @param[inout] vec2 The vector @f${\bf v_2}@f$. On exit the result of the
 * operation.
 */
void add_vecs(double alpha, const std::vector<double> &vec1,
		std::vector<double> &vec2);

/**
 * @ingroup grp_lin_alg
 *
 * Computes the dot product of two vectors.
 *
 * @param[in] vec1 One of the vectors
 * @param[in] vec2 The other vector
 * @return The dot product of the two inputs
 */
double get_dot_prod(const std::vector<double> &vec1,
		const std::vector<double> &vec2);

/**
 * @ingroup grp_lin_alg
 *
 * Computes Cholesky factorization of the provided sparse matrix. The matrix is
 * assumed to be given in symmetric lower triangular band representation with
 * the given number of sub-diagonals.
 *
 * @param[in] num_subdiag The number of sub-diagonals of the matrix (not
 * including the main diagonal)
 * @param[inout] mat The lower triangular band storage of the matrix; on
 * return the lower triangular Cholesky factor.
 * @return 0 on successful execution, a non-zero value on error.
 */
int compute_cholesky_band(int num_subdiag, std::vector<double> &mat);

/**
 * @ingroup grp_lin_alg
 *
 * Solves a linear system using a band representation of a lower triangular
 * matrix.
 *
 * @param[in] trans If `'N'` the function solves @f$L {\bf x} = {\bf b}@f$,
 * if `'T'` it solves @f$L^T {\bf x} = {\bf b}@f$
 * @param[in] num_subdiag The number of sub-diagonals of the matrix (not
 * including the main diagonal)
 * @param[in] l_mat The lower triangular Cholesky factor of the matrix of the
 * linear system
 * @param[inout] vec The right hand side vector. On exit the result of the solve
 * @return 0 on successful execution, a non-zero value on error
 */
int solve_band_trian(char trans, int num_subdiag,
		const std::vector<double> &l_mat, std::vector<double> &vec);

/**
 * @ingroup grp_lin_alg
 *
 * Multiplies a vetor by a band lower triangular matrix.
 *
 * @param[in] trans If `'N'` the function solves @f$L {\bf x} = {\bf b}@f$,
 * if `'T'` it solves @f$L^T {\bf x} = {\bf b}@f$
 * @param[in] num_subdiag The number of sub-diagonals of the matrix (not
 * including the main diagonal)
 * @param [in] l_mat The lower triangular matrix in band storage
 * @param [inout] vec The vector to multiply by the matrix. On exit the result of the
 * multiplication.
 */
void mult_band_trian_mat_vec(char trans, int num_subdiag,
		const std::vector<double> &l_mat, std::vector<double> &vec);

/**
 * @ingroup grp_lin_alg
 *
 * Multiplies a vector by a full storage matrix (column-major storage) in the
 * form defined by the BLAS DGEMV function.
 *
 * @param[in] trans If `'N'` the function computes
 * @f${\bf y} = \alpha A{\bf x} + \beta {\bf y}@f$
 * if `'T'` the formula is @f${\bf y} = \alpha A^T{\bf x} + \beta {\bf y}@f$
 * @param[in] alpha The @f$\alpha@f$ multiplier in the above formulae
 * @param[in] A The matrix in full, column-major storage
 * @param[in] x The vector to multiply by the matrix
 * @param[in] beta The @f$\beta@f$ multiplier in the above formulae
 * @param[inout] y The input-output vector
 */
void mult_full_mat_vec(char trans, double alpha, const std::vector<double> &A,
		const std::vector<double> &x, double beta, std::vector<double> &y);

/**
 * @defgroup grp_sim_help Simulation helpers
 *
 * @ingroup grp_helpers
 */

/**
 * @ingroup grp_sim_help
 *
 * Prints a usage information - instructions of how to run the binary.
 *
 * @param binary The command to execute the binary
 * @param options An array of characters that define the standard POSIX `getopt`
 * options
 * @param num_opts The number of options
 * @param num_opts_required The number of options required (first
 * @p num_opts_required from @p options have to be provided for each run, the
 * rest is optional)
 * @param options_msg Information to be printed about each option
 * @param options_example An example argument for each option (for an example
 * of how to run the binary) to be printed out
 */
void print_usage_info(std::string binary, char options[], int num_opts,
		int num_opts_required, const std::string options_msg[],
		const std::string options_example[]);

/**
 * Checks whether all the required command line options were provided. If not
 * prints information about the missing ones.
 *
 * @param options The array of characters defining each option.
 * @param num_opts_required The number of options that are required (it is
 * assumed required options are put in @ options before the optional ones)
 * @param options_found The entries define whether the respective option have
 * been found on command line
 * @param options_msg Messages for all the options that should be printed if
 * a option is not found on command line
 *
 * @return 'true' if all the required option were found on command line
 * (according to @p options_found), `false` if any required option is missing
 */
bool check_missing_opts(char options[], int num_opts_required,
		const std::vector<bool> &options_found, std::string options_msg[]);

/**
 * @ingroup grp_sim_help
 *
 * Return a suffix that should be used for file name for all output.
 *
 * @param num_dropped_bp The number of base pairs dropped from each end
 * @param use_jac A flag indicating if the Jacobian was used for the run or not
 * @param random_seed The seed used for the random number generator
 * @param intrinsic A flag indicating whether the suffix for an intrinsic data
 * file (`true`) or simulation file (`false`) should be returned
 *
 * @return A suffix that should be used for file name for all output.
 */
std::string get_file_suffix(int num_dropped_bp, bool use_jac, long random_seed,
		bool intrinsic);

/**
 * @ingroup grp_sim_help
 *
 * Generates a random move according to the multivariate Gaussian distribution
 * with the given mean and inverse covariance given by the Cholesky factor
 * @f$L@f$ in triangular band storage
 *
 * @param[in] chol_L_band The band storage of the Cholesky factor @f$L@f$ of the
 * inverse covariance
 * @param[in] mean The given mean of the distribution
 * @param[out] move Storage for the generated move
 *
 * @return 0 on successful execution, a non-zero value on error
 */
int generate_random_move_band(const std::vector<double> &chol_L_band,
		const std::vector<double> &mean, std::vector<double> &move);

/**
 * @ingroup grp_sim_help
 *
 * Generates a random move according to the multivariate Gaussian distribution
 * with the given mean and a full matrix @f$PD^{-\frac{1}{2}}@f$ that comes from
 * the spectral decomposition of the inverse covariance matrix.
 *
 * @param[in] p_d_sqrt_inv The full storage of the @f$PD^{-\frac{1}{2}}@f$ from
 * the spectral decomposition of the inverse covariance matrix
 * @param[in] mean The given mean of the distribution
 * @param[out] move Storage for the generated move
 *
 * @return 0 on successful execution, a non-zero value on error
 */
int generate_random_move_spect(const std::vector<double> &p_d_sqrt_inv,
		const std::vector<double> &mean, std::vector<double> &move);

/**
 * @defgroup grp_comp_expect Computation of expectations
 *
 * @ingroup grp_helpers
 */

/**
 * @ingroup grp_comp_expect
 *
 * Computes the the @f${\bf t}^{[0]}_i@f$ tangents at each base pair @f$i@f$
 * defined as base pair normals (the third column of the matrix representations)
 * of base pair frames @f${\bf R}_i@f$ provided as quaternions.
 *
 * The function ignores @p num_dropped_bp base pairs at either end to avoid
 * end effects. The base pair index @f$i@f$ in the formulae is defined so that
 * the @p num_dropped_bp'th base pair (the first one that is not dropped) has
 * index 0.
 *
 * @param[in] num_dropped_bp Number of base pairs to ignore at each end.
 * @param[in] R_bp The orientations of base pairs provided as quaternions.
 * @param[out] t0 Storage for the computed tangents.
 *
 * @return 0 on successful execution, a non-zero value on error
 */
int compute_t0(int num_dropped_bp,
		const std::vector<algebra3d::QuaternionD> &R_bp,
		std::vector<algebra3d::Vector3D> &t0);

/**
 * @ingroup grp_comp_expect
 *
 * Computes the @f${\bf t}^{[1]}@f$ tangents at each base pair @f$i@f$, defined
 * as unit chords between the provided points:
 * @f$\displaystyle\frac{{\bf r}_{i+1} - {\bf r}_{i}} {|| {\bf r}_{i+1} - {\bf r}_{i} ||}@f$.
 *
 * The function ignores @p num_dropped_bp base pairs at either end to avoid
 * end effects. The base pair index @f$i@f$ in the formulae is defined so that
 * the @p num_dropped_bp'th base pair (the first one that is not dropped) has
 * index 0.
 *
 * @param[in] num_dropped_bp Number of base pairs to ignore at each end.
 * @param[in] r_bp The base pair positions that define the the chords.
 * @param[out] t1 Storage for the computed tangents.
 *
 * @return 0 on successful execution, a non-zero value on error
 *
 * @note @f${\bf t}^{[1]}@f$ could be included in @ref compute_tk, but it is
 * kept separate for efficiency.
 */
int compute_t1(int num_dropped_bp, const std::vector<algebra3d::Vector3D> &r_bp,
		std::vector<algebra3d::Vector3D> &t1);

/**
 * @ingroup grp_comp_expect
 *
 * Computes the generalized tangents @f${\bf t}^{[k]}_i@f$ at each base pair
 * @f$i@f$, for a given @f$k@f$. The tangent @f${\bf t}^{[k]}_i@f$ is computed as
 * unit eigenvector (with positive projection on the chord
 * @f${\bf r}_{i+ \lceil k / 2 \rceil} - {\bf r}_{i - \lfloor k / 2 \rfloor}@f$)
 * corresponding to the largest eigenvalue of the local gyration matrix
 * @f$\displaystyle\sum_{j = i - \lfloor k / 2 \rfloor}^
 {i + \lceil k / 2 \rceil}
 ({\bf r}_{j}-{\bf c}_i^k) \otimes ({\bf r}_{j}-{\bf c}_i^k)@f$
 * where @f${\bf c}_i^k = \displaystyle\frac{1}{k+1}
 \displaystyle\sum_{j = i - \lfloor k / 2 \rfloor}^
 {i + \lceil k / 2 \rceil}{\bf r}_{j}@f$ is the geometrical centre of mass.
 *
 * The function ignores @p num_dropped_bp base pairs at either end to avoid
 * end effects. The base pair index @f$i@f$ in the formulae is defined so that
 * the @p num_dropped_bp'th base pair (the first one that is not dropped) has
 * index 0.
 *
 * @param[in] num_dropped_bp Number of base pairs to ignore at each end.
 * @param[in] r_bp The base pair positions that define the the generalized chords.
 * @param[in] k The "index" of the tangent (> 1)
 * @param[out] tk Storage for the computed tangents
 *
 * @return 0 on successful execution, a non-zero value on error
 */
int compute_tk(int num_dropped_bp, const std::vector<algebra3d::Vector3D> &r_bp,
		int k, std::vector<algebra3d::Vector3D> &tk);

/**
 * @ingroup grp_comp_expect
 *
 * Computes the (generalized) arclengths @f${\bf s}^{[k]}_i@f$ for a given
 * @f$k>0@f$ at each base pair @f$i@f$. Let @f$n@f$ be the number of base pairs
 * considered and @f$q@f$ the reminder of the division @f$n / k@f$. The arclength
 * @f${\bf s}^{[k]}@f$ for indices @f$i \in \{0, k, 2k, \ldots, n - q, n\}@f$
 * is computed as the a sum
 * @f${\bf s}_i^{[k]} = \displaystyle\sum_{j = 0}^{i/k - 1}
 \delta{\bf s}_{j \cdot k}^{[k]}@f$
 * of distances of base pairs separated by @f$k@f$ junctions
 * @f$\delta{\bf s}_p^{[k]} := ||{\bf r}_{p + k} - {\bf r}_{p}||@f$.
 * For index @f$n@f$ the last "step" is defined as
 * @f$\delta{\bf s}_n^{[k]} := ||{\bf r}_{n} - {\bf r}_{n-q}||@f$. Note that
 * @f${\bf s}^{[1]}@f$ is simply the "standard" arclenth of the polyline defined
 * by @f${\bf r}_i@f$
 *
 * The function ignores @p num_dropped_bp base pairs at either end to avoid
 * end effects. The base pair index @f$i@f$ in the formulae is defined so that
 * the @p num_dropped_bp'th base pair (the first one that is not dropped) has
 * index 0.
 *
 * @note This function only computes the value of @f${\bf s}^{[k]}@f$ for the
 * subset of indices mentioned above (no other elements of the argument @p ts
 * are accessed). For all other indices the values can be computed as linear
 * interpolation of the values of the closest 2 indices from the above
 * definition as in @ref interpolate_sk(). The interpolation is not computed
 * here for efficiency. Thanks to linearity of expected values this computation
 * can be left for when the data is used (e.g. written to a file).
 *
 * @note The computed values of @f${\bf s}^{[k]}@f$ are expressed in Angstroms
 * as these are the units used in the
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> model.
 *
 * @param[in] num_dropped_bp Number of base pairs to ignore at each end.
 * @param[in] r_bp The base pair positions that define the the generalized chords.
 * @param[in] k The "index" of the arclength (> 0)
 * @param[out] sk Storage for the computed arclengths
 *
 * @return 0 on successful execution, a non-zero value on error
 *
 * @see interpolate_sk
 */
int compute_sk(int num_dropped_bp, const std::vector<algebra3d::Vector3D>&r_bp,
		int k, std::vector<double> &sk);

/**
 * @ingroup grp_comp_expect
 *
 * For an array of values of arclength @f${\bf s}^{[k]}_i@f$ computed for indices
 * @f$i \in \{0, k, 2k, \ldots, n\}@f$ (e.g. an output from @ref compute_sk)
 * computes the values for all intermediate points using linear interpolation.
 *
 * @param[in] k The "index" of arclength
 * @param[inout] sk The storage for the arclengths with values at indices
 * @f$i \in \{0, k, 2k, \ldots, n - r, n\}@f$ already computed; on exit all
 * entries are filled with linear interpolation between the provided values.
 *
 * @return 0 on successful execution, a non-zero value on error
 *
 * @see compute_sk
 */
int interpolate_sk(int k, std::vector<double> &sk);

/**
 * @ingroup grp_comp_expect
 *
 * Computes data for the cord vectors @f${\bf R}_0^T({\bf r}_i - {\bf r}_0)@f$
 * for computing the Flory persistence vector at each base pair @f$i@f$.
 *
 * The function ignores @p num_dropped_bp base pairs at either end to avoid
 * end effects. The base pair index @f$i@f$ in the formulae is defined so that
 * the @p num_dropped_bp'th base pair (the first one that is not dropped) has
 * index 0.
 *
 * @note The computed Flory persistence vectors are expressed in Angstroms as
 * these are the units used in the
 * <a href="http://lcvmwww.epfl.ch/cgDNA">cgDNA</a> model.
 *
 * @param[in] num_dropped_bp Number of points to ignore at each end.
 * @param[in] r_bp The base pair positions that define the the chords.
 * @param[in] R_bp The base pair orientations; the function requires only
 * @f${\bf R}_0@f$ -- the orientation of the @p num_dropped_bp'th base pair.
 * @param[out] r Storage for the computed vectors
 *
 * @return 0 on successful execution, a non-zero value on error
 */
int compute_r(int num_dropped_bp, const std::vector<algebra3d::Vector3D> &r_bp,
		const std::vector<algebra3d::QuaternionD> &R_bp,
		std::vector<algebra3d::Vector3D> &r);

#endif /* HELPERS_H_ */

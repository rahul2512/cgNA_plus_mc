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

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>

#include <algebra3d/Vector3.h>
#include <algebra3d/Quaternion.h>
#include <algebra3d/Matrix3.h>

#include <helpers.h>

using namespace std;
using namespace algebra3d;

/*****************************************************************************/
uint_fast64_t get_time() {
	struct timeval time;
	uint_fast64_t t;

	gettimeofday(&time, NULL);
	t = time.tv_sec;
	t *= 1000000;
	t += time.tv_usec;

	return t;
}

/*****************************************************************************/
uint_fast64_t get_random_seed() {
	int64_t seed;
	int fd = open("/dev/urandom", O_RDONLY);

	if (fd == -1) {
		cerr << "Couldn't open /dev/urandom to read a read a random seed"
				<< endl;
		close(fd);
		_exit(1);
	}
	if (read(fd, &seed, sizeof(seed)) != sizeof(seed)) {
		cerr << "Couldn't read a random seed from /dev/urandom" << endl;
		close(fd);
		_exit(1);
	}

	close(fd);

	// Return only positive values
	// Doing this exactly that way makes it easy to handle such number
	// Using standard 'string to number' functions without errors
	return (seed < 0 ? -seed : seed);
}

/*****************************************************************************/
void init_rand_generator(uint_fast64_t seed) {
	// Seed the generator
	seed_xorshift1024(seed);
	// Make sure the ZIGNOR data is precomputed
	precompute_zig_nor_data();
}

/*****************************************************************************/
void print_usage_info(string binary, char options[], int num_opts,
		int num_opts_required, const string options_msg[],
		const string options_example[]) {
	cerr << "Usage:\n"
			"    " << binary << " options\n\n"
			"Where options (given in any order) are:\n";
	for (int i = 0; i < num_opts; ++i) {
		cerr << "    -" << options[i] << "  ";
		if (i < num_opts_required) {
			cerr << "required";
		} else {
			cerr << "optional";
		}
		cerr << "\n" << options_msg[i] << endl;
	}
	cerr << "Example:\n" << "    " << binary;
	for (int i = 0; i < num_opts; ++i) {
		if (options_example[i].length() > 0) {
			cerr << "  -" << options[i] << " " << options_example[i];
		}
	}
	cerr << endl;
}

/*****************************************************************************/
bool check_missing_opts(char options[], int num_opts_required,
		const vector<bool> &options_found, string options_msg[]) {
	for (int i = 0; i < num_opts_required; ++i) {
		if (!options_found[i]) {
			cerr << "Error!\nMissing command line options:\n";
			for (int j = 0; j < num_opts_required; ++j) {
				if (!options_found[j]) {
					cerr << "    -" << options[j] << "  required\n"
							<< options_msg[j] << endl;
				}
			}
			return false;
		}
	}
	return true;
}

/*****************************************************************************/
string get_file_suffix(int num_dropped_bp, bool use_jac, long random_seed,
		bool intrinsic) {
	stringstream str_stream;

	if (intrinsic) {
		str_stream << "intr_";
	}

	str_stream << num_dropped_bp;
	// Any of the remaining information does not affect the intrinsic shape data
	// hence is not included in the suffix
	if (!intrinsic) {
		str_stream << (use_jac ? "_j" : "_nj") << "_" << random_seed;
	}
	return str_stream.str();
}

/*****************************************************************************/
void  add_vecs(double alpha, const vector<double> &vec1, vector<double> &vec2) {
	int dim = vec1.size();
	int one = 1;

	// Compute vec2 = alpha * vec1 + vec2.
	daxpy_(&dim, &alpha, &(vec1[0]), &one, &(vec2[0]), &one);
}

/*****************************************************************************/
double get_dot_prod(const std::vector<double> &vec1,
		const std::vector<double> &vec2) {
	int dim = vec1.size();
	int one = 1;

	return ddot_(&dim, &(vec1[0]), &one, &(vec2[0]), &one);
}

/*****************************************************************************/
int compute_cholesky_band(int num_subdiag, vector<double> &mat) {
	int dim = mat.size() / (num_subdiag + 1);

	char up_lo = 'L';				// lower triangular

	int kd = num_subdiag;			// number of subdiagonas
									// (or superdiagonals if 1st param. is 'U')
	int lead_dim_A = kd + 1;		// leading dimension of the storage of the matrix

	int info;

	dpbtrf_(&up_lo, &dim, &kd, &(mat[0]), &lead_dim_A, &info);

	return info;
}

/*****************************************************************************/
int solve_band_trian(char trans, int num_subdiag, const vector<double> &l_mat,
		vector<double> &vec) {
	int dim = l_mat.size() / (num_subdiag + 1);

	char up_lo = 'L';				// lower triangular
	char diag = 'N';				// non-unit diagonal

	int kd = num_subdiag;			// number of subdiagonas
									// (or superdiagonals if 1st param. is 'U')
	int num_rhs = 1;				// number of right rand sides
	int lead_dim_A = kd + 1;		// leading dimension of the storage of l_mat
	int lead_dim_B = dim;			// leading dimension of the storage of vec

	int info;

	// Solve L * x = b.
	dtbtrs_(&up_lo, &trans, &diag, &dim, &kd, &num_rhs, &(l_mat[0]),
			&lead_dim_A, &(vec[0]), &lead_dim_B, &info);

	return 2 * info;
}

/*****************************************************************************/
void mult_band_trian_mat_vec(char trans, int num_subdiag,
		const std::vector<double> &l_mat, std::vector<double> &vec) {
	int dim = l_mat.size() / (num_subdiag + 1);

	int one = 1;
	char up_lo = 'L';				// lower triangular
	char unit_diag = 'N';			// lower triangular

	int kd = num_subdiag;			// number of subdiagonas
									// (or superdiagonals if 1st param. is 'U')
	int lead_dim_A = kd + 1;		// leading dimension of the storage of l_mat

	dtbmv_(&up_lo, &trans, &unit_diag, &dim, &kd, &(l_mat[0]), &lead_dim_A,
			&(vec[0]), &one);
}

/*****************************************************************************/
void mult_full_mat_vec(char trans, double alpha, const std::vector<double> &A,
		const std::vector<double> &x, double beta, std::vector<double> &y) {
	int m = x.size();
	int one = 1;

	dgemv_(&trans, &m, &m, &alpha, &(A[0]), &m, &(x[0]), &one, &beta, &(y[0]),
			&one);
}

/*****************************************************************************/
int generate_random_move_band(const vector<double> &chol_L_band,
		const vector<double> &mean, vector<double> &move) {
	int dim = mean.size();
	int num_subdiag = chol_L_band.size() / dim - 1;
	int result;

	// Check if the sizes of the banded matrix storage and the mean vector
	// and the number of sub-diagonals are consistent
	if (chol_L_band.size() % dim != 0) {
		return 1;
	}

	move.resize(dim);

	// Generate random numbers from the normal distribution
	for (int i = 0; i < dim; ++i) {
		move[i] = get_normal_rand_d();
	}

	// Solve the linear system L^{T}(x - mean) = move o obtain the correct
	// covariance:
	// move <- x - mean = L^{-T} * move
	result = solve_band_trian('T', num_subdiag, chol_L_band, move);
	if (result != 0) {
		return result;
	}

	// Add the shift: move <- x = move + mean
	for (int i = 0; i < dim; ++i) {
		move[i] += mean[i];
	}

	return 0;
}

/*****************************************************************************/
int generate_random_move_spect(const vector<double> &p_d_sqrt_inv,
		const vector<double> &mean, vector<double> &move) {
	int dim = mean.size();
	vector<double> std_move;

	// Check if the sizes of the P*D^{-1/2} matrix storage and the mean vector
	if ((int) p_d_sqrt_inv.size() != dim * dim) {
		return 1;
	}

	std_move.resize(dim);
	move.resize(dim);

	// Generate random numbers from the normal distribution
	for (int i = 0; i < dim; ++i) {
		std_move[i] = get_normal_rand_d();
	}

	// Multiply the move by the P*D^{-1/2} matrix to obtain the correct
	// covariance:
	// move <- x - mean = P*D^{-1/2} * move
	mult_full_mat_vec('N', 1.0, p_d_sqrt_inv, std_move, 0.0, move);

	// Add the shift: move <- x = move + mean
	for (int i = 0; i < dim; ++i) {
		move[i] += mean[i];
	}

	return 0;
}

/*****************************************************************************/
int compute_t0(int num_dropped_bp, const vector<QuaternionD> &R_bp,
		vector<Vector3D> &t0) {

	int ind;
	int num_bp = R_bp.size();
	int num_left = num_bp - 2 * num_dropped_bp;

	// Check whether the number of dropped base pairs is reasonable
	if (num_left < 1) {
		t0.resize(0);
		return 1;
	}

	t0.resize(num_left);

	// Compute the tangents
	for (int i = 0; i < (int) t0.size(); ++i) {
		ind = i + num_dropped_bp;
		t0[i] = ((Matrix3D) R_bp[ind]).getColumn(2);

	}
	return 0;
}

/*****************************************************************************/
int compute_t1(int num_dropped_bp, const vector<Vector3D> &r_bp,
		vector<Vector3D> &t1) {

	int ind;
	int num_bp = r_bp.size();
	int num_left = num_bp - 2 * num_dropped_bp;

	// Check whether the number of dropped base pairs is reasonable
	if (num_left < 1) {
		t1.resize(0);
		return 1;
	}

	// If no base pairs are dropped the chord from the last base pair on
	// cannot be computed
	if (num_dropped_bp == 0) {
		t1.resize(num_left - 1);
	} else {
		t1.resize(num_left);
	}

	// Compute the tangents
	for (int i = 0; i < (int) t1.size(); ++i) {
		ind = i + num_dropped_bp;
		t1[i] = (r_bp[ind + 1] - r_bp[ind]).getUnit();
	}

	return 0;
}

/*****************************************************************************/
int compute_tk(int num_dropped_bp, const vector<Vector3D> &r_bp, int k,
		vector<Vector3D> &tk) {
	Vector3D avg;
	Vector3D curr = Vector3D(0.0, 0.0, 1.0);
	Vector3D prev = Vector3D(0.0, 0.0, 0.0);
	double margin = 1.0e-10;
	Matrix3D cov;
	Vector3D diff;
	int num_left = r_bp.size() - 2 * num_dropped_bp;
	int back;
	int front;
	int start_r;
	int end_r;
	int num_pts;
	// Make sure anything can be done with the base pairs left
	if (k < num_left) {
		// Resize and the array if necessary
		tk.resize(num_left);
		for (int i = 0; i < num_left; ++i) {
			back = (k / 2);
			front = (0.5 * k + 0.5);

			start_r = max(0, num_dropped_bp + i - back);
			end_r = min((int) r_bp.size() - 1, num_dropped_bp + i + front);
			num_pts = end_r - start_r + 1;

			// Compute the average
			avg = r_bp[start_r];
			for (int k = start_r - num_dropped_bp - i + back + 1;
					num_dropped_bp + i - back + k <= end_r; ++k) {
				avg += r_bp[num_dropped_bp + i - back + k];
			}
			avg /= num_pts;

			// Compute covariance
			diff = r_bp[start_r] - avg;
			cov = diff.outer(diff);
			for (int k = start_r - num_dropped_bp - i + back + 1;
					num_dropped_bp + i - back + k <= end_r; ++k) {
				diff = r_bp[num_dropped_bp + i - back + k] - avg;
				// Outer porduct
				cov += diff.outer(diff);
			}
			cov /= num_pts;

			// Find principal component by power iterations
			// A reasonable approximation
			curr = r_bp[start_r] - r_bp[end_r];
			curr.normalize();
			while ((curr - prev).getNorm() > margin) {
				prev = curr;
				curr = cov * curr;
				curr.normalize();
			}
			// Sign correction is unnecessary due to the choice
			// of starting point for power iterations
			tk[i] = curr;
		}
	} else {
		return 1;
	}

	return 0;
}

/*****************************************************************************/
int compute_sk(int num_dropped_bp, const vector<Vector3D> &r_bp, int k,
		vector<double> &sk) {
	int num_left = r_bp.size() - 2 * num_dropped_bp;
	int i;

	// Make sure anything can be done with the base pairs left
	if (k < num_left) {
		// Resize and the array if necessary
		sk.resize(num_left, 0.0);
		sk[0] = 0.0;
		for (i = k; i < (int) sk.size(); i += k) {
			sk[i] =
					sk[i - k]
							+ (r_bp[num_dropped_bp + i]
									- r_bp[num_dropped_bp + i - k]).getNorm();
		}
		// In most cases arclngth for the last point is not assigned
		i -= k;
		sk[sk.size() - 1] = sk[i]
				+ (r_bp[num_dropped_bp + (sk.size() - 1)]
						- r_bp[num_dropped_bp + i]).getNorm();
	} else {
		return 1;
	}
	return 0;
}

/*****************************************************************************/
int interpolate_sk(int k, std::vector<double> &sk) {
	double offset;
	double delta;
	int i;
	int size = sk.size();
	// For each ith point
	for (i = k; i < size; i += k) {
		// Interpolate values between point (i - k) and i
		offset = sk[i - k];
		delta = (sk[i] - sk[i - k]) / k;
		for (int j = 1; j <= k; ++j) {
			sk[i - k + j] = offset + delta * j;
		}
	}
	// Interpolate the reminder that might be shorter than k
	// (the last point is always computed)
	offset = sk[i - k];
	delta = (sk[size - 1] - sk[i - k]) / (size - 1 - i + k);
	for (int j = 1; i - k + j < size; ++j) {
		sk[i - k + j] = offset + delta * j;
	}

	return 0;
}

/*****************************************************************************/
int compute_r(int num_dropped_bp, const vector<Vector3D> &r_bp,
		const vector<QuaternionD> &R_bp, vector<Vector3D> &r) {
	// Resize the array if necessary
	r.resize(r_bp.size() - 2 * num_dropped_bp);
	// Get the inverse of the orientation of the first bas pair
	Matrix3D R_first = R_bp[num_dropped_bp].inv();
	for (unsigned int i = 0; i < r.size(); ++i) {
		r[i] = R_first * (r_bp[i + num_dropped_bp] - r_bp[num_dropped_bp]);
	}

	return 0;
}

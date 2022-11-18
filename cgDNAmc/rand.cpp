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

#include <rand.h>
#include <cmath>

using namespace std;

uint_fast64_t XORSHIFT1024_STATE[16] = { 0xaaaaaaaaaaaaaaaaLL,
		0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL,
		0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL,
		0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL,
		0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL,
		0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL, 0xaaaaaaaaaaaaaaaaLL };
int XORSHIFT1024_INDEX;

/** The number of block of the ZIGNOR algorithm */
static const int ZIGNOR_C = 128;
/** The beginning of the right tail of the ZIGNOR algorithm */
static const double ZIGNOR_R = 3.442619855899;
/**
 * The value of the area of the 0th block - the tail:
 * V = (R * phi(R) + Pr(X>=R)) * sqrt(2\pi)
 * (see the Doornik paper).
 */
static const double ZIGNOR_V = 9.91256303526217e-3;
/**
 * An array to store the coordinates of the blocks of equal area.
 * This is precomputed for performance.
 */
static double s_adZigX[ZIGNOR_C + 1];
/**
 * An array to store the rations of subsequent coordinates of the blocks of
 * equal area.
 * This is precomputed for performance.
 */
static double s_adZigR[ZIGNOR_C];

/*****************************************************************************/
void seed_xorshift1024(uint_fast64_t seed) {
	// Make ruse each call of this function with the same value of seed
	// restarts exqctly the same sequence of pseudo-random numbers
	XORSHIFT1024_INDEX = 0;

	// Initialize the entire state with the seed
	for (int i = 0; i < 8; ++i) {
		// Make sure the state is not all bits 0 and not all bits 1
		// To do that every element of the state vector has one half filled with
		// alternationg 0s and 1s

		XORSHIFT1024_STATE[2 * i] =
				// Higher 4 bytes: higher 4 bytes of the seed
				(seed & 0xffffffff00000000LL)
				// Lower 4 bytes: alternating 0s and 1s
				| 0xaaaaaaaaLL;
		XORSHIFT1024_STATE[2 * i + 1] =
				// Higher 4 bytes: alternating 0s and 1s
				0xaaaaaaaa00000000LL
				// Lower 4 bytes: lower 4 bytes of the seed
				| (seed & 0x00000000ffffffffLL) ;
	}
}

/*****************************************************************************/
void precompute_zig_nor_data() {
	int i;
	// f(x) is the normal distribution density
	// Here we do not need to normalize
	double f = exp(-0.5 * ZIGNOR_R * ZIGNOR_R);
	// The lowest block's corner is is V / f(R)
	s_adZigX[0] = ZIGNOR_V / f;
	// The lowest block's corner is is V / f(R)
	s_adZigX[1] = ZIGNOR_R;
	// Last x is 0
	s_adZigX[ZIGNOR_C] = 0;
	for (i = 2; i < ZIGNOR_C; ++i) {
		s_adZigX[i] = sqrt(-2 * log(ZIGNOR_V / s_adZigX[i - 1] + f));
		f = exp(-0.5 * s_adZigX[i] * s_adZigX[i]);
	}
	for (i = 0; i < ZIGNOR_C; ++i) {
		s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
	}
}

/*****************************************************************************/
double get_normal_rand_zignor_d() {
	unsigned int i;
	double x;
	double y;
	double u;
	double f0;
	double f1;

	// This is a rejection method - try until a value is found
	while (true) {
		// Get a uniform random number i [-1, 1]
		u = 2 * get_uni_rand_xorshift1024star_d() - 1;
		// Get a unifom integer number
		// Get the index of a box chosen from a single selected byte
		i = get_uni_rand_xorshift1024star_i() & 0x7F;
		// If the number is trivially 'inside' the rectangular box
		// we have a result
		if (fabs(u) < s_adZigR[i]) {
			return u * s_adZigX[i];
		}
		// Otherwise if the level with the tail is selected - do the tail test
		if (i == 0) {
			do {
				x = log(get_uni_rand_xorshift1024star_d()) / ZIGNOR_R;
				y = log(get_uni_rand_xorshift1024star_d());
			} while (-2 * y < x * x);
			return (u < 0) ? x - ZIGNOR_R : ZIGNOR_R - x;
		}
		// Are we are inside the in the wedges?
		x = u * s_adZigX[i];
		f0 = exp(-0.5 * (s_adZigX[i] * s_adZigX[i] - x * x));
		f1 = exp(-0.5 * (s_adZigX[i + 1] * s_adZigX[i + 1] - x * x));
		if (f1 + get_uni_rand_xorshift1024star_d() * (f0 - f1) < 1.0) {
			return x;
		}
	}
}

/*****************************************************************************/
double get_normal_rand_polar_d() {
	// A generated value stored between calls
	// This will store the generated value that is not returned in the
	// generation call
	static double stored = 0.0;
	double v1;
	double v2;
	double radiusSq;
	double frac;

	// Every other call draws two numbers, each from a normal distribution
	if (stored == 0.0) {
		// Find coordinates within a unit circle (uniformly distributed)
		do {
			v1 = 2.0 * get_uni_rand_xorshift1024star_d() - 1.0;
			v2 = 2.0 * get_uni_rand_xorshift1024star_d() - 1.0;
			radiusSq = v1 * v1 + v2 * v2;
		} while (radiusSq >= 1.0 || radiusSq == 0.0);
		// Use the Box-Muller transformation and return one of the two possible
		// store one of them for next call and return the other one
		frac = sqrt(-2.0 * log(radiusSq) / radiusSq);
		// Remember one of the values for next call...
		stored = v1 * frac;
		// ... and return the other one
		return v2 * frac;
	}
	// The other call simply returns a value already computed in the
	// preceding call
	else {
		// Use the value computed in previous call...
		v1 = stored;
		// .. and make sure generation is done in the next call
		stored = 0.0;
		return v1;
	}
}

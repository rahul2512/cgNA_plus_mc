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

/*
 * This code uses the ZIGNOR implementation of the Ziggurat algorithm.
 * Below is the copyright statement of that implementation.
 */
/*==========================================================================
 *  This code is Copyright (C) 2005, Jurgen A. Doornik.
 *  Permission to use this code for non-commercial purposes
 *  is hereby given, provided proper reference is made to:
 *		Doornik, J.A. (2005), "An Improved Ziggurat Method to Generate Normal
 *          Random Samples", mimeo, Nuffield College, University of Oxford,
 *			and www.doornik.com/research.
 *		or the published version when available.
 *	This reference is still required when using modified versions of the code.
 *  This notice should be maintained in modified versions of the code.
 *	No warranty is given regarding the correctness of this code.
 *==========================================================================*/

/**
 * @file
 * Pseudo-random number generation functions
 */

#ifndef RAND_H_
#define RAND_H_

#include <stdint.h>

/**
 * @defgroup grp_random Pseudo-random number generation
 *
 * These wrappers makes it easy to change the random number generator globally.
 * Originally the functions from @ref grp_random are used. Everything is defined
 * in global scope to avoid any potential overhead of encapsulation.
 */

/**
 * @ingroup grp_random
 *
 * A global state required by the xorshift1024* random number generator.
 * It is initialized by default with altering zeros and ones.
 *
 * @note A state of all bits 0 and all bits 1 are a fixed points of the algorithm
 * -- this needs to be avoided.
 *
 * see @ref seed_xorshift1024()
 */
extern uint_fast64_t XORSHIFT1024_STATE[16];
/**
 * @ingroup grp_random
 *
 * A global index required by the xorshift1024* random number generator
 */
extern int XORSHIFT1024_INDEX;

/**
 * @ingroup grp_random
 *
 * Seeds the xorshift1024* uniform random number generator with the provided
 * seed. Resets the @ref XORSHIFT1024_INDEX to 0 and fills in the
 * @ref XORSHIFT1024_STATE array. The procedure that fills in the state array
 * can be seen as setting every element of the array to the provided value of
 * @p seed and subsequently replacing the higher or lower 4 bytes of each
 * element with alternating 0 and 1 bits (see the code for details).
 *
 * @note Here the seed can be all 0 or all 1 bits. The procedure described
 * above makes sure that the state is not all bits 0 and not all bits 1 (which
 * are fixed points of the XORSHIFT algorithm). It also ensures the state is
 * different for every value of @p seed.
 *
 * @param seed The seed to use.
 *
 * @see @ref XORSHIFT1024_STATE
 * @see @ref XORSHIFT1024_INDEX
 */
void seed_xorshift1024(uint_fast64_t seed);

/**
 * @ingroup grp_random
 *
 * Precomputes the data necessary for the ZIGNOR algorithm.
 *
 * @note This initialization is necessary for the Ziggurat method to work.
 * Please always run this function before generating any pseudo-random numbers
 */
void precompute_zig_nor_data();

/**
 * @ingroup grp_random
 *
 * Generates a random 64-bit unsigned integer number with uniform distribution
 * using the xorshift1024* algorithm.
 * The definition is in the header to make the function inline for performance.
 *
 * @return A random 64-bit unsigned integer number with uniform distribution.
 *
 * @see Vigna, Sebastiano, "An experimental exploration of Marsaglia's
 * xorshift generators, scrambled", 2014, arXiv:1402.6246
 *
 */
inline uint_fast64_t get_uni_rand_xorshift1024star_i() {
	uint_fast64_t s0 = XORSHIFT1024_STATE[XORSHIFT1024_INDEX];
	uint_fast64_t s1 = XORSHIFT1024_STATE[XORSHIFT1024_INDEX = (XORSHIFT1024_INDEX
			+ 1) & 15];
	s1 ^= s1 << 31; // a
	s1 ^= s1 >> 11; // b
	s0 ^= s0 >> 30; // c
	return (XORSHIFT1024_STATE[XORSHIFT1024_INDEX] = s0 ^ s1)
			* 1181783497276652981LL;
}

/**
 * @ingroup grp_random
 *
 * Generates a random number in [0, 1) with uniform distribution using the
 * xorshift1024* algorithm by Vigna.
 *
 * @return A random 64-bit unsigned integer number with uniform distribution.
 *
 * @note The definition is in the header to make the function inline for
 * performance.
 */
inline double get_uni_rand_xorshift1024star_d() {
	// Divide by 2^64
	return 5.42101086242752217E-20 * get_uni_rand_xorshift1024star_i();
}

/**
 * @ingroup grp_random
 *
 * Generates a random  number from Gaussian distribution with 0 mean and
 * 1 std. dev. Uses the ZIGNOR variant of the Marsaglia ziggurat method,
 * proposed by Doornik. This is the faster of the two exact methods.
 *
 * @return A random  number from Gaussian distribution with 0 mean and
 * 1 std. dev.
 *
 * @ see Numerical Recipes in C, The Art of Scientiﬁc Computing
 */
double get_normal_rand_zignor_d();

/**
 * @ingroup grp_random
 *
 * Generates a random  number from Gaussian distribution with 0 mean and
 * 1 std. dev. Uses the Marsaglia polar method (a version of the Box-Muller
 * transformation) for two uniformly distributed coordinates of a point within
 * a unit circle.
 *
 * @return A random  number from Gaussian distribution with 0 mean and
 * 1 std. dev.
 *
 * @ see Numerical Recipes in C, The Art of Scientiﬁc Computing
 */
double get_normal_rand_polar_d();

#endif /* RAND_H_ */

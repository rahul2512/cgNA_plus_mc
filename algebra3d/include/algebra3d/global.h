/*
 * Copyright 2015 Jaroslaw Glowacki
 * jarek (dot) glowacki (at) gmail (dot) com
 *
 * This file is part of algebra3d.
 *
 * algebra3d is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * algebra3d is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with algebra3d.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file global.h
 *
 * Global functions and declarations.
 */

#ifndef _ALGEBRA3D_GLOBAL_H_
#define _ALGEBRA3D_GLOBAL_H_

#include<limits>
#include<cmath>
#include<ostream>
#include<istream>

/**
 * The floating type to use (float or double)
 */

/**
 * Default precision of comparison floating point type given as parameter.
 * It is 1.0-12 for doubles and 1.0e-05 for floats.
 */
#define _PRECISION(FP_TYPE) ( \
	(sizeof(FP_TYPE) == sizeof(float)) * 1.0e-5 \
	+ (sizeof(FP_TYPE) > sizeof(float)) * 1.0e-12\
	)
/**
 * Default precision of comparison for different floating point types.
 * Macro used within templates where fpType is the template parameter.
 */
#define PRECISION _PRECISION(fpType)
/**
 * Default precision for double.
 */
#define PRECISION_D _PRECISION(double)
/**
 * Default precision for float.
 */
#define PRECISION_F _PRECISION(float)

/**
 * Representation of infinity in floating point type given as parameter.
 */
#define _INF(FP_TYPE) (std::numeric_limits<FP_TYPE>::infinity())
/**
 * Representation of infinity in different floating point types.
 * Macro used within templates where fpType is the template parameter.
 */
#define INF _INF(fpType)
/**
 * Representation of infinity in for double.
 */
#define INF_D _INF(double)
/**
 * Representation of infinity in for float.
 */
#define INF_F _INF(float)

/**
 * Representation of NaN in different floating point types (INF - INF).
 * Macro used within templates where fpType is the template parameter.
 */
#define N_A_N (INF - INF)
/**
 * Representation of NaN in for double (INF_D - INF_D).
 */
#define N_A_N_D (INF_D - INF_D)
/**
 * Representation of NaN in for float (INF_F - INF_F).
 */
#define N_A_N_F (INF_F - INF_F)

/**
 * @namespace algebra3d
 * @brief A namespace for the whole library.
 */
namespace algebra3d {

/**
 * Returns the sign of the given number.
 *
 * @param[in] val The number to return the sign for.
 * @return The sign of the given number.
 *
 * @see signStar
 */
template<class fpType>
inline int sign(const fpType &val) {
	// A branchless version
	return (val > 0) - (val < 0);
}

/**
 * Returns the sign for non-negative numbers and 1 for 0.
 *
 * @param[in] val The number to return the sgn* for.
 * @return The sign for non-negative numbers and 1 for 0.
 *
 * @see sign
 */
template<class fpType>
inline int signStar(const fpType &val) {
	// A branchless version
	return (val >= 0) - (val < 0);
}

/**
 * Checks if the provided number is Inf.
 *
 * Two function pointers have been added to make it easier to use.
 *
 * @param[in] number The number to check.
 * @return `true` if the numbers is Inf, `false`
 * otherwise.
 */
template<class fpType>
inline bool isInf(const fpType &number) {
	return number == std::numeric_limits<fpType>::infinity()
			|| -number == std::numeric_limits<fpType>::infinity();
}

/**
 * Checks if the provided number is NaN.
 *
 * Two function pointers have been added to make it easier to use.
 *
 * @param[in] number The number to check.
 * @return `true` if the numbers is NaN, `false`
 * otherwise.
 */
template<class fpType>
inline bool isNaN(const fpType &number) {
	return number != number;
}

/**
 * Compares the two fpType numbers up to the precision. Returns 0 if the
 * numbers are equal up to precision, 1 if the first number is greater, -1
 * if the second number is greater, 2 if one of the numbers is `NaN`.
 *
 * Two function pointers have been added to make it easier to use.
 *
 * @param[in] d1 A number to compare.
 * @param[in] d2 The other number to compare.
 * @param[in] precision The precision to use; by default the global
 * PRECISION is used.
 * @return `0` if the numbers are equal up to PRECISION,
 * `1` if the first one is greater, `-1` if the second
 * number is greater, `2` if one of the numbers is `NaN`.
 *
 * @see compareD
 * @see compareF
 */
template<class fpType>
inline int compare(const fpType &d1, const fpType&d2, const fpType &precision =
PRECISION) {
	// A branchless version
	return (isNaN<fpType>(d1) || isNaN<fpType>(d2)) * 2 + (d1 - d2 > precision)
			- (d2 - d1 > precision);
}

/**
 * "Normalizes" the angle i.e. transforms it to (-PI, PI].
 *
 * Two function pointers have been added to make it easier to use.
 *
 * @param[in] angle The angle to "normalize".
 * @return The "normalized" angle.
 *
 * @see normalizeAngleD
 * @see normalizeAngleF
 */
template<class fpType>
fpType normalizeAngle(const fpType &angle);

/**
 * Computes a normalized angle for the given values of sine and cosine.
 *
 * @param sinAlpha sine of the angle to calculate
 * @param cosAlpha cosine of the angle to calculate.
 * @return The computed normalized angle.
 */
template<class fpType>
inline fpType getAngleFromSinCose(fpType sinAlpha, fpType cosAlpha) {
	return acos(cosAlpha) * ((sinAlpha > 0.0) - (sinAlpha < 0.0));
}

} /* namespace algebra3d */

#endif /* _ALGEBRA3D_GLOBAL_H_ */

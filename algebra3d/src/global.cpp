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

/*
 * Definitions of the global functions.
 */
#include <algebra3d/global.h>

namespace algebra3d {

/****************************************************************************/
template<class fpType>
fpType normalizeAngle(const fpType &angle) {
	int sign = compare<fpType>(angle, 0.0);
	fpType angleNorm = sign * angle;
	int div = (int) (angleNorm / M_PI * 0.5);
	angleNorm -= (div << 1) * M_PI;
	angleNorm -= 2 * M_PI * (angleNorm > M_PI);

	// The special case of exact equality to both pi and -pi should give pi
	int isPi = compare<fpType>(angleNorm, M_PI) == 0;
	sign = isPi + !isPi * sign;

	return sign * angleNorm;
}

// Define the specializations that will be usable
template
double normalizeAngle<double>(const double &angle);
template
float normalizeAngle<float>(const float &angle);

} /* namespace algebra3d */

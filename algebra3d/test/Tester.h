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
 * Header file for all tests. Defines the precision depending on compilation
 * parameters.
 */
#ifndef _TESTER_H_
#define _TESTER_H_

#ifdef FPTYPE_FLOAT
typedef float fpType;
#else
typedef double fpType;
#endif

// A flag indicating if a "heavy" test should be done
// This is a test scanning through different axes and a range of angles
extern bool DO_HEAVY_TESTS;

#endif // _TESTER_H_

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
 * A class to test global functions
 */

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <algebra3d/global.h>
#include <iostream>

#include "Tester.h"

using namespace std;
using namespace algebra3d;

class GlobalTest: public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE (GlobalTest);
	CPPUNIT_TEST (testSign);
	CPPUNIT_TEST (testCompare);
	CPPUNIT_TEST (testInf);
	CPPUNIT_TEST (testNaN);
	CPPUNIT_TEST (testAngleNormalization);
	CPPUNIT_TEST (testGetAngleSC);CPPUNIT_TEST_SUITE_END();

public:
	void setUp() {
	}

	void tearDown() {
	}

	void testSign() {
		CPPUNIT_ASSERT(sign<fpType>(1.2345) == 1);
		CPPUNIT_ASSERT(sign<fpType>(-98.765) == -1);
		CPPUNIT_ASSERT(sign<fpType>(0.0) == 0);
		CPPUNIT_ASSERT(sign<fpType>(-1.0e-26) == -1);
		CPPUNIT_ASSERT(sign<fpType>(1.0e-26) == 1);
		CPPUNIT_ASSERT(sign<fpType>(N_A_N) == 0);
		CPPUNIT_ASSERT(sign<fpType>(-INF) == -1);
		CPPUNIT_ASSERT(sign<fpType>(INF) == 1);
	}

	void testCompare() {
		CPPUNIT_ASSERT(compare<fpType>(0.0, 0.0) == 0);
		CPPUNIT_ASSERT(compare<fpType>(1.0, 1.0) == 0);
		CPPUNIT_ASSERT(compare<fpType>(0.5, 1.0) == -1);
		CPPUNIT_ASSERT(compare<fpType>(0.5, 0.3) == 1);
		CPPUNIT_ASSERT(compare<fpType>(PRECISION, 0.0) == 0);
		CPPUNIT_ASSERT(compare<fpType>(PRECISION, PRECISION) == 0);
		CPPUNIT_ASSERT(compare<fpType>(1.1 * PRECISION, -9.0 * PRECISION) == 1);
		CPPUNIT_ASSERT(
				compare<fpType>(-0.9 * PRECISION, 1.1 * PRECISION) == -1);

		CPPUNIT_ASSERT(compare<fpType>(INF, -INF) == 1);
		CPPUNIT_ASSERT(compare<fpType>(-INF, INF) == -1);
		CPPUNIT_ASSERT(compare<fpType>(INF, INF) == 0);
		CPPUNIT_ASSERT(compare<fpType>(-INF, -INF) == 0);

		CPPUNIT_ASSERT(compare<fpType>(N_A_N, -INF) == 2);
		CPPUNIT_ASSERT(compare<fpType>(INF, N_A_N) == 2);
		CPPUNIT_ASSERT(compare<fpType>(N_A_N, 0.5) == 2);
		CPPUNIT_ASSERT(compare<fpType>(-4.5, N_A_N) == 2);
	}

	void testInf() {
		CPPUNIT_ASSERT(isInf(INF));
		CPPUNIT_ASSERT(isInf(-INF));
		CPPUNIT_ASSERT(!isInf(1.0e16));
		CPPUNIT_ASSERT(!isInf(-1.0e16));
		CPPUNIT_ASSERT(!isInf(4.5));
		CPPUNIT_ASSERT(!isInf(-1.4));
	}

	void testNaN() {
		CPPUNIT_ASSERT(!isNaN(INF));
		CPPUNIT_ASSERT(!isNaN(-INF));
		CPPUNIT_ASSERT(!isNaN(1.0e16));
		CPPUNIT_ASSERT(!isNaN(-1.0e16));
		CPPUNIT_ASSERT(!isNaN(4.5));
		CPPUNIT_ASSERT(!isNaN(-1.4));
		CPPUNIT_ASSERT(isNaN(INF - INF));
	}

	void testAngleNormalization() {
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(0.0), 0.0) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(0.1), 0.1) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(-0.1), -0.1) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(2.0 * M_PI), 0.0) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(7.0 * M_PI), M_PI) == 0);
		CPPUNIT_ASSERT(
				compare<fpType>(normalizeAngle(1.7 * M_PI), -0.3 * M_PI) == 0);
		CPPUNIT_ASSERT(
				compare<fpType>(normalizeAngle(-1.6 * M_PI), 0.4 * M_PI) == 0);

		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(-M_PI + PRECISION), -M_PI + PRECISION) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(-M_PI - 0.5 * PRECISION), M_PI - 0.5 * PRECISION) == 0);

		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(-5.0 * M_PI), M_PI) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(-M_PI), M_PI) == 0);
		CPPUNIT_ASSERT(compare<fpType>(normalizeAngle(-8.0 * M_PI), 0.0) == 0);
		CPPUNIT_ASSERT(
				compare<fpType>(normalizeAngle(0.25 * M_PI), 0.25 * M_PI) == 0);
		CPPUNIT_ASSERT(
				compare<fpType>(normalizeAngle(-0.25 * M_PI - 4 * M_PI), -0.25 * M_PI) == 0);
		CPPUNIT_ASSERT(
				compare<fpType>(normalizeAngle(0.5 * M_PI + 3 * M_PI), -0.5 * M_PI) == 0);

		for (fpType angle = -10.0 * M_PI; angle < 10.0 * M_PI; angle += 0.1) {
			fpType normalized = normalizeAngle(angle);
			fpType sumOverPI = (abs(angle) + abs(normalized)) / M_PI;
			fpType diffOverPI = (abs(angle) - abs(normalized)) / M_PI;
			CPPUNIT_ASSERT(
					compare<fpType>(diffOverPI, round(diffOverPI)) == 0 || //
					compare<fpType>(sumOverPI, round(sumOverPI)) == 0);
		}
	}

	void testGetAngleSC() {
		for (fpType angle = -10.0 * M_PI; angle < 10.0 * M_PI; angle += 0.1) {
			fpType sinAlpha = sin(angle);
			fpType cosAlpha = cos(angle);
			fpType normalized = normalizeAngle(angle);
			CPPUNIT_ASSERT(
					compare<fpType>(normalized, getAngleFromSinCose(sinAlpha, cosAlpha)) == 0);
		}
	}

} ;

CPPUNIT_TEST_SUITE_REGISTRATION (GlobalTest);

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
 * A class to test vector class
 */

#include<cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include<cmath>

#include<algebra3d/global.h>
#include<algebra3d/Vector3.h>
#include<iostream>
#include<string>
#include<sstream>

#include "Tester.h"

using namespace std;
using namespace algebra3d;

class VectorTest: public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE (VectorTest);
	CPPUNIT_TEST (testInitEquality);
	CPPUNIT_TEST (testReferences);
	CPPUNIT_TEST (testInitDataAccess);
	CPPUNIT_TEST (testDataAccess);
	CPPUNIT_TEST (testDataAccessEquality);
	CPPUNIT_TEST (testEpsEquality);
	CPPUNIT_TEST (testIsZero);
	CPPUNIT_TEST (testAxisAng);
	CPPUNIT_TEST (testCopying);
	CPPUNIT_TEST (testCopyConstructor);
	CPPUNIT_TEST (testScale);
	CPPUNIT_TEST (testSetGet);
	CPPUNIT_TEST (testCross);
	CPPUNIT_TEST (testDot);
	CPPUNIT_TEST (testNorm);
	CPPUNIT_TEST (testAddSubt);
	CPPUNIT_TEST (testInput);CPPUNIT_TEST_SUITE_END();

public:
	void setUp() {
	}

	void tearDown() {
	}

	void testInitEquality() {
		Vector3<fpType> v1;
		Vector3<fpType> v2;
		Vector3<fpType> v3;

		// Initially matrices should be identical
		CPPUNIT_ASSERT(v1 == v1);
		CPPUNIT_ASSERT(!(v1 != v1));

		CPPUNIT_ASSERT(v1 == v2);
		CPPUNIT_ASSERT(!(v1 != v2));

		CPPUNIT_ASSERT(v2 == v3);
		CPPUNIT_ASSERT(!(v2 != v3));

	}

	void testReferences() {
		Vector3<fpType> v1;
		Vector3<fpType> v2(3.0, 4.0, 5.0);
		Vector3<fpType> v3(v2);

		// Initially matrices should be identical
		CPPUNIT_ASSERT(v1.x == 0.0);
		CPPUNIT_ASSERT(v1.y == 0.0);
		CPPUNIT_ASSERT(v1.z == 0.0);

		CPPUNIT_ASSERT(v2.x == 3.0);
		v2.x *= 2.0;
		CPPUNIT_ASSERT(v2.x == v2[0]);
		v2[0] *= 0.5;
		CPPUNIT_ASSERT(v2.x == v2[0]);

		CPPUNIT_ASSERT(v3.x == 3.0);
		v3.x *= 2.0;
		CPPUNIT_ASSERT(v3.x == v3[0]);
		v3[0] *= 0.5;
		CPPUNIT_ASSERT(v3.x == v3[0]);

		CPPUNIT_ASSERT(v2.y == 4.0);
		v2.y *= 2.0;
		CPPUNIT_ASSERT(v2.y == v2[1]);
		v2[1] *= 0.5;
		CPPUNIT_ASSERT(v2.y == v2[1]);

		CPPUNIT_ASSERT(v2.z == 5.0);
		v2.z *= 2.0;
		CPPUNIT_ASSERT(v2.z == v2[2]);
		v2[2] *= 0.5;
		CPPUNIT_ASSERT(v2.z == v2[2]);
	}

	void testInitDataAccess() {
		fpType a1[3] = { 1.0, -2.0, 3.0 };

		Vector3<fpType> v1;
		Vector3<fpType> v2(a1);
		Vector3<fpType> v3(a1[0], a1[1], a1[2]);
		Vector3<fpType> v4(v2);
		Vector3<fpType> v5(INF, INF, INF);
		Vector3<fpType> v6(NAN, NAN, NAN);

		for (int ind = 0; ind < 3; ++ind) {
			CPPUNIT_ASSERT(v1[ind] == 0.0);
			CPPUNIT_ASSERT(v2[ind] == a1[ind]);
			CPPUNIT_ASSERT(v3[ind] == a1[ind]);
			CPPUNIT_ASSERT(v4[ind] == a1[ind]);
			CPPUNIT_ASSERT(v5[ind] == INF);
			CPPUNIT_ASSERT(v6[ind] != v6[ind]);
		}

		CPPUNIT_ASSERT(v2 == v3);
		CPPUNIT_ASSERT(v2 == v4);
	}

	void testDataAccess() {
		Vector3<fpType> v1;
		Vector3<fpType> v2;
		const Vector3<fpType> v3(1.0, 2.0, 3.0);

		// Set some values in the matrices and compare them
		for (int i = 0; i < 3; ++i) {
			v1[i] = i;
			v2[i] = i;
			CPPUNIT_ASSERT(v1[i] == i);
			CPPUNIT_ASSERT(v2[i] == i);
			CPPUNIT_ASSERT(v1[i] == v2[i]);
			CPPUNIT_ASSERT(v3[i] == i + 1);
		}

		CPPUNIT_ASSERT(v1[1] == 1.0);
		CPPUNIT_ASSERT(v1[-1] == 0.0);
		CPPUNIT_ASSERT(v1[3] == 2.0);

		v1[1] = 2.0;
		v1[-1] = 1.0;
		v1[3] = 3.0;
		CPPUNIT_ASSERT(v1 == v3);
	}

	void testDataAccessEquality() {
		fpType tmp;
		Vector3<fpType> v1;
		Vector3<fpType> v2;

		for (int k = 0; k < 3; ++k) {
			for (int ind = 0; ind < 3; ++ind) {
				tmp = ind * 3.0;
				v1[ind] = tmp;
				v2[ind] = tmp + (ind != k);
			}
			CPPUNIT_ASSERT(v1 != v2);
		}
	}

	void testEpsEquality() {
		fpType a[3] = { 1.0, -2.0, 3.0 };
		fpType prec = PRECISION * 2.1;
		Vector3<fpType> v1;
		Vector3<fpType> v2;

		// Machine epsilon shouldn't make difference
		v1.set(a);
		for (int ind = 0; ind < 3; ++ind) {
			v2.set(a);
			v2[ind] += PRECISION * 0.9;
			CPPUNIT_ASSERT(v1 == v2);
			v2[ind] -= prec;
			CPPUNIT_ASSERT(v1 != v2);
			CPPUNIT_ASSERT(v1.compare(v2, prec));
		}
	}

	void testIsZero() {
		CPPUNIT_ASSERT(Vector3<fpType>(0.0, 0.0, 0.0).isZero());
		CPPUNIT_ASSERT(Vector3<fpType>(0.1 * PRECISION, 0.0, 0.0).isZero());
		CPPUNIT_ASSERT(Vector3<fpType>(0.0, 0.1 * PRECISION, 0.0).isZero());
		CPPUNIT_ASSERT(Vector3<fpType>(0.0, 0.0, 0.1 * PRECISION).isZero());

		CPPUNIT_ASSERT(!Vector3<fpType>(1.1 * PRECISION, 0.0, 0.0).isZero());
		CPPUNIT_ASSERT(!Vector3<fpType>(0.0, 1.1 * PRECISION, 0.0).isZero());
		CPPUNIT_ASSERT(!Vector3<fpType>(0.0, 0.0, 1.1 * PRECISION).isZero());

		CPPUNIT_ASSERT(!Vector3<fpType>(0.0, 2.0, 1.0).isZero());
		CPPUNIT_ASSERT(!Vector3<fpType>(0.0, -0.1, 0.0).isZero());
		CPPUNIT_ASSERT(!Vector3<fpType>(0.3, 0.0, 9.0).isZero());
	}

	void testAxisAng() {
		fpType a1[][3] = { { 1.0, -2.0, 3.0 }, //
				{ 4.0, -5.0, 3.0 }, //
				{ -8.0, 9.0, -10.0 }, //
				{ 0.0, 18.0, 0.0 }, //
				{ -1.0, 2.0, -3.0 } };
		fpType an1[] = { 0.0, 3.0 * M_PI + PRECISION, -5.0 * M_PI - PRECISION,
				M_PI / 3.0, M_PI * 0.5 };

		for (int i = 0; i < 5; ++i) {
			Vector3<fpType> v1(a1[i], an1[i]);
			Vector3<fpType> v2 = Vector3<fpType>(a1[i]).getUnit() * 2
					* tan(an1[i] / 2);

			CPPUNIT_ASSERT(v1 == v2);
			if (compare<fpType>(normalizeAngle<fpType>(an1[i]), 0.0) != 0) {
				CPPUNIT_ASSERT(
						v1.getUnit()
								== (fpType) sign<fpType>(tan(an1[i] / 2))
										* Vector3<fpType>(a1[i]).getUnit());
			}
			CPPUNIT_ASSERT(
					compare<fpType>(v1.getAngle(),
							fabs(normalizeAngle<fpType>(an1[i]))) == 0);
		}
	}

	void testCopying() {
		Vector3<fpType> v1;
		Vector3<fpType> v2;

		// Set some values in one matrix
		for (int row = 0; row < 3; ++row) {
			v1[row] = row * 2.0;
		}

		CPPUNIT_ASSERT(v2 != v1);
		v2 = v1;
		CPPUNIT_ASSERT(v2 == v1);
		v1[1] = 1.0;
		CPPUNIT_ASSERT(v2 != v1);
	}

	void testCopyConstructor() {
		Vector3<fpType> v1;
		Vector3<fpType> v2;

		// Set some values in one matrix
		for (int ind = 0; ind < 3; ++ind) {
			v1[ind] = ind * 8.0;
			v2[ind] = ind * 7.0;
		}

		Vector3<fpType> v3(v1);

		CPPUNIT_ASSERT(v3 == v1);
		CPPUNIT_ASSERT(v3 != v2);
	}

	void testSetGet() {
		fpType a1[3] = { 1.0, -2.0, 3.0 };
		fpType a2[3];

		Vector3<fpType> v1;
		Vector3<fpType> v2(a1);

		CPPUNIT_ASSERT(v1 != v2);
		v1.set(a1);
		CPPUNIT_ASSERT(v1 == v2);
		v2.get(a2);
		Vector3<fpType> v3(a2);
		CPPUNIT_ASSERT(v3 == v2);

	}

	void testScale() {
		Vector3<fpType> v1;
		Vector3<fpType> v2;

		// Set some values in one matrix
		for (int row = 0; row < 3; ++row) {
			v1[row] = row;
			v2[row] = row * 7.0;
		}

		CPPUNIT_ASSERT(v1 != v2);
		CPPUNIT_ASSERT(v1 * 7.0 == v2);
		CPPUNIT_ASSERT(v1 == (1.0 / 7.0) * v2);
		CPPUNIT_ASSERT(v1 == v2 / 7.0);
		CPPUNIT_ASSERT(v1 / (1.0 / 7.0) == v2);
		v1 *= 7;
		CPPUNIT_ASSERT(v1 == v2);
		v1 /= 2.0;
		v2 *= 0.5;
		CPPUNIT_ASSERT(v1 == v2);
		v1 = 3.0 * v2;
		CPPUNIT_ASSERT(v1 == v2 * 3.0);
	}

	void testCross() {
		fpType a1[3][3] = { { 1.0, 0.0, 0.0 }, //
				{ 0.0, 1.0, 0.0 }, //
				{ 0.0, 0.0, 1.0 } };
		fpType a2[3][3] = { { 3.0, -4.0, -5.0 }, //
				{ -3.0, 4.0, -5.0 }, //
				{ 0.8, 0.6, 0.0 } };

		Vector3<fpType> z;
		Vector3<fpType> e1(a1[0]);
		Vector3<fpType> e1m(a1[0]);
		Vector3<fpType> e2(a1[1]);
		Vector3<fpType> e3(a1[2]);

		Vector3<fpType> r1(a2[0]);
		Vector3<fpType> r2(a2[1]);
		Vector3<fpType> r3(a2[2]);
		r1 /= sqrt(50);
		r2 /= sqrt(50);

		CPPUNIT_ASSERT(z == z.cross(z));
		CPPUNIT_ASSERT(z == e1.cross(z));

		CPPUNIT_ASSERT(e1 == e2.cross(e3));
		CPPUNIT_ASSERT(e2 == e3.cross(e1));
		CPPUNIT_ASSERT(e3 == e1.cross(e2));

		CPPUNIT_ASSERT(r1 == r2.cross(r3));
		CPPUNIT_ASSERT(r2 == r3.cross(r1));
		CPPUNIT_ASSERT(r3 == r1.cross(r2));

		e1m[0] = -2.0;
		CPPUNIT_ASSERT(z == e1m.cross(e1));
	}

	void testDot() {
		fpType a1[3][3] = { { 1.0, 0.0, 0.0 }, //
				{ 0.0, 1.0, 0.0 }, //
				{ 0.0, 0.0, 1.0 } };
		fpType a2[3][3] = { { 2.0, -1.0, -3.0 }, //
				{ 4.0, 1.0, 2.0 } };
		fpType a3[3][3] = { { 3.0, -4.0, -5.0 }, //
				{ -3.0, 4.0, -5.0 } };

		Vector3<fpType> z;
		Vector3<fpType> e1(a1[0]);
		Vector3<fpType> e2(a1[1]);
		Vector3<fpType> e3(a1[2]);

		Vector3<fpType> v1(a2[0]);
		Vector3<fpType> v2(a2[1]);

		Vector3<fpType> v3(a3[0]);
		Vector3<fpType> v4(a3[1]);

		CPPUNIT_ASSERT(0.0 == z.dot(e1));
		CPPUNIT_ASSERT(0.0 == e1.dot(e2));
		CPPUNIT_ASSERT(0.0 == e2.dot(e3));
		CPPUNIT_ASSERT(0.0 == e3.dot(e1));
		CPPUNIT_ASSERT(1.0 == e1.dot(e1));
		CPPUNIT_ASSERT(1.0 == e2.dot(e2));
		CPPUNIT_ASSERT(1.0 == e3.dot(e3));

		CPPUNIT_ASSERT(1.0 == v1.dot(v2));
		CPPUNIT_ASSERT(0.0 == v3.dot(v4));
	}

	void testNorm() {
		fpType a[4][3] = { { -2.0, 1.0, 2.0 }, //
				{ -1.0, 2.0, -2.0 }, //
				{ 0.5, 1.0, 1.0 }, //
				{ -0.8, -0.6, 0.0 } };

		Vector3<fpType> v1(a[0]);
		Vector3<fpType> v2(a[1]);
		Vector3<fpType> v3(a[2]);
		Vector3<fpType> v4(a[3]);
		Vector3<fpType> v4a(a[3]);
		Vector3<fpType> v5;

		CPPUNIT_ASSERT(v1.getNorm() == 3.0);
		CPPUNIT_ASSERT(v1 == 3.0 * v1.getUnit());
		CPPUNIT_ASSERT(v2.getNorm() == 3.0);
		CPPUNIT_ASSERT(v2 == 3.0 * v2.getUnit());
		CPPUNIT_ASSERT(v3.getNorm() == v2.getNorm() * 0.5);
		CPPUNIT_ASSERT(v3 == v3.getUnit() * 1.5);

		CPPUNIT_ASSERT(compare<fpType>(v4.getNorm(), 1.0) == 0);

		v4.normalize();
		CPPUNIT_ASSERT(v4 == v4a);
		CPPUNIT_ASSERT(v4.getUnit() == v4a);

		v1.normalize();
		CPPUNIT_ASSERT(compare<fpType>(v1.getNorm(), 1.0) == 0);
		CPPUNIT_ASSERT(v1[0] == -2.0 * v1[1]);
		CPPUNIT_ASSERT(v1[0] == -v1[2]);

		CPPUNIT_ASSERT(v5.getNorm() == 0);
		v5.normalize();
		CPPUNIT_ASSERT(isNaN<fpType>(v5.getNorm()));
	}

	void testAddSubt() {
		fpType a1[3][3] = { { 1.0, 2.0, 3.0 }, //
				{ -1.0, -2.0, -3.0 }, //
				{ 0.5, 1.0, 1.5 } };

		Vector3<fpType> z;
		Vector3<fpType> v1(a1[0]);
		Vector3<fpType> v1a(a1[0]);
		Vector3<fpType> v2(a1[1]);
		Vector3<fpType> v3(a1[2]);

		CPPUNIT_ASSERT(z == z + z);
		CPPUNIT_ASSERT(v1 + v1 == 2.0 * v1);
		CPPUNIT_ASSERT(v1 + v2 == z);
		CPPUNIT_ASSERT(z - v1 == -v1);
		CPPUNIT_ASSERT(z - v2 == v2 * -1.0);
		v1a += v1a;
		CPPUNIT_ASSERT(v1 * 2.0 == v1a);

		CPPUNIT_ASSERT(z == z - z);
		CPPUNIT_ASSERT(v2 + 2.0 * v3 == z);
		CPPUNIT_ASSERT(v1 - v3 == v3);
		v1a -= v1a;
		CPPUNIT_ASSERT(z == v1a);
	}

	void testInput() {
		string data = "5 6	\n 7.1\n	8.2 123.456";
		istringstream input(data);

		fpType tmp;
		Vector3<fpType> v;
		input >> tmp >> v >> tmp;
		CPPUNIT_ASSERT(v == Vector3<fpType>(6.0, 7.1, 8.2));
		CPPUNIT_ASSERT(compare<fpType>(tmp, 123.456) == 0);
	}
};

CPPUNIT_TEST_SUITE_REGISTRATION (VectorTest);

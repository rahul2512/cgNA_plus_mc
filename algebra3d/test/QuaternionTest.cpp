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
#include<algebra3d/Quaternion.h>
#include<iostream>

#include "Tester.h"

using namespace std;
using namespace algebra3d;

class QuaternionTest: public CppUnit::TestFixture {

CPPUNIT_TEST_SUITE(QuaternionTest);
		CPPUNIT_TEST(testInitEquality);
		CPPUNIT_TEST(testReferences);
		CPPUNIT_TEST(testInitDataAccess);
		CPPUNIT_TEST(testDataAccess);
		CPPUNIT_TEST(testDataAccessEquality);
		CPPUNIT_TEST(testEpsEquality);
		CPPUNIT_TEST(testCopying);
		CPPUNIT_TEST(testCopyConstructor);
		CPPUNIT_TEST(testSetGet);
		CPPUNIT_TEST(testNorm);
		CPPUNIT_TEST(testNeg);
		CPPUNIT_TEST(testMult);
		CPPUNIT_TEST(testInput);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp() {
	}

	void tearDown() {
	}

	void testInitEquality() {
		Quaternion<fpType> q1;
		Quaternion<fpType> q2;
		Quaternion<fpType> q3;

		// Initially quaternions should be identical
		CPPUNIT_ASSERT(q1 == q1);
		CPPUNIT_ASSERT(!(q1 != q1));

		CPPUNIT_ASSERT(q1 == q2);
		CPPUNIT_ASSERT(!(q1 != q2));

		CPPUNIT_ASSERT(q2 == q3);
		CPPUNIT_ASSERT(!(q2 != q3));

	}

	void testReferences() {
		Quaternion<fpType> q1;
		Quaternion<fpType> q2(3.0, 4.0, 5.0, 6.0);

		// Initially quaternions should be identical
		CPPUNIT_ASSERT(q1.x == 0.0);
		CPPUNIT_ASSERT(q1.y == 0.0);
		CPPUNIT_ASSERT(q1.z == 0.0);
		CPPUNIT_ASSERT(q1.w == 1.0);

		CPPUNIT_ASSERT(q2.x == 3.0);
		q2.x *= 2.0;
		CPPUNIT_ASSERT(q2.x == q2[0]);
		q2[0] *= 0.5;
		CPPUNIT_ASSERT(q2.x == q2[0]);

		CPPUNIT_ASSERT(q2.y == 4.0);
		q2.y *= 2.0;
		CPPUNIT_ASSERT(q2.y == q2[1]);
		q2[1] *= 0.5;
		CPPUNIT_ASSERT(q2.y == q2[1]);

		CPPUNIT_ASSERT(q2.z == 5.0);
		q2.z *= 2.0;
		CPPUNIT_ASSERT(q2.z == q2[2]);
		q2[2] *= 0.5;
		CPPUNIT_ASSERT(q2.z == q2[2]);

		CPPUNIT_ASSERT(q2.w == 6.0);
		q2.z *= 2.0;
		CPPUNIT_ASSERT(q2.w == q2[3]);
		q2[2] *= 0.5;
		CPPUNIT_ASSERT(q2.w == q2[3]);
	}

	void testInitDataAccess() {
		fpType a1[4] = { 1.0, -2.0, 3.0, -4.0 };

		Quaternion<fpType> q1;
		Quaternion<fpType> q2(a1);
		Quaternion<fpType> q3(a1[0], a1[1], a1[2], a1[3]);
		Quaternion<fpType> q4(q2);

		for (int ind = 0; ind < 4; ++ind) {
			CPPUNIT_ASSERT(q1[ind] == (ind == 3));
			CPPUNIT_ASSERT(q2[ind] == a1[ind]);
			CPPUNIT_ASSERT(q3[ind] == a1[ind]);
			CPPUNIT_ASSERT(q4[ind] == a1[ind]);
		}

		CPPUNIT_ASSERT(q2 == q3);
		CPPUNIT_ASSERT(q2 == q4);
	}

	void testDataAccess() {
		Quaternion<fpType> q1;
		Quaternion<fpType> q2;
		const Quaternion<fpType> q3(1.0, 2.0, 3.0, 4.0);

		// Set some values in the matrices and compare them
		for (int i = 0; i < 4; ++i) {
			q1[i] = i;
			q2[i] = i;
			CPPUNIT_ASSERT(q1[i] == i);
			CPPUNIT_ASSERT(q2[i] == i);
			CPPUNIT_ASSERT(q1[i] == q2[i]);
			CPPUNIT_ASSERT(q3[i] == i + 1);
		}

		CPPUNIT_ASSERT(q1[1] == 1.0);
		CPPUNIT_ASSERT(q2[2] == 2.0);
		CPPUNIT_ASSERT(q1[-1] == 0.0);
		CPPUNIT_ASSERT(q1[4] == 3.0);
		CPPUNIT_ASSERT(q3[-1] == 1.0);
		CPPUNIT_ASSERT(q3[4] == 4.0);

		q1[1] = 2.0;
		q1[2] = 3.0;
		q1[-1] = 1.0;
		q1[4] = 4.0;
		CPPUNIT_ASSERT(q1 == q3);
	}

	void testDataAccessEquality() {
		fpType tmp;
		Quaternion<fpType> q1;
		Quaternion<fpType> q2;

		for (int k = 0; k < 4; ++k) {
			for (int ind = 0; ind < 4; ++ind) {
				tmp = ind * 3.0;
				q1[ind] = tmp;
				q2[ind] = tmp + (ind != k);
			}
			CPPUNIT_ASSERT(q1 != q2);
		}
	}

	void testEpsEquality() {
		fpType a[4] = { 1.0, -2.0, 3.0, -4.0 };
		fpType prec = PRECISION * 2.1;
		Quaternion<fpType> q1;
		Quaternion<fpType> q2;

		// Machine epsilon shouldn't make difference
		q1.set(a);
		for (int ind = 0; ind < 4; ++ind) {
			q2.set(a);
			q2[ind] += PRECISION * 0.9;
			CPPUNIT_ASSERT(q1 == q2);
			q2[ind] -= prec;
			CPPUNIT_ASSERT(q1 != q2);
			CPPUNIT_ASSERT(q1.compare(q2, prec));
		}
	}

	void testCopying() {
		Quaternion<fpType> q1;
		Quaternion<fpType> q2;

		// Set some values in one matrix
		for (int row = 0; row < 4; ++row) {
			q1[row] = row * 2.0;
		}

		CPPUNIT_ASSERT(q2 != q1);
		q2 = q1;
		CPPUNIT_ASSERT(q2 == q1);
		q1[1] = 1.0;
		CPPUNIT_ASSERT(q2 != q1);
	}

	void testCopyConstructor() {
		Quaternion<fpType> q1;
		Quaternion<fpType> q2;

		// Set some values
		for (int ind = 0; ind < 3; ++ind) {
			q1[ind] = ind * 8.0;
			q2[ind] = ind * 7.0;
		}

		Quaternion<fpType> q3(q1);

		CPPUNIT_ASSERT(q3 == q1);
		CPPUNIT_ASSERT(q3 != q2);
	}

	void testSetGet() {
		fpType a1[4] = { 1.0, -2.0, 3.0, -4.0 };
		fpType a2[4];

		Quaternion<fpType> q1;
		Quaternion<fpType> q2(a1);

		CPPUNIT_ASSERT(q1 != q2);
		q1.set(a1);
		CPPUNIT_ASSERT(q1 == q2);
		q2.get(a2);

		Quaternion<fpType> v3(a2);
		CPPUNIT_ASSERT(v3 == q2);

	}

	void testNorm() {
		Quaternion<fpType> q1;
		Quaternion<fpType> q2(1.0, 1.0, 1.0, 1.0);
		Quaternion<fpType> q3(-0.4, -0.5, -0.6, -0.7);
		Quaternion<fpType> q4;

		CPPUNIT_ASSERT(q1.getNorm() == 1.0);
		CPPUNIT_ASSERT(q1.isUnit());
		CPPUNIT_ASSERT(q2.getNorm() == 2.0);
		CPPUNIT_ASSERT(!q2.isUnit());
		CPPUNIT_ASSERT(
				compare<fpType>(q3.getNorm(),
						sqrt(0.4 * 0.4 + 0.5 * 0.5 + 0.6 * 0.6 + 0.7 * 0.7))
						== 0);
		CPPUNIT_ASSERT(!q3.isUnit());

		q4 = q1;
		q4.normalize();
		CPPUNIT_ASSERT(q4.getNorm() == 1.0);
		CPPUNIT_ASSERT(q4.isUnit());
		CPPUNIT_ASSERT(q4 == q1);

		q4 = q2;
		q4.normalize();
		CPPUNIT_ASSERT(q4.getNorm() == 1.0);
		CPPUNIT_ASSERT(q4.isUnit());

		q4 = q3;
		q4.normalize();
		CPPUNIT_ASSERT(compare<fpType>(q4.getNorm(), 1.0) == 0);
		CPPUNIT_ASSERT(q4.isUnit());
	}

	void testNeg() {
		Quaternion<fpType> id;
		fpType coshalf = cos(M_PI / 6.0);
		fpType comp = 1.0 / sqrt(3) * sin(M_PI / 6.0);
		Quaternion<fpType> q2(comp, comp, comp, coshalf);

		CPPUNIT_ASSERT(id.inv() == id);
		CPPUNIT_ASSERT(q2.inv().getAngle() == q2.getAngle());
		CPPUNIT_ASSERT(id.getNorm() == id.getNorm());
		CPPUNIT_ASSERT(q2.getNorm() == q2.inv().getNorm());
	}

	void testMult() {
		Quaternion<fpType> id;
		fpType coshalf = cos(M_PI / 6.0);
		fpType mult = sin(M_PI / 6.0) / sqrt(14);
		Quaternion<fpType> q2(mult, 2.0 * mult, 3.0 * mult, coshalf);
		Quaternion<fpType> q3(q2);
		coshalf = cos(M_PI / 3.0);
		mult = sin(M_PI / 3.0) / sqrt(14);
		Quaternion<fpType> q4(mult, 2.0 * mult, 3.0 * mult, coshalf);
		mult = 1.0 / sqrt(14);
		Quaternion<fpType> q5(mult, 2.0 * mult, 3.0 * mult, 0.0);

		CPPUNIT_ASSERT(id * q3 == q3);
		CPPUNIT_ASSERT(q4 * id == q4);
		CPPUNIT_ASSERT(q2 * q2.inv() == id);
		CPPUNIT_ASSERT(q3 * q3 == q4);
		CPPUNIT_ASSERT(q2 * q2 * q2 == q5);
		q3 *= q2;
		CPPUNIT_ASSERT(q3 == q4);
		q4 *= q2;
		CPPUNIT_ASSERT(q4 == q5);

		// Got from matlab
		q2 = Quaternion<fpType>(1, 2, 3, 5);
		q3 = Quaternion<fpType>(7, 11, 13, 17);
		q4 = Quaternion<fpType>(45, 97, 113, 17);
		CPPUNIT_ASSERT(q2 * q3 == q4);

	}

	void testInput() {
		string data = "5 6	\n 7.1\n	8.2 9.3 123.456";
		istringstream input(data);

		fpType tmp;
		Quaternion<fpType> q;
		input >> tmp >> q >> tmp;
		CPPUNIT_ASSERT(q == Quaternion<fpType>(6.0, 7.1, 8.2, 9.3));
		CPPUNIT_ASSERT(compare<fpType>(tmp, 123.456) == 0);
	}
};

CPPUNIT_TEST_SUITE_REGISTRATION(QuaternionTest);

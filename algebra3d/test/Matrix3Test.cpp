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
 * A class to test matrix class
 */

#include<cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include<algebra3d/global.h>
#include<algebra3d/Matrix3.h>
#include<algebra3d/Vector3.h>
#include<iostream>

#include "Tester.h"

using namespace std;
using namespace algebra3d;

class MatrixTest: public CppUnit::TestFixture {

CPPUNIT_TEST_SUITE(MatrixTest);
		CPPUNIT_TEST(testInitEquality);
		CPPUNIT_TEST(testInitDataAccess);
		CPPUNIT_TEST(testDataAccess);
		CPPUNIT_TEST(testDataAccessEquality);
		CPPUNIT_TEST(testEpsEquality);
		CPPUNIT_TEST(testCopying);
		CPPUNIT_TEST(testCopyConstructor);
		CPPUNIT_TEST(testScale);
		CPPUNIT_TEST(testSetGet);
		CPPUNIT_TEST(testSetGetCol);
		CPPUNIT_TEST(testTranspose);
		CPPUNIT_TEST(testMult);
		CPPUNIT_TEST(testAddSubt);
		CPPUNIT_TEST(testDet);
		CPPUNIT_TEST(testOrtho);
		CPPUNIT_TEST(testInv);
		CPPUNIT_TEST(testInput);
	CPPUNIT_TEST_SUITE_END();

public:
	void setUp() {
	}

	void tearDown() {
	}

	void testInitEquality() {
		Matrix3<fpType> m1;
		Matrix3<fpType> m2;
		Matrix3<fpType> m3;

		// Initially matrices should be identical
		CPPUNIT_ASSERT(m1 == m1);
		CPPUNIT_ASSERT(!(m1 != m1));

		CPPUNIT_ASSERT(m1 == m2);
		CPPUNIT_ASSERT(!(m1 != m2));

		CPPUNIT_ASSERT(m2 == m3);
		CPPUNIT_ASSERT(!(m2 != m3));

	}

	void testInitDataAccess() {
		fpType a[3][3] = { { 10.0, -2.0, 7.0 }, //
				{ -11.0, 7.0, -3.0 }, //
				{ 82.0, -14.0, 56.0 } };

		Matrix3<fpType> m1;
		Matrix3<fpType> m2(true);
		Matrix3<fpType> m3(a);
		Matrix3<fpType> m4(m3);

		// Initially the matrix should be I or 0
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				CPPUNIT_ASSERT(m1[col][row] == 0.0);
				CPPUNIT_ASSERT(m3[col][row] == a[col][row]);
				CPPUNIT_ASSERT(m4[col][row] == a[col][row]);
				if (row == col) {
					CPPUNIT_ASSERT(m2[col][row] == 1.0);
				} else {
					CPPUNIT_ASSERT(m2[col][row] == 0.0);
				}
			}
		}

		CPPUNIT_ASSERT(m3 == m4);
	}

	void testDataAccess() {
		fpType a[3][3] = { { 9.0, 8.0, 7.0 }, //
				{ 4.0, 5.0, 6.0 }, //
				{ 3.0, 2.0, 1.0 } };

		fpType tmp;
		Matrix3<fpType> m1;
		Matrix3<fpType> m2;
		const Matrix3<fpType> m3(a);

		// Set some values in the matrices and compare them
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				tmp = (col - row) * 3;
				m1[col][row] = tmp;
				m2[col][row] = tmp;
				CPPUNIT_ASSERT(m1[col][row] == tmp);
				CPPUNIT_ASSERT(m2[col][row] == tmp);
				CPPUNIT_ASSERT(m1[col][row] == m2[col][row]);
				CPPUNIT_ASSERT(m3[col][row] == a[col][row]);
			}
		}

		CPPUNIT_ASSERT(m1[1][1] == 0.0);
		CPPUNIT_ASSERT(m2[2][2] == 0.0);
		CPPUNIT_ASSERT(m1[-1][1] == -3.0);
		CPPUNIT_ASSERT(m1[3][0] == 6.0);

		CPPUNIT_ASSERT(m1 == m2);
	}

	void testDataAccessEquality() {
		fpType tmp;
		Matrix3<fpType> m1;
		Matrix3<fpType> m2;

		// Set some values in the matrices and compare them

		for (int k = 0; k < 9; ++k) {
			for (int row = 0; row < 3; ++row) {
				for (int col = 0; col < 3; ++col) {
					tmp = (col + row) * 3.0;
					m1[col][row] = tmp;
					m2[col][row] = tmp + ((col * 3 + row) != k);
				}
			}
			CPPUNIT_ASSERT(m1 != m2);
		}
	}

	void testEpsEquality() {
		fpType a[3][3] = { { 1.0, -2.0, 3.0 }, //
				{ -4.0, 5.0, -6.0 }, //
				{ 7.0, -8.0, 9.0 } };

		fpType prec = PRECISION * 2.1;
		Matrix3<fpType> m1;
		Matrix3<fpType> m2(true);

		// Default precision shouldn't make difference
		m1.set(a);
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				m2.set(a);
				m2[col][row] += PRECISION * 0.9;
				CPPUNIT_ASSERT(m1 == m2);
				m2[col][row] -= prec;
				CPPUNIT_ASSERT(m1 != m2);
				CPPUNIT_ASSERT(m1.compare(m2, prec));
			}
		}
	}

	void testCopying() {
		Matrix3<fpType> m1;
		Matrix3<fpType> m2;

		// Set some values in one matrix
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				m1[col][row] = (col - row) * 2.0;
			}
		}

		CPPUNIT_ASSERT(m2 != m1);
		m2 = m1;
		CPPUNIT_ASSERT(m2 == m1);
		m1[1][1] = 1.0;
		CPPUNIT_ASSERT(m2 != m1);
	}

	void testCopyConstructor() {
		Matrix3<fpType> m1;
		Matrix3<fpType> m2;

		// Set some values in one matrix
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				m1[col][row] = (col - row) * 8.0;
				m2[col][row] = (col - row) * 7.0;
			}
		}

		Matrix3<fpType> m3(m1);

		CPPUNIT_ASSERT(m3 == m1);
		CPPUNIT_ASSERT(m3 != m2);
	}

	void testSetGet() {
		fpType a1[3][3] = { { 1.0, -2.0, 3.0 }, //
				{ -4.0, 5.0, -6.0 }, //
				{ 7.0, -8.0, 9.0 } };
		fpType a2[3][3];
		fpType a3[11] = { 22.2, //
				1.0, -2.0, 3.0, //
				-4.0, 5.0, -6.0, //
				7.0, -8.0, 9.0, //
				33.3};

		Vector3<fpType> v1(a1[0]);
		Vector3<fpType> v2(a1[1]);
		Vector3<fpType> v3(a1[2]);

		Matrix3<fpType> m1;
		Matrix3<fpType> m2(a1);
		Matrix3<fpType> m3(a3 + 1);
		Matrix3<fpType> m4;
		Matrix3<fpType> m5;

		CPPUNIT_ASSERT(m1 != m2);
		CPPUNIT_ASSERT(m2 == m3);
		m1.set(a1);
		CPPUNIT_ASSERT(m1 == m2);

		CPPUNIT_ASSERT(m4 != m2);
		m2.get(a2);
		m4.set(a2);
		CPPUNIT_ASSERT(m4 == m2);

		CPPUNIT_ASSERT(m5 != m2);
		m5.set(v1, v2, v3);
		CPPUNIT_ASSERT(m5 == m2);
	}

	void testSetGetCol() {
		fpType a1[3][3] = { { 1.0, -2.0, 3.0 }, //
				{ -4.0, 5.0, -6.0 }, //
				{ 7.0, -8.0, 9.0 } };
		fpType a2[3][3] = { { -111.0, 222.0, -333.0 }, //
				{ 444.0, -555.0, 666.0 }, //
				{ -777.0, 888.0, -999.0 } };
		fpType a3[3] = {0.0, 0.0, 0.0};

		Matrix3<fpType> m1(a1);
		Matrix3<fpType> m2(a2);
		Matrix3<fpType> m3;

		CPPUNIT_ASSERT(m1 != m2);

		m1.setColumn(0, m2[0]);
		m1.getColumn(0, a3);
		CPPUNIT_ASSERT(m1 != m2);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[0]) == m2[0]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[1]) != m2[1]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[2]) != m2[2]);

		m1.setColumn(1, m2[1]);
		m1.getColumn(1, a3);
		CPPUNIT_ASSERT(m1 != m2);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[0]) == m2[0]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[1]) == m2[1]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[2]) != m2[2]);

		m1.setColumn(2, m2[2]);
		m1.getColumn(2, a3);
		CPPUNIT_ASSERT(m1 == m2);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[0]) == m2[0]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[1]) == m2[1]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[2]) == m2[2]);

		m1.setColumn(1, m2[2]);
		m1.getColumn(1, a3);
		CPPUNIT_ASSERT(m1 != m2);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[0]) == m2[0]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[1]) == m2[2]);
		CPPUNIT_ASSERT(Vector3<fpType>(m1[2]) == m2[2]);
	}

	void testScale() {
		Matrix3<fpType> m1;
		Matrix3<fpType> m2;

		// Set some values in one matrix
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				m1[col][row] = (col - row);
				m2[col][row] = (col - row) * 7.0;
			}
		}

		CPPUNIT_ASSERT(m1 != m2);
		CPPUNIT_ASSERT(m1 * 7.0 == m2);
		CPPUNIT_ASSERT(m1 == (1.0 / 7.0) * m2);
		CPPUNIT_ASSERT(m1 == m2 / 7.0);
		CPPUNIT_ASSERT(m1 / (1.0 / 7.0) == m2);
		m1 *= 7;
		CPPUNIT_ASSERT(m1 == m2);
		m1 /= 2.0;
		m2 *= 0.5;
		CPPUNIT_ASSERT(m1 == m2);
		m1 = 3.0 * m2;
		CPPUNIT_ASSERT(m1 == m2 * 3.0);
	}

	void testTranspose() {
		fpType a1[3][3] = { { 1.0, -2.0, 3.0 }, //
				{ -2.0, 5.0, -6.0 }, //
				{ 3.0, -6.0, 9.0 } };
		fpType a2[3][3] = { { 1.5, -4.5, 7.5 }, //
				{ -2.5, 5.5, -8.5 }, //
				{ 3.5, -6.5, 9.5 } };
		fpType a3[3][3] = { { 1.5, -2.5, 3.5 }, //
				{ -4.5, 5.5, -6.5 }, //
				{ 7.5, -8.5, 9.5 } };

		Matrix3<fpType> z;
		Matrix3<fpType> id(true);
		Matrix3<fpType> m1(a1);
		Matrix3<fpType> m2(a2);
		Matrix3<fpType> m3(a3);

		CPPUNIT_ASSERT(z == z.getTranspose());
		CPPUNIT_ASSERT(id == id.getTranspose());
		CPPUNIT_ASSERT(m1 == m1.getTranspose());
		CPPUNIT_ASSERT(m2 == m3.getTranspose());
		CPPUNIT_ASSERT(m3 == m2.getTranspose());

		m3.transpose();
		CPPUNIT_ASSERT(m2 == m3);

		m2.transpose();
		CPPUNIT_ASSERT(m2 == m3.getTranspose());
	}

	void testMult() {
		fpType a1[3][3] = { { 1.0, -2.0, 3.0 }, //
				{ -4.0, 5.0, -6.0 }, //
				{ 7.0, -8.0, 9.0 } };
		fpType a2[3][3] = { { 1.5, -4.5, 7.5 }, //
				{ -2.5, 5.5, -8.5 }, //
				{ 3.5, -6.5, 9.5 } };
		fpType a3[3][3] = { { 72.0, -85.5, 99.0 }, //
				{ -84.0, 100.5, -117.0 }, //
				{ 96.0, -115.5, 135.0 } };
		fpType a4[3][3] = { { 17.0, -35.0, 53.0, }, //
				{ -39.5, 84.5, -129.5 }, //
				{ 62.0, -134.0, 206.000 } };

		Matrix3<fpType> z;
		Matrix3<fpType> id(true);
		Matrix3<fpType> m1(a1);
		Matrix3<fpType> m1a(a1);
		Matrix3<fpType> m2(a2);
		Matrix3<fpType> m3(a3);
		Matrix3<fpType> m4(a4);

		CPPUNIT_ASSERT(z == z * z);
		CPPUNIT_ASSERT(z == z * id);
		CPPUNIT_ASSERT(z == id * z);
		CPPUNIT_ASSERT(id == id * id);
		CPPUNIT_ASSERT(z == z * m1);
		CPPUNIT_ASSERT(z == m1 * z);
		CPPUNIT_ASSERT(m1 == id * m1);
		CPPUNIT_ASSERT(m1 == m1 * id);
		CPPUNIT_ASSERT(m3 == m1 * m2);
		CPPUNIT_ASSERT(m2 * m1 == m4);
		m1a *= m2;
		CPPUNIT_ASSERT(m1a == m3);
	}

	void testAddSubt() {
		fpType a1[3][3] = { { 1.0, 2.0, 3.0 }, //
				{ 2.0, 1.0, 2.0 }, //
				{ 3.0, 2.0, 1.0 } };
		fpType a2[3][3] = { { 3.0, 0.0, 0.0 }, //
				{ 2.0, 1.0, -1.0 }, //
				{ 2.0, 0.0, 2.0 } };

		Matrix3<fpType> z;
		Matrix3<fpType> id(true);
		Matrix3<fpType> m1(a1);
		Matrix3<fpType> m1a(a1);
		Matrix3<fpType> m2(a2);

		CPPUNIT_ASSERT(z == z + z);
		CPPUNIT_ASSERT(id + id == 2.0 * id);
		CPPUNIT_ASSERT(z - m1 == -m1);
		CPPUNIT_ASSERT(z - m2 == m2 * -1.0);
		CPPUNIT_ASSERT(m1 + m1.getTranspose() == 2.0 * m1);
		m1a += m1a;
		CPPUNIT_ASSERT(m1 * 2.0 == m1a);

		CPPUNIT_ASSERT(z == z - z);
		CPPUNIT_ASSERT(id - id == z);
		CPPUNIT_ASSERT(m1 - m1 == z);
		CPPUNIT_ASSERT(m1 - m1.getTranspose() == z);
		CPPUNIT_ASSERT(-m1 + m1.getTranspose() == z);
		m1a -= m1a;
		CPPUNIT_ASSERT(z == m1a);
	}

	void testDet() {
		fpType a1[3][3] = { { 1.0, -2.0, 3.0 }, //
				{ -4.0, 5.0, -6.0 }, //
				{ 7.0, -8.0, 9.0 } };
		fpType a2[3][3] = { { 3.0, 0.0, 0.0 }, //
				{ 2.0, 1.0, -1.0 }, //
				{ 2.0, 0.0, 2.0 } };

		Matrix3<fpType> z;
		Matrix3<fpType> id(true);
		Matrix3<fpType> m1(a1);
		Matrix3<fpType> m2(a2);

		CPPUNIT_ASSERT(z.det() == 0.0);
		CPPUNIT_ASSERT(id.det() == 1.0);
		CPPUNIT_ASSERT(m1.det() == 0.0);
		CPPUNIT_ASSERT(m2.det() == 6.0);
	}

	void testOrtho() {
		fpType a1[3][3] = { { 3.0, 0.0, 0.0 }, //
				{ 2.0, 1.0, -1.0 }, //
				{ 2.0, 0.0, 2.0 } };
		fpType a2[3][3] = { { 0.0, 1.0, 0.0 }, //
				{ -1.0, 0.0, 0.0 }, //
				{ 0.0, 0.0, 1.0 } };
		fpType cosp = cos(M_PI * 0.25);
		fpType a3[3][3] = { { cosp, cosp, 0.0 }, //
				{ -cosp, cosp, 0.0 }, //
				{ 0.0, 0.0, 1.0 } };

		Matrix3<fpType> z;
		Matrix3<fpType> id(true);
		Matrix3<fpType> m1(a1);
		Matrix3<fpType> m2(a2);
		Matrix3<fpType> m3(a3);

		CPPUNIT_ASSERT(!z.isOrthogonal());
		CPPUNIT_ASSERT(id.isOrthogonal());
		CPPUNIT_ASSERT(!m1.isOrthogonal());
		CPPUNIT_ASSERT(m1 * m1.getTranspose() != id);
		CPPUNIT_ASSERT(m2.isOrthogonal());
		CPPUNIT_ASSERT(m2 * m2.getTranspose() == id);
		CPPUNIT_ASSERT(m3.isOrthogonal());
		CPPUNIT_ASSERT(m3 * m3.getTranspose() == id);
	}

	void testInv() {
		fpType a1[3][3] = { { 3.0, 0.0, 0.0 }, //
				{ 2.0, 1.0, -1.0 }, //
				{ 2.0, 0.0, 2.0 } };
		fpType a1a[3][3] = { { 1.0 / 3.0, 0.0, 0.0 }, //
				{ -1.0, 1.0, 0.5 }, //
				{ -1.0 / 3.0, 0.0, 0.5 } };
		fpType a2[3][3] = { { 0.0, 0.0, -1.0 }, //
				{ 0.0, -1.0, 0.0 }, //
				{ -1.0, 0.0, 0.0 } };

		Matrix3<fpType> z;
		Matrix3<fpType> id(true);
		Matrix3<fpType> m(a1);
		Matrix3<fpType> mInv(a1a);
		Matrix3<fpType> m2(a2);

		CPPUNIT_ASSERT(z.getInverse() == z);
		CPPUNIT_ASSERT(id.getInverse() == id);
		CPPUNIT_ASSERT(m.getInverse() == mInv);
		CPPUNIT_ASSERT(mInv.getInverse() == m);
		id.invert();
		CPPUNIT_ASSERT(id.getInverse() == id);
		m.invert();
		mInv.invert();
		CPPUNIT_ASSERT(m.getInverse() == mInv);
		CPPUNIT_ASSERT(mInv.getInverse() == m);

		CPPUNIT_ASSERT(m2.getInverse() == m2.getTranspose());
	}

	void testInput() {
		fpType a[3][3] = { { 1, 2.6, 3.5 }, //
				{ 4, 5.4, 6 }, //
				{ 7.3, 8.2, 9.1 } };

		string data = "5 \t\n 1\n\t4.0 7.3\n\n2.600 5.4 8.2" //
						"\n\t3.5\n\t 6\n\t\t9.1 123.456";
		istringstream input(data);

		fpType tmp;
		Matrix3<fpType> m;
		input >> tmp >> m >> tmp;
		CPPUNIT_ASSERT(m == Matrix3<fpType>(a));
		CPPUNIT_ASSERT(compare<fpType>(tmp, 123.456) == 0);
	}
};

CPPUNIT_TEST_SUITE_REGISTRATION(MatrixTest);

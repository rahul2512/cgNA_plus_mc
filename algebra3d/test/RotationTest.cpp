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
 * A class to test operations involving rotations using matrices, vectors and
 * quaternions
 */
#include<cstdlib>
#include<sys/time.h>

#include<iostream>
#include<vector>

#include<cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include<algebra3d/global.h>
#include<algebra3d/Matrix3.h>
#include<algebra3d/Vector3.h>
#include<algebra3d/Quaternion.h>

#include "Tester.h"

using namespace std;
using namespace algebra3d;

class RotationTest: public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE (RotationTest);
	CPPUNIT_TEST (testCreation);
	CPPUNIT_TEST (testOuter);
	CPPUNIT_TEST (testRotVec);
	CPPUNIT_TEST (testMultMatVec);
	CPPUNIT_TEST (testCrossMatrix);
	CPPUNIT_TEST (testConversion);
	CPPUNIT_TEST (testAngleMult);
	CPPUNIT_TEST (testSqrtRot);
	CPPUNIT_TEST (testRotTowards);
	CPPUNIT_TEST (testMatCols);CPPUNIT_TEST_SUITE_END();

public:
	void setUp() {
		// Make the test reproducible
		srand(0);
	}

	void tearDown() {
		// Make sure the output comes out as soon as it is available
		cout.flush();
	}

	void checkAngleAxis(Vector3<fpType> axis, fpType angle) {
		Vector3<fpType> unitAxis = axis.getUnit();
		fpType normalizedAngle = normalizeAngle(angle);
		Vector3<fpType> axisAngle = axis * normalizedAngle;

		Vector3<fpType> v = Vector3<fpType>(axis, angle);
		Quaternion<fpType> qvSqrt = v.getSqrtRotation();
		Quaternion<fpType> q = Quaternion<fpType>(axis, angle);
		Quaternion<fpType> qNeg(-q.x, -q.y, -q.z, -q.w);
		Quaternion<fpType> q2 = Quaternion<fpType>(v);
		Quaternion<fpType> qSqrt = q.getSqrtRotation();
		Quaternion<fpType> qNegSqrt = qNeg.getSqrtRotation();
		Matrix3<fpType> m = Matrix3<fpType>(axis, angle);
		Matrix3<fpType> m2 = Matrix3<fpType>(v);
		Matrix3<fpType> m3 = Matrix3<fpType>(m[0], m[1], m[2]);
		Matrix3<fpType> m4;
		m4.set(m[0], m[1], m[2]);

		CPPUNIT_ASSERT(m.isOrthogonal());
		CPPUNIT_ASSERT(compare<fpType>(q.getNorm(), 1.0) == 0);

		CPPUNIT_ASSERT(q.compareRotations(Quaternion<fpType>(m)));

		CPPUNIT_ASSERT(m3 == m);
		CPPUNIT_ASSERT(m3 == m4);

//		cerr << q << "   ->  " << qSqrt * qSqrt << "  " << qSqrt << endl;
//		cerr << qNeg << "  ->  " << qNegSqrt * qNegSqrt << "  " << qNegSqrt
//				<< endl;
//				cerr << "_________" << (-0.0 > 0.0) << endl;

		CPPUNIT_ASSERT(q == qSqrt * qSqrt);
		CPPUNIT_ASSERT(q.w >= 0);

		// Rotation by pi or -pi
		if (compare<fpType>(normalizedAngle, M_PI) == 0
				|| compare<fpType>(normalizedAngle, -M_PI) == 0) {
			CPPUNIT_ASSERT(qNeg == qNegSqrt * qNegSqrt);

			CPPUNIT_ASSERT(Matrix3<fpType>(q) == m);
			CPPUNIT_ASSERT(
					(m.getAxis() * m.getAngle() == //
							unitAxis * normalizedAngle)
							|| (m.getAxis() * m.getAngle() == //
									-unitAxis * normalizedAngle));
		}
		// Not pi nor -pi
		else {
			CPPUNIT_ASSERT(q == qNegSqrt * qNegSqrt);

			CPPUNIT_ASSERT(Vector3<fpType>(q).compare(v, PRECISION * 5));
			CPPUNIT_ASSERT(Vector3<fpType>(qNeg).compare(v, PRECISION * 5));
			CPPUNIT_ASSERT(Vector3<fpType>(m).compare(v, PRECISION * 10));

			CPPUNIT_ASSERT(Quaternion<fpType>(v).compareRotations(q));
			CPPUNIT_ASSERT(Quaternion<fpType>(v).compareRotations(qNeg));
			CPPUNIT_ASSERT(Matrix3<fpType>(v) == m);

			CPPUNIT_ASSERT(q.compareRotations(q2));
			CPPUNIT_ASSERT(m == m2);
			//cerr << v << Vector3<fpType>(qvSqrt * qvSqrt) << endl;
			CPPUNIT_ASSERT(v.compare(qvSqrt * qvSqrt, PRECISION * 2));

			// Not pi, -pi nor 0
			if (compare<fpType>(normalizedAngle, 0.0) != 0) {
				CPPUNIT_ASSERT(v.getUnit() * v.getAngle() == //
						unitAxis * normalizedAngle);
				v = Vector3<fpType>(q);
				CPPUNIT_ASSERT(v.getUnit() == q.getAxis());
				CPPUNIT_ASSERT(v.getUnit() == qNeg.getAxis());
				v = Vector3<fpType>(m);
				CPPUNIT_ASSERT(v.getUnit() == m.getAxis());

				CPPUNIT_ASSERT(m.getAxis() == q.getAxis());
				CPPUNIT_ASSERT(m.getAxis() == qNeg.getAxis());

				CPPUNIT_ASSERT(
						q.getAxis() * q.getAngle()
								== unitAxis * normalizedAngle);
				CPPUNIT_ASSERT(
						qNeg.getAxis() * qNeg.getAngle()
								== unitAxis * normalizedAngle);
			}
		}
	}

	void testCreation() {
		int n = 11;
		fpType a[][3] = { { 3.9, -0.5, 0.3 }, //
				{ 2.0, 1.0, -1.0 }, //
				{ -2.1, 0.8, 2.2 }, //
				//
				{ 3.9, 0.8, -0.9 }, //
				{ -0.1, 4.8, 0.7 }, //
				{ -2.1, 0.8, 2.2 }, //
				{ 2.9, 0.8, 0.15 }, //
				{ -0.5, 3.0, -0.5 }, //
				{ 0.1, 0.2, 0.3 }, //
				{ 0.0, 0.0, 2.2 }, //
				{ 0.123, 0.3, 8.2 } };
		fpType an[] = { 0.0, M_PI, 3.4 * M_PI, //
				0.32 * M_PI, -0.2 * M_PI, 0.15 * M_PI, //
				0.80 * M_PI, -0.85 * M_PI, 0.75 * M_PI, //
				-0.9 * M_PI, -0.8 * M_PI };

		fpType a3[5][4] = { { 1.1, -1.2, 1.3, -1.4 }, //
				{ 1.0, -2.0, 3.0, 9.9 }, //
				{ -4.0, 5.0, -6.0, 8.8 }, //
				{ 7.0, -8.0, 9.0, 7.7 }, //
				{ 6.6, -5.5, -4.4, 3.3 } };

		Vector3<fpType> b1 = Vector3<fpType>(a3[1]);
		Vector3<fpType> b2 = Vector3<fpType>(a3[2]);
		Vector3<fpType> b3 = Vector3<fpType>(a3[3]);
		Matrix3<fpType> m(b1, b2, b3);
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				CPPUNIT_ASSERT(m[i][j] == a3[i + 1][j]);
			}
		}

		Vector3<fpType> axis;
		Vector3<fpType> v;
		Quaternion<fpType> q;
		Quaternion<fpType> q2;
		Matrix3<fpType> m2;
		Matrix3<fpType> m3;
		Matrix3<fpType> m4;

		for (int i = 0; i < n; ++i) {
			checkAngleAxis(Vector3<fpType>(a[i]), an[i]);
		}

		if (DO_HEAVY_TESTS) {
			char point = '*';
			double rangeMax = 2.0;
			double rangeMin = -2.0;
			double step = 0.1;
			int steps = 1 + static_cast<int>((rangeMax - rangeMin) / step);

			cout << "\nStarting a very long testing loop!\n" << steps
					<< " lines of " << steps << " '" << point
					<< "' will be printed:" << endl;
			for (fpType x = rangeMin; x < rangeMax; x += step) {
				for (fpType y = rangeMin; y < rangeMax; y += step) {
					cout << point;
					cout.flush();
					for (fpType z = rangeMin; z < rangeMax; z += step) {
						if ((compare<fpType>(x, 0.0) != 0) //
						|| (compare<fpType>(y, 0.0) != 0) //
								|| (compare<fpType>(z, 0.0) != 0)) {
							for (fpType a = rangeMin; a < rangeMax; a += step) {
								checkAngleAxis(Vector3<fpType>(x, y, z),
										a * M_PI);
							}
						}
					}
				}
				cout << endl;
			}
			cout << "Heavy testing done!" << endl;
		}
	}

	void testOuter() {
		fpType a[][3] = { { 1.0, 2.0, 3.0 }, //
				{ 5.0, 7.0, 11.0 } };
		fpType a2[][3][3] = { { { 5.0, 10.0, 15.0 }, //
				{ 7.0, 14.0, 21.0 }, //
				{ 11.0, 22.0, 33.0 } //
		} //
		};
		Vector3<fpType> v1(a[0]);
		Vector3<fpType> v2(a[1]);
		Matrix3<fpType> m(a2[0]);

		CPPUNIT_ASSERT(v1.outer(v2) == m);
	}

	void testMultMatVec() {
		fpType a1[3][3] = { { 3.9, -0.5, 0.3 }, //
				{ 2.0, 1.0, -1.0 }, //
				{ -2.1, 0.8, 2.2 } };
		fpType a1v[3][3] = { { 1.1, 2.2, 3.3 }, //
				{ 1.76, 4.29, 5.39 }, //
				{ 4.18, 1.1, 6.71 } };

		fpType a2[3][3] = { { -0.1, 0.5, 0.2 }, //
				{ -1.0, 1.0, -1.0 }, //
				{ -2.0, -1.2, -3.2 } };
		fpType a2v[3][3] = { { 4.0, 8.0, 4.0 }, //
				{ -16.4, 5.2, -20.0 }, //
				{ 4.4, 0.0, -30.4 } };

		fpType cosp = cos(M_PI * 0.25);
		fpType a3[3][3] = { { cosp, cosp, 0.0 }, //
				{ -cosp, cosp, 0.0 }, //
				{ 0.0, 0.0, 1.0 } };
		fpType sqrt2 = sqrt(2);
		fpType a3v[3][3] = { { 1.0, 1.0, 1.0 }, //
				{ 0.0, sqrt2, 1.0 }, //
				{ sqrt2, 0.0, 1.0 } };

		Vector3<fpType> zv;
		Matrix3<fpType> zm;
		Matrix3<fpType> id(true);

		Matrix3<fpType> m1(a1);
		Vector3<fpType> v1(a1v[0]);
		Vector3<fpType> m1v1(a1v[1]);
		Vector3<fpType> v1m1(a1v[2]);

		Matrix3<fpType> m2(a2);
		Vector3<fpType> v2(a2v[0]);
		Vector3<fpType> m2v2(a2v[1]);
		Vector3<fpType> v2m2(a2v[2]);

		Matrix3<fpType> m3(a3);
		Vector3<fpType> v3(a3v[0]);
		Vector3<fpType> m3v3(a3v[1]);
		Vector3<fpType> v3m3(a3v[2]);

		CPPUNIT_ASSERT(zm * v1 == zv);
		CPPUNIT_ASSERT(v1 * zm == zv);
		CPPUNIT_ASSERT(id * v2 == v2);
		CPPUNIT_ASSERT(v2 * id == v2);

		CPPUNIT_ASSERT(v1 * m1 == v1m1);
		CPPUNIT_ASSERT(m1 * v1 == m1v1);

		CPPUNIT_ASSERT(m2 * v2 == m2v2);
		CPPUNIT_ASSERT(v2 * m2 == v2m2);

		CPPUNIT_ASSERT(m3 * v3 == m3v3);
		CPPUNIT_ASSERT(v3 * m3 == v3m3);
	}

	void testCrossMatrix() {
		fpType a[3][3] = { { 0.0, 5.0, -3.0 }, //
				{ -5.0, 0.0, 2.0 }, //
				{ 3.0, -2.0, 0.0 } };

		Vector3<fpType> v1(2.0, 3.0, 5.0);
		Vector3<fpType> v2(7.0, 11.0, 13.0);
		Matrix3<fpType> m1(a);

		CPPUNIT_ASSERT(m1 == v1.crossMatrix());
		CPPUNIT_ASSERT(v1.crossMatrix() * v2 == v1.cross(v2));
	}

	void testRotVec() {
		Vector3<fpType> id;
		Vector3<fpType> axis(1.0, 1.0, 1.0);
		fpType comp = 1.0 / sqrt(3.0);
		Vector3<fpType> axisn(comp, comp, comp);
		Vector3<fpType> e1(1.0, 0.0, 0.0);
		Vector3<fpType> e2(0.0, 1.0, 0.0);
		Vector3<fpType> e3(0.0, 0.0, 1.0);
		fpType angle = M_PI / 6.0;

		Vector3<fpType> rotV1(axis, angle);
		Vector3<fpType> rotV2(axis, 2.0 * angle);
		Vector3<fpType> rotV3(axis, 3.0 * angle);
		Vector3<fpType> rotV4(axis, 4.0 * angle);

		Quaternion<fpType> rotQ1(axis, angle);
		Quaternion<fpType> rotQ2(axis, 2.0 * angle);
		Quaternion<fpType> rotQ3(axis, 3.0 * angle);
		Quaternion<fpType> rotQ4(axis, 4.0 * angle);

		Matrix3<fpType> rotM1(axis, angle);
		Matrix3<fpType> rotM2(axis, 2.0 * angle);
		Matrix3<fpType> rotM3(axis, 3.0 * angle);
		Matrix3<fpType> rotM4(axis, 4.0 * angle);

		CPPUNIT_ASSERT(id * id == id);
		CPPUNIT_ASSERT(rotV3 * id == id);
		CPPUNIT_ASSERT(id * rotV3 == rotV3);

		CPPUNIT_ASSERT(rotV4 * e1 == e2);
		CPPUNIT_ASSERT(rotV4 * e2 == e3);
		CPPUNIT_ASSERT(rotV4 * e3 == e1);
		CPPUNIT_ASSERT(rotV4 * (rotV4 * e1) == e3);
		CPPUNIT_ASSERT(rotV4 * (rotV4 * e2) == e1);
		CPPUNIT_ASSERT(rotV4 * (rotV4 * e3) == e2);
		CPPUNIT_ASSERT(rotV4 * (rotV4 * (rotV4 * e1)) == e1);
		CPPUNIT_ASSERT(rotV4 * (rotV4 * (rotV4 * e2)) == e2);
		CPPUNIT_ASSERT(rotV4 * (rotV4 * (rotV4 * e3)) == e3);

		CPPUNIT_ASSERT(rotQ4 * e1 == e2);
		CPPUNIT_ASSERT(rotQ4 * e2 == e3);
		CPPUNIT_ASSERT(rotQ4 * e3 == e1);
		CPPUNIT_ASSERT(rotQ4 * rotQ4 * e1 == e3);
		CPPUNIT_ASSERT(rotQ4 * rotQ4 * e2 == e1);
		CPPUNIT_ASSERT(rotQ4 * rotQ4 * e3 == e2);
		CPPUNIT_ASSERT(rotQ4 * rotQ4 * rotQ4 * e1 == e1);
		CPPUNIT_ASSERT(rotQ4 * rotQ4 * rotQ4 * e2 == e2);
		CPPUNIT_ASSERT(rotQ4 * rotQ4 * rotQ4 * e3 == e3);

		CPPUNIT_ASSERT(rotM4 * e1 == e2);
		CPPUNIT_ASSERT(rotM4 * e2 == e3);
		CPPUNIT_ASSERT(rotM4 * e3 == e1);
		CPPUNIT_ASSERT(rotM4 * rotM4 * e1 == e3);
		CPPUNIT_ASSERT(rotM4 * rotM4 * e2 == e1);
		CPPUNIT_ASSERT(rotM4 * rotM4 * e3 == e2);
		CPPUNIT_ASSERT(rotM4 * rotM4 * rotM4 * e1 == e1);
		CPPUNIT_ASSERT(rotM4 * rotM4 * rotM4 * e2 == e2);
		CPPUNIT_ASSERT(rotM4 * rotM4 * rotM4 * e3 == e3);

		// Vectors co-linear with rotation axis do not change
		CPPUNIT_ASSERT(rotV3 * rotV1 == rotV1);
		CPPUNIT_ASSERT(rotV1 * rotV3 == rotV3);
		CPPUNIT_ASSERT(rotV3 * rotV2 == rotV2);
		CPPUNIT_ASSERT(rotV1 * rotV4 == rotV4);

		CPPUNIT_ASSERT(rotQ3 * rotV1 == rotV1);
		CPPUNIT_ASSERT(rotQ1 * rotV3 == rotV3);
		CPPUNIT_ASSERT(rotQ3 * rotV2 == rotV2);
		CPPUNIT_ASSERT(rotQ1 * rotV4 == rotV4);

		CPPUNIT_ASSERT(rotM3 * rotV1 == rotV1);
		CPPUNIT_ASSERT(rotM1 * rotV3 == rotV3);
		CPPUNIT_ASSERT(rotM3 * rotV2 == rotV2);
		CPPUNIT_ASSERT(rotM1 * rotV4 == rotV4);
	}

	void testConversion() {
		int numTests = 1000;
		Vector3<fpType> axis(1.987654321987654321, 2.48163264128,
				7.654321987654321);
		fpType angle = 0.123456789123456789123456789;

		for (int i = 0; i < numTests; ++i) {
			Vector3<fpType> rotV(axis, angle);

			Quaternion<fpType> rotQ(axis, angle);

			Matrix3<fpType> rotM(axis, angle);

			CPPUNIT_ASSERT(rotV == Vector3<fpType>(rotM));
			CPPUNIT_ASSERT(rotV == Vector3<fpType>(rotQ));

			CPPUNIT_ASSERT(rotM == Matrix3<fpType>(rotV));
			CPPUNIT_ASSERT(rotM == Matrix3<fpType>(rotQ));

			CPPUNIT_ASSERT(rotQ == Quaternion<fpType>(rotM));
			CPPUNIT_ASSERT(rotQ == Quaternion<fpType>(rotV));

			CPPUNIT_ASSERT(rotQ == Matrix3<fpType>(Quaternion<fpType>(rotM)));
			CPPUNIT_ASSERT(rotQ == Vector3<fpType>(Quaternion<fpType>(rotV)));

			angle = -0.5 * M_PI + (M_PI * i) / numTests;
			axis.x = (rand() / (double) RAND_MAX);
			axis.y = (rand() / (double) RAND_MAX);
			axis.z = (rand() / (double) RAND_MAX);
			axis.normalize();
		}
	}

	void testAngleMult() {
		fpType angle = 1.38692759366;
		Vector3<fpType> axis = Vector3<fpType>(2.3, 4.5, 9.3);
		Quaternion<fpType> id;
		Quaternion<fpType> id2(0.0, 0.0, 0.0, -1.0);
		Quaternion<fpType> q1(axis, angle);
		Quaternion<fpType> q2(axis, angle * 0.25);
		Quaternion<fpType> q3(axis, angle * 0.5);
		Quaternion<fpType> q4 = q1;
		q4.scaleAngle(1 / 3.0);

		Quaternion<fpType> q = q1.getScaledAngleRotation(0.25);

		CPPUNIT_ASSERT(id.getScaledAngleRotation(0.4) == id);
		CPPUNIT_ASSERT(id2.getScaledAngleRotation(8.6) == id2);
		CPPUNIT_ASSERT(q1.getScaledAngleRotation(0.25) == q2);
		CPPUNIT_ASSERT(q1 == q2.getScaledAngleRotation(4.0));
		CPPUNIT_ASSERT(q1.getScaledAngleRotation(0.5) == q3);
		CPPUNIT_ASSERT(q1 == q3.getScaledAngleRotation(2.0));

		CPPUNIT_ASSERT(q4 * q4 * q4 == q1);
	}

	void testSqrtRot() {
		int numTests = 1000;
		double halfAngle;
		Vector3<fpType> cay(0.0, 0.0, 0.0);
		Vector3<fpType> cayHalf = cay.getSqrtRotation();

		Quaternion<fpType> quat(0.0, 0.0, 0.0, 1.0);
		Quaternion<fpType> quatHalf = quat.getSqrtRotation();

		// Identity
		CPPUNIT_ASSERT(cay == cayHalf);
		CPPUNIT_ASSERT(quat == quatHalf);

		for (int i = 1; i < numTests; ++i) {
			halfAngle = -0.5 * M_PI + (M_PI * i) / numTests;
			cay.x = (rand() / (double) RAND_MAX);
			cay.y = (rand() / (double) RAND_MAX);
			cay.z = (rand() / (double) RAND_MAX);
			cay.normalize();

			quat.x = cay.x * sin(halfAngle);
			quat.y = cay.y * sin(halfAngle);
			quat.z = cay.z * sin(halfAngle);
			quat.w = (rand() % 2 ? 1.0 : -1.0) * cos(halfAngle);

			quatHalf = quat.getSqrtRotation();

			cay *= 2.0 * tan(halfAngle);
			cayHalf = cay.getSqrtRotation();

			if (compare<fpType>(halfAngle, 0.0) == 0) {
				CPPUNIT_ASSERT(compare<fpType>(cay.getNorm(), 0.0) == 0);
				CPPUNIT_ASSERT(compare<fpType>(cayHalf.getNorm(), 0.0) == 0);
			} else {
				CPPUNIT_ASSERT(cay.getUnit() == cayHalf.getUnit());
				CPPUNIT_ASSERT(quat.getAxis() == quatHalf.getAxis());
			}
			CPPUNIT_ASSERT(
					compare<fpType>(cay.getAngle(), 2.0 * cayHalf.getAngle())
							== 0);
			CPPUNIT_ASSERT(
					compare<fpType>(quat.getAngle(), 2.0 * quatHalf.getAngle(),
							PRECISION * 10) == 0);
			CPPUNIT_ASSERT(quat.compareRotations(quatHalf * quatHalf));

			quat = Quaternion<fpType>(cayHalf);
			CPPUNIT_ASSERT(Quaternion<fpType>(cay) == quat * quat);
		}
	}

	void testRotTowards() {
		fpType angle = 1.38692759366;
		Quaternion<fpType> id;
		Quaternion<fpType> q1(Vector3<fpType>(2.3, 4.5, 9.3), angle);
		Quaternion<fpType> q2(Vector3<fpType>(6.5, 4.3, 2.1), -0.3421);
		Quaternion<fpType> q3(Vector3<fpType>(1.5, 2.3, 3.1), 1.1111);
		Quaternion<fpType> trans = q1.getSqrtRotationTowards(q2);
		Quaternion<fpType> trans2 = q1.getSqrtRotationTowards(q3);
		Quaternion<fpType> trans3 = q2.getSqrtRotationTowards(q3);

		CPPUNIT_ASSERT(id.getSqrtRotationTowards(id) == id);
		CPPUNIT_ASSERT(q1.getSqrtRotationTowards(q1) == id);
		CPPUNIT_ASSERT(
				q1.getSqrtRotationTowards(id)
						== q1.inv().getScaledAngleRotation(0.5));
		CPPUNIT_ASSERT(
				id.getSqrtRotationTowards(q1)
						== q1.getScaledAngleRotation(0.5));
		CPPUNIT_ASSERT(trans == q2.getSqrtRotationTowards(q1).inv());

		CPPUNIT_ASSERT(q1 * trans * trans == q2);
		CPPUNIT_ASSERT(q2 * trans.inv() * trans.inv() == q1);

		CPPUNIT_ASSERT(q1 * trans2 * trans2 == q3);
		CPPUNIT_ASSERT(q3 * trans2.inv() * trans2.inv() == q1);

		CPPUNIT_ASSERT(q2 * trans3 * trans3 == q3);
		CPPUNIT_ASSERT(q3 * trans3.inv() * trans3.inv() == q2);
	}

	void testMatCols() {
		fpType a[3][3] = { { 0.1, 5.0, -3.0 }, //
				{ -5.0, 0.2, 2.0 }, //
				{ 3.0, -2.0, 0.3 } };

		Matrix3<fpType> m(a);
		Vector3<fpType> col0(a[0]);
		Vector3<fpType> col1(a[1]);
		Vector3<fpType> col2(a[2]);

		CPPUNIT_ASSERT(m.getColumn(0) == col0);
		CPPUNIT_ASSERT(m.getColumn(1) == col1);
		CPPUNIT_ASSERT(m.getColumn(2) == col2);
	}
};

CPPUNIT_TEST_SUITE_REGISTRATION (RotationTest);

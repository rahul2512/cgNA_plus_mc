/*
 * Copyright 2015 Jaroslaw Glowacki
 *
 * This file is part of cgDNArecon.
 *
 * cgDNArecon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cgDNArecon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cgDNArecon.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * A class to test parameter computation and shape reconstruction.
 */

#include <unistd.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <algebra3d/global.h>
#include <algebra3d/Matrix3.h>
#include <algebra3d/Vector3.h>
#include <algebra3d/Quaternion.h>

#include <cgDNArecon/reconstruct.h>

#include "Tester.h"

using namespace std;
using namespace algebra3d;
using namespace cgdna_recon;

class ReconstructTest: public CppUnit::TestFixture {

CPPUNIT_TEST_SUITE(ReconstructTest);
	CPPUNIT_TEST(testParamComp);
	CPPUNIT_TEST(testParamCompArr);
	CPPUNIT_TEST(testReconstr);
	CPPUNIT_TEST(testReconstrArr);
	CPPUNIT_TEST(testMultiple);
	CPPUNIT_TEST(testBegMidEnd);
	CPPUNIT_TEST(testBegMidEndBP);
	CPPUNIT_TEST(testBegMidEndArr);
	CPPUNIT_TEST(testGeneratedData);CPPUNIT_TEST_SUITE_END()
	;

protected:
	const static char * dataFileName;
	const static char * paramsFileName;

	int nBP;
	int nP;
	int extraPar;
	fpType extraVal;
	fpType *paramsRef;
	fpType *paramsC;
	fpType *paramsR;
	fpType *paramsBP;

	Vector3<fpType>* rBP;
	Quaternion<fpType> *dBP;
	Vector3<fpType>* rBP2;
	Quaternion<fpType> *dBP2;
	Vector3<fpType>* rBPM;
	Matrix3<fpType> *dBPM;
	Vector3<fpType>* rBPM2;
	Matrix3<fpType> *dBPM2;

	Vector3<fpType>* r1Ref;
	Vector3<fpType>* r2Ref;
	Quaternion<fpType> *d1Ref;
	Quaternion<fpType> *d2Ref;

	Vector3<fpType>* r1C;
	Vector3<fpType>* r2C;
	Quaternion<fpType> *d1C;
	Quaternion<fpType> *d2C;

	Vector3<fpType>* r1R;
	Vector3<fpType>* r2R;
	Quaternion<fpType> *d1R;
	Quaternion<fpType> *d2R;

	fpType **rBPArr;
	fpType ***dBPArr;
	fpType **rBPArr2;
	fpType ***dBPArr2;

	fpType **r1ArrC;
	fpType **r2ArrC;
	fpType ***d1ArrC;
	fpType ***d2ArrC;

	fpType **r1ArrR;
	fpType **r2ArrR;
	fpType ***d1ArrR;
	fpType ***d2ArrR;

public:
	void setUp() {
		extraPar = 20;
		extraVal = 12345.678;
		int tmp;
		int tmp2;
		fpType tmpFP;
		Matrix3<fpType> tmpM;
		Vector3<fpType> tmpV1;
		Vector3<fpType> tmpV2;
		Vector3<fpType> tmpV3;

		ifstream paramsStream;
		paramsStream.open(paramsFileName);
		ifstream dataStream(dataFileName);

		if (dataStream.fail() || paramsStream.fail()) {
			cerr << "Cannot read from files: '" << dataFileName << "', '"
					<< paramsFileName << "'" << endl;
			dataStream.close();
			paramsStream.close();
			_exit(1);
		}

		dataStream >> nBP;
		paramsStream >> nP;

		if (12 * nBP - 6 != nP) {
			cerr << "Testing data file mismatch: " << 12 * nBP - 6 << " != "
					<< nP << endl;
			_exit(2);
		}

		rBP = new Vector3<fpType> [nBP];
		dBP = new Quaternion<fpType> [nBP];
		rBP2 = new Vector3<fpType> [nBP];
		dBP2 = new Quaternion<fpType> [nBP];
		rBPM = new Vector3<fpType> [nBP];
		dBPM = new Matrix3<fpType> [nBP];
		rBPM2 = new Vector3<fpType> [nBP];
		dBPM2 = new Matrix3<fpType> [nBP];

		r1Ref = new Vector3<fpType> [nBP];
		r2Ref = new Vector3<fpType> [nBP];
		d1Ref = new Quaternion<fpType> [nBP];
		d2Ref = new Quaternion<fpType> [nBP];

		r1C = new Vector3<fpType> [nBP];
		r2C = new Vector3<fpType> [nBP];
		d1C = new Quaternion<fpType> [nBP];
		d2C = new Quaternion<fpType> [nBP];

		r1R = new Vector3<fpType> [nBP];
		r2R = new Vector3<fpType> [nBP];
		d1R = new Quaternion<fpType> [nBP];
		d2R = new Quaternion<fpType> [nBP];

		rBPArr = new fpType *[nBP];
		dBPArr = new fpType **[nBP];
		rBPArr2 = new fpType *[nBP];
		dBPArr2 = new fpType **[nBP];

		r1ArrC = new fpType *[nBP];
		r2ArrC = new fpType *[nBP];
		d1ArrC = new fpType **[nBP];
		d2ArrC = new fpType **[nBP];

		r1ArrR = new fpType *[nBP];
		r2ArrR = new fpType *[nBP];
		d1ArrR = new fpType **[nBP];
		d2ArrR = new fpType **[nBP];

		for (int i = 0; i < nBP; ++i) {
			rBPArr[i] = new fpType[3];
			dBPArr[i] = new fpType *[3];
			rBPArr2[i] = new fpType[3];
			dBPArr2[i] = new fpType *[3];

			r1ArrC[i] = new fpType[3];
			r2ArrC[i] = new fpType[3];
			d1ArrC[i] = new fpType *[3];
			d2ArrC[i] = new fpType *[3];

			r1ArrR[i] = new fpType[3];
			r2ArrR[i] = new fpType[3];
			d1ArrR[i] = new fpType *[3];
			d2ArrR[i] = new fpType *[3];

			for (int j = 0; j < 3; ++j) {
				dBPArr[i][j] = new fpType[3];
				dBPArr2[i][j] = new fpType[3];

				d1ArrC[i][j] = new fpType[3];
				d2ArrC[i][j] = new fpType[3];

				d1ArrR[i][j] = new fpType[3];
				d2ArrR[i][j] = new fpType[3];
			}
		}

		paramsRef = new fpType[nP + extraPar];
		paramsC = new fpType[nP + extraPar];
		paramsR = new fpType[nP + extraPar];
		paramsBP = new fpType[6 * nBP + extraPar];

		for (int i = 0; i < nBP; ++i) {
			dataStream >> tmp >> tmp2;
			dataStream >> r1Ref[i];
			dataStream >> tmpM;
			tmpM.transpose();
			// Orthogonalise the Matrix3
			tmpV1 = Vector3<fpType>(tmpM[0]).getUnit();
			tmpV2 = Vector3<fpType>(tmpM[1]).getUnit();
			tmpV3 = tmpV1.cross(tmpV2).getUnit();
			tmpV1 = tmpV2.cross(tmpV3).getUnit();
			d1Ref[i] = Matrix3<fpType>(tmpV1, tmpV2, tmpV3);
		}
		for (int i = 0; i < nBP; ++i) {
			dataStream >> tmp >> tmp2;
			dataStream >> r2Ref[i];
			dataStream >> tmpM;
			tmpM.transpose();
			// Orthogonalise the matrix
			tmpV1 = Vector3<fpType>(tmpM[0]).getUnit();
			tmpV2 = -Vector3<fpType>(tmpM[1]).getUnit();
			tmpV3 = tmpV1.cross(tmpV2).getUnit();
			tmpV1 = tmpV2.cross(tmpV3).getUnit();
			d2Ref[i] = Matrix3<fpType>(tmpV1, tmpV2, tmpV3);
		}
		for (int i = nP; i < nP + extraPar - 1; ++i) {
			paramsRef[i] = extraVal;
			paramsC[i] = extraVal;
			paramsR[i] = extraVal;
		}
		for (int i = 6 * nBP; i < 6 * nBP + extraPar - 1; ++i) {
			paramsBP[i] = extraVal;
		}
		for (int i = 0; i < nP; ++i) {
			paramsStream >> tmpFP;
			paramsRef[i] = tmpFP;
		}
		for (int i = 0; i < nBP; ++i) {
			for (int j = 0; j < 6; ++j) {
				paramsBP[j + 6 * i] = paramsRef[j + 12 * i + 6];
			}
		}
		paramsStream.close();
		dataStream.close();
	}

	void checkExtraParams() {
		for (int i = nP; i < nP + extraPar - 1; ++i) {
			CPPUNIT_ASSERT(paramsRef[i] == extraVal);
			CPPUNIT_ASSERT(paramsC[i] == extraVal);
			CPPUNIT_ASSERT(paramsR[i] == extraVal);
		}
		for (int i = 6 * nBP; i < 6 * nBP + extraPar - 1; ++i) {
			CPPUNIT_ASSERT(paramsBP[i] == extraVal);
		}
	}

	void tearDown() {
		for (int i = 0; i < nBP; ++i) {
			for (int j = 0; j < 3; ++j) {
				delete[] dBPArr[i][j];
				delete[] dBPArr2[i][j];

				delete[] d1ArrC[i][j];
				delete[] d2ArrC[i][j];

				delete[] d1ArrR[i][j];
				delete[] d2ArrR[i][j];
			}
			delete[] rBPArr[i];
			delete[] dBPArr[i];
			delete[] rBPArr2[i];
			delete[] dBPArr2[i];

			delete[] r1ArrC[i];
			delete[] r2ArrC[i];
			delete[] d1ArrC[i];
			delete[] d2ArrC[i];

			delete[] r1ArrR[i];
			delete[] r2ArrR[i];
			delete[] d1ArrR[i];
			delete[] d2ArrR[i];
		}
		delete[] rBP;
		delete[] dBP;

		delete[] r1Ref;
		delete[] r2Ref;
		delete[] d1Ref;
		delete[] d2Ref;

		delete[] rBPArr;
		delete[] dBPArr;
		delete[] rBPArr2;
		delete[] dBPArr2;

		delete[] r1ArrC;
		delete[] r2ArrC;
		delete[] d1ArrC;
		delete[] d2ArrC;

		delete[] r1ArrR;
		delete[] r2ArrR;
		delete[] d1ArrR;
		delete[] d2ArrR;

		delete[] r1C;
		delete[] r2C;
		delete[] d1C;
		delete[] d2C;

		delete[] r1R;
		delete[] r2R;
		delete[] d1R;
		delete[] d2R;

		delete[] paramsRef;
		delete[] paramsC;
		delete[] paramsR;
		delete[] paramsBP;
	}

	void copyToArr(Vector3<fpType> *r1, Vector3<fpType> *r2,
			Quaternion<fpType> *d1, Quaternion<fpType> *d2, fpType** r1Arr,
			fpType** r2Arr, fpType***d1Arr, fpType***d2Arr) {
		for (int i = 0; i < nBP; ++i) {
			r1[i].get(r1Arr[i]);
			r2[i].get(r2Arr[i]);
			Matrix3<fpType>(d1[i]).get(d1Arr[i][0], d1Arr[i][1], d1Arr[i][2]);
			Matrix3<fpType>(d2[i]).get(d2Arr[i][0], d2Arr[i][1], d2Arr[i][2]);
		}
	}

	void testParamComp() {
		memset(paramsC, 0, nP * sizeof(fpType));

		computeParams(nBP, r1Ref, r2Ref, d1Ref, d2Ref, paramsC);
		checkExtraParams();

		for (int i = 0; i < nP; i += 3) {
			CPPUNIT_ASSERT(
					Vector3<fpType>(paramsRef + i).compare(paramsC + i,
							5.0e-3));
		}
	}

	void testParamCompArr() {
		memset(paramsC, 0, nP * sizeof(fpType));

		copyToArr(r1Ref, r2Ref, d1Ref, d2Ref, r1ArrC, r2ArrC, d1ArrC, d2ArrC);

		computeParamsArr(nBP, r1ArrC, r2ArrC, d1ArrC, d2ArrC, paramsC);
		checkExtraParams();

		for (int i = 0; i < nP; i += 3) {
			CPPUNIT_ASSERT(
					Vector3<fpType>(paramsRef + i).compare(paramsC + i,
							5.0e-3));
		}
	}

	void testReconstr() {
		vector<Vector3<fpType> > rBPVec;
		vector<Quaternion<fpType> > dBPVec;

		vector<Vector3<fpType> > r1Vec;
		vector<Vector3<fpType> > r2Vec;
		vector<Quaternion<fpType> > d1Vec;
		vector<Quaternion<fpType> > d2Vec;
		vector<fpType> paramsVec;

		paramsVec.resize(nP + extraPar);
		for (int i = 0; i < nP; ++i) {
			paramsVec[i] = paramsRef[i];
		}
		for (int i = nP; i < nP + extraPar; ++i) {
			paramsVec[i] = 13.0;
		}

		reconstruct(nBP, paramsRef, r1R, r2R, d1R, d2R);
		checkExtraParams();
		reconstructVec(paramsVec, r1Vec, r2Vec, d1Vec, d2Vec);
		reconstructBP(nBP, paramsRef, rBP, dBP, true);
		checkExtraParams();
		reconstructBPM(nBP, paramsRef, rBPM, dBPM, true);
		checkExtraParams();

		for (int i = nP; i < nP + extraPar; ++i) {
			CPPUNIT_ASSERT(paramsVec[i] == 13.0);
		}
		reconstructBP(nBP, paramsBP, rBP2, dBP2, false);
		checkExtraParams();

		reconstructBPM(nBP, paramsBP, rBPM2, dBPM2, false);
		checkExtraParams();

		Quaternion<fpType> offRot = d2Ref[0]
				* d2Ref[0].getSqrtRotationTowards(d1Ref[0]); //d2Ref[0]
				//* (d2Ref[0].inv() * d1Ref[0]).getScaledAngleRotation(0.5);
		Vector3<fpType> offTran = (r1Ref[0] + r2Ref[0]) * 0.5;

		for (int i = 0; i < nBP; ++i) {
			d1Ref[i] = offRot.inv() * d1Ref[i];
			d2Ref[i] = offRot.inv() * d2Ref[i];
			r1Ref[i] = offRot.inv() * (r1Ref[i] - offTran);
			r2Ref[i] = offRot.inv() * (r2Ref[i] - offTran);

			CPPUNIT_ASSERT(r1Ref[i].compare(r1R[i], 5.0e-2));
			CPPUNIT_ASSERT(r2Ref[i].compare(r2R[i], 5.0e-2));
			CPPUNIT_ASSERT(Matrix3<fpType>(d1Ref[i]).compare(d1R[i], 5.0e-3));
			CPPUNIT_ASSERT(Matrix3<fpType>(d2Ref[i]).compare(d2R[i], 5.0e-3));

			CPPUNIT_ASSERT(r1Vec[i].compare(r1R[i]));
			CPPUNIT_ASSERT(r2Vec[i].compare(r2R[i]));
			CPPUNIT_ASSERT(d1Vec[i].compare(d1R[i]));
			CPPUNIT_ASSERT(d2Vec[i].compare(d2R[i]));

			CPPUNIT_ASSERT(rBP[i].compare((r1R[i] + r2R[i]) * 0.5));
			CPPUNIT_ASSERT(rBP[i] == rBP2[i]);
			CPPUNIT_ASSERT(rBP[i] == rBPM[i]);
			CPPUNIT_ASSERT(rBP[i] == rBPM2[i]);
			CPPUNIT_ASSERT(
					dBP[i] == d1R[i] * d1R[i].getSqrtRotationTowards(d2R[i]));
			CPPUNIT_ASSERT(dBP[i] == dBP2[i]);
			// The comparison done is this order to compare Matrices not
			// quaternions
			CPPUNIT_ASSERT(dBPM[i] == dBP[i]);
			CPPUNIT_ASSERT(dBPM2[i] == dBP[i]);
		}

		// This can be used as sample DNA data
//		Vector3<fpType> v;
//		Quaternion<fpType> q;
//		for (int i = 0; i < nB; ++i) {
//			v = r1Ref[i];
//			cout << "{" << v.x << ", " << v.y << ", " << v.z << "}, //" << endl;
//		}
//		for (int i = 0; i < nB; ++i) {
//			q = d1Ref[i];
//			cout << "{" << q.x << ", " << q.y << ", " << q.z << ", " << q.w
//					<< "}, //" << endl;
//		}
//		cout << endl << endl;
//		for (int i = 0; i < nB; ++i) {
//			v = r2Ref[i];
//			cout << "{" << v.x << ", " << v.y << ", " << v.z << "}, //" << endl;
//		}
//		for (int i = 0; i < nB; ++i) {
//			q = d2Ref[i];
//			cout << "{" << q.x << ", " << q.y << ", " << q.z << ", " << q.w
//					<< "}, //" << endl;
//		}
	}

	void testReconstrArr() {
		Vector3<fpType> rTmp;
		Quaternion<fpType> dTmp;

		reconstructArr(nBP, paramsRef, r1ArrR, r2ArrR, d1ArrR, d2ArrR);
		checkExtraParams();

		reconstructBPArr(nBP, paramsRef, rBPArr, dBPArr, true);
		checkExtraParams();

		reconstructBPArr(nBP, paramsBP, rBPArr2, dBPArr2, false);
		checkExtraParams();

		Quaternion<fpType> offRot = d2Ref[0]
				* ((d2Ref[0].inv() * d1Ref[0]).getScaledAngleRotation(0.5));
		Vector3<fpType> offTran = (r1Ref[0] + r2Ref[0]) * 0.5;

		for (int i = 0; i < nBP; ++i) {
			d1Ref[i] = offRot.inv() * d1Ref[i];
			d2Ref[i] = offRot.inv() * d2Ref[i];
			r1Ref[i] = offRot.inv() * (r1Ref[i] - offTran);
			r2Ref[i] = offRot.inv() * (r2Ref[i] - offTran);

			CPPUNIT_ASSERT(r1Ref[i].compare(r1ArrR[i], 5.0e-2));
			CPPUNIT_ASSERT(r2Ref[i].compare(r2ArrR[i], 5.0e-2));
			CPPUNIT_ASSERT(
					Matrix3<fpType>(d1Ref[i]).compare(
							Matrix3<fpType>(d1ArrR[i][0], //
									d1ArrR[i][1], //
									d1ArrR[i][2]), 5.0e-3));
			CPPUNIT_ASSERT(
					Matrix3<fpType>(d2Ref[i]).compare(
							Matrix3<fpType>(d2ArrR[i][0], //
									d2ArrR[i][1], //
									d2ArrR[i][2]), 5.0e-3));

			rTmp = r1ArrR[i];
			rTmp = (rTmp + r2ArrR[i]) * 0.5;

			dTmp = Matrix3<fpType>(d2ArrR[i][0], d2ArrR[i][1], d2ArrR[i][2]);
			dTmp = dTmp
					* dTmp.getSqrtRotationTowards(
							Matrix3<fpType>(d1ArrR[i][0], d1ArrR[i][1],
									d1ArrR[i][2]));

			CPPUNIT_ASSERT(rTmp.compare(rBPArr[i]));
			CPPUNIT_ASSERT(dTmp.compare(Matrix3<fpType>(dBPArr[i][0], //
					dBPArr[i][1], //
					dBPArr[i][2])));

			CPPUNIT_ASSERT(rTmp.compare(rBPArr2[i]));
			CPPUNIT_ASSERT(dTmp.compare(Matrix3<fpType>(dBPArr2[i][0], //
					dBPArr2[i][1], //
					dBPArr2[i][2])));
		}
	}

	void testMultiple() {
		memset(paramsC, 0, nP * sizeof(fpType));

		copyToArr(r1Ref, r2Ref, d1Ref, d2Ref, r1ArrC, r2ArrC, d1ArrC, d2ArrC);

		computeParamsArr(nBP, r1ArrC, r2ArrC, d1ArrC, d2ArrC, paramsC);

		checkExtraParams();
		reconstruct(nBP, paramsRef, r1R, r2R, d1R, d2R);
		checkExtraParams();
		computeParams(nBP, r1R, r2R, d1R, d2R, paramsR);
		checkExtraParams();

		for (int i = 0; i < nP; i += 3) {
			CPPUNIT_ASSERT(
					Vector3<fpType>(paramsR + i).compare(paramsC + i, 5.0e-3));
		}

		reconstructArr(nBP, paramsR, r1ArrR, r2ArrR, d1ArrR, d2ArrR);
		checkExtraParams();

		for (int i = 0; i < nBP; ++i) {
			CPPUNIT_ASSERT(r1R[i].compare(r1ArrR[i], PRECISION * 5));
			CPPUNIT_ASSERT(r2R[i].compare(r2ArrR[i], PRECISION * 5));
			CPPUNIT_ASSERT(
					Matrix3<fpType>(d1R[i]) == Matrix3<fpType>(d1ArrR[i][0], //
							d1ArrR[i][1], //
							d1ArrR[i][2]));
			CPPUNIT_ASSERT(
					Matrix3<fpType>(d2R[i]) == Matrix3<fpType>(d2ArrR[i][0], //
							d2ArrR[i][1], //
							d2ArrR[i][2]));
		}

		computeParamsArr(nBP, r1ArrR, r2ArrR, d1ArrR, d2ArrR, paramsC);
		checkExtraParams();

		for (int i = 0; i < nP; ++i) {
			CPPUNIT_ASSERT(
					compare < fpType > (paramsR[i], paramsC[i], PRECISION * 5) == 0);
		}
	}

	void testBegMidEnd() {
		Vector3<fpType> rBegin;
		Quaternion<fpType> dBegin;
		Vector3<fpType> rMid;
		Quaternion<fpType> dMid;
		Vector3<fpType> rEnd;
		Quaternion<fpType> dEnd;

		int m = (nBP % 2 == 0) ? nBP / 2 : (nBP + 1) / 2;
		int e = nBP - 1;

		reconstruct(nBP, paramsRef, r1R, r2R, d1R, d2R);
		checkExtraParams();
		reconstructFML(nBP, paramsRef, rBegin, dBegin, rMid, dMid, rEnd, dEnd,
				true);
		checkExtraParams();

		CPPUNIT_ASSERT(rBegin == 0.5 * (r1R[0] + r2R[0]));
		CPPUNIT_ASSERT(rEnd == 0.5 * (r1R[e] + r2R[e]));
		CPPUNIT_ASSERT(rMid == 0.5 * (r1R[m] + r2R[m]));

		CPPUNIT_ASSERT(
				dBegin.compareRotations(
						d1R[0] * d1R[0].getSqrtRotationTowards(d2R[0])));
		CPPUNIT_ASSERT(
				dMid.compareRotations(
						d1R[m] * d1R[m].getSqrtRotationTowards(d2R[m])));
		CPPUNIT_ASSERT(
				dEnd.compareRotations(
						d1R[e] * d1R[e].getSqrtRotationTowards(d2R[e])));
	}

	void testBegMidEndBP() {
		// Index of the middle BP
		int mid = (nBP % 2 == 0) ? nBP / 2 : (nBP + 1) / 2;
		Vector3<fpType> rBegin;
		Quaternion<fpType> dBegin;
		Vector3<fpType> rMid;
		Quaternion<fpType> dMid;
		Vector3<fpType> rEnd;
		Quaternion<fpType> dEnd;

		Vector3<fpType> rBeginBP;
		Quaternion<fpType> dBeginBP;
		Vector3<fpType> rMidBP;
		Quaternion<fpType> dMidBP;
		Vector3<fpType> rEndBP;
		Quaternion<fpType> dEndBP;

		reconstructBP(nBP, paramsRef, rBP, dBP, true);

		reconstructFML(nBP, paramsRef, rBegin, dBegin, rMid, dMid, rEnd, dEnd,
				true);
		checkExtraParams();
		reconstructFML(nBP, paramsBP, rBeginBP, dBeginBP, rMidBP, dMidBP,
				rEndBP, dEndBP, false);
		checkExtraParams();

		CPPUNIT_ASSERT(rBP[0] == rBegin);
		CPPUNIT_ASSERT(rBegin == rBeginBP);
		CPPUNIT_ASSERT(rBP[mid] == rMid);
		CPPUNIT_ASSERT(rMid == rMidBP);
		CPPUNIT_ASSERT(rBP[nBP - 1] == rEnd);
		CPPUNIT_ASSERT(rEnd == rEndBP);

		CPPUNIT_ASSERT(dBP[0].compareRotations(dBegin));
		CPPUNIT_ASSERT(dBegin.compareRotations(dBeginBP));
		CPPUNIT_ASSERT(dBP[mid].compareRotations(dMid));
		CPPUNIT_ASSERT(dMid.compareRotations(dMidBP));
		CPPUNIT_ASSERT(dBP[nBP - 1].compareRotations(dEnd));
		CPPUNIT_ASSERT(dEnd.compareRotations(dEndBP));
	}

	void testBegMidEndArr() {
		Vector3<fpType> rBegin;
		Quaternion<fpType> dBegin;
		Vector3<fpType> rMid;
		Quaternion<fpType> dMid;
		Vector3<fpType> rEnd;
		Quaternion<fpType> dEnd;

		int m = (nBP % 2 == 0) ? nBP / 2 : (nBP + 1) / 2;
		int e = nBP - 1;

		reconstruct(nBP, paramsRef, r1R, r2R, d1R, d2R);
		checkExtraParams();
		reconstructFMLArr(nBP, paramsRef, r1ArrR[0], d1ArrR[0], r1ArrR[1],
				d1ArrR[1], r1ArrR[2], d1ArrR[2], true);
		checkExtraParams();

		rBegin = r1ArrR[0];
		rMid = r1ArrR[1];
		rEnd = r1ArrR[2];

		dBegin = Matrix3<fpType>(d1ArrR[0][0], d1ArrR[0][1], d1ArrR[0][2]);
		dMid = Matrix3<fpType>(d1ArrR[1][0], d1ArrR[1][1], d1ArrR[1][2]);
		dEnd = Matrix3<fpType>(d1ArrR[2][0], d1ArrR[2][1], d1ArrR[2][2]);

		CPPUNIT_ASSERT(rBegin == 0.5 * (r1R[0] + r2R[0]));
		CPPUNIT_ASSERT(rMid == 0.5 * (r1R[m] + r2R[m]));
		CPPUNIT_ASSERT(rEnd == 0.5 * (r1R[e] + r2R[e]));

		CPPUNIT_ASSERT(
				dBegin.compareRotations(
						d1R[0] * d1R[0].getSqrtRotationTowards(d2R[0])));
		CPPUNIT_ASSERT(
				dMid.compareRotations(
						d1R[m] * d1R[m].getSqrtRotationTowards(d2R[m])));
		CPPUNIT_ASSERT(
				dEnd.compareRotations(
						d1R[e] * d1R[e].getSqrtRotationTowards(d2R[e])));

	}

	void testGeneratedData() {
#include "gen/testData.cpp"

		int n = 12 * nbp - 6;

		// See comment below
		fpType errorScale = 0.1;

		fpType shape[n];

		Vector3<fpType> zero(0.0, 0.0, 0.0);
		Vector3<fpType> r1[nbp];
		Vector3<fpType> r2[nbp];
		Quaternion<fpType> D1[nbp];
		Quaternion<fpType> D2[nbp];

		Vector3<fpType> r1m[nbp];
		Vector3<fpType> r2m[nbp];
		Matrix3<fpType> D1m[nbp];
		Matrix3<fpType> D2m[nbp];

		fpType **r1Arr = new fpType *[nbp];
		fpType **r2Arr = new fpType *[nbp];
		fpType ***D1Arr = new fpType **[nbp];
		fpType ***D2Arr = new fpType **[nbp];

		for (int i = 0; i < nbp; ++i) {
			r1Arr[i] = new fpType[3];
			r2Arr[i] = new fpType[3];
			D1Arr[i] = new fpType *[3];
			D2Arr[i] = new fpType *[3];
			for (int j = 0; j < 3; ++j) {
				D1Arr[i][j] = new fpType[3];
				D2Arr[i][j] = new fpType[3];
			}
		}

		for (int t = 0; t < numTests; ++t) {
			// Test reconstruction
			for (int i = 0; i < n; ++i) {
				shape[i] = shapesCorr[t][i];
			}

			reconstruct<fpType>(nbp, shape, r1, r2, D1, D2);
			reconstructM<fpType>(nbp, shape, r1m, r2m, D1m, D2m);
			reconstructArr<fpType>(nbp, shape, r1Arr, r2Arr, D1Arr, D2Arr);

			for (int i = 0; i < nbp; ++i) {
				Vector3<fpType> r1CorrV(r1Corr[t][i]);
				Vector3<fpType> r2CorrV(r2Corr[t][i]);

				Vector3<fpType> r1b = Vector3<fpType>(r1Arr[i]);
				Vector3<fpType> r2b = Vector3<fpType>(r2Arr[i]);

				Matrix3<fpType> m1Corr(D1Corr[t][i]);
				Matrix3<fpType> m2Corr(D2Corr[t][i]);

				Matrix3<fpType> m1 = D1[i];
				Matrix3<fpType> m2 = D2[i];

				Matrix3<fpType> m1b(D1Arr[i][0], D1Arr[i][1], D1Arr[i][2]);
				Matrix3<fpType> m2b(D2Arr[i][0], D2Arr[i][1], D2Arr[i][2]);

				// Comparisons are made up to PRECISION
				// The error might be grater for longer sequences
				// hence extra division
				CPPUNIT_ASSERT(
						(r1CorrV - r1[i]) / r1CorrV.getNorm() * errorScale
								== zero);
				CPPUNIT_ASSERT(
						(r2CorrV - r2[i]) / r2CorrV.getNorm() * errorScale
								== zero);

				CPPUNIT_ASSERT(
						(r1CorrV - r1b) / r1CorrV.getNorm() * errorScale
								== zero);
				CPPUNIT_ASSERT(
						(r2CorrV - r2b) / r2CorrV.getNorm() * errorScale
								== zero);

				CPPUNIT_ASSERT(
						(r1CorrV - r1m[i]) / r1CorrV.getNorm() * errorScale
								== zero);
				CPPUNIT_ASSERT(
						(r2CorrV - r2m[i]) / r2CorrV.getNorm() * errorScale
								== zero);

				CPPUNIT_ASSERT(m1Corr.compare(m1));
				CPPUNIT_ASSERT(m2Corr.compare(m2));

				CPPUNIT_ASSERT(m1Corr.compare(m1b));
				CPPUNIT_ASSERT(m2Corr.compare(m2b));

				CPPUNIT_ASSERT(m1Corr.compare(D1m[i]));
				CPPUNIT_ASSERT(m2Corr.compare(D2m[i]));
			}

			// Test param computation
			computeParams<fpType>(nbp, r1, r2, D1, D2, shape);

			for (int i = 0; i < n; ++i) {
				CPPUNIT_ASSERT(
						compare<fpType>(shapesCorr[t][i], shape[i], PRECISION * 50.0) == 0);
			}

			computeParamsArr<fpType>(nbp, r1Arr, r2Arr, D1Arr, D2Arr, shape);

			for (int i = 0; i < n; ++i) {
				CPPUNIT_ASSERT(
						compare<fpType>(shapesCorr[t][i], shape[i], PRECISION * 50.0) == 0);
			}
		}
	}
};

const char * ReconstructTest::dataFileName = "test/data.txt";
const char * ReconstructTest::paramsFileName = "test/params.txt";

CPPUNIT_TEST_SUITE_REGISTRATION(ReconstructTest);

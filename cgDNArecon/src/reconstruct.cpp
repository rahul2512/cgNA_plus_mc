/*
 * Copyright 2015 Jaroslaw Glowacki
 * jarek (dot) glowacki (at) gmail (dot) com
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
 * Definitions of the global functions.
 */
#include <vector>

#include <algebra3d/global.h>
#include <algebra3d/Matrix3.h>
#include <algebra3d/Vector3.h>
#include <algebra3d/Quaternion.h>

using namespace std;
using namespace algebra3d;

/**
 * @namespace cgdna_recon
 * @brief A namespace for the whole library.
 */
namespace cgdna_recon {

/****************************************************************************/
// A Helper function used for reconstruction
template<class fpType>
void computeParamsHelper(int i,
		fpType *params, //
		Vector3<fpType> &r11, Vector3<fpType> &r21, Vector3<fpType> &r12,
		Vector3<fpType> &r22, //
		Quaternion<fpType> &d11, Quaternion<fpType> &d21,
		Quaternion<fpType> &d12,
		Quaternion<fpType> &d22, //
		Quaternion<fpType> &bpTran, Quaternion<fpType> &bpQuat,
		Quaternion<fpType> &bpTran2, Quaternion<fpType> &bpQuat2,
		Vector3<fpType> &off, Vector3<fpType> &off2) {

	// Compute intra
	off = (r11 - r21);
	off = bpQuat.inv() * off;

	// Store params
	(Vector3<fpType>(bpTran) * 5.0).get(params + 3 * (4 * i));
	off.get(params + 3 * (4 * i + 1));

	// Compute inter
	bpTran2 = bpTran;
	bpQuat2 = bpQuat;
	bpTran = (d22.inv() * d12);
	bpQuat = d22 * bpTran.getSqrtRotation();
	// Junction frame wrt. lab
	bpTran2 = (bpQuat2.inv() * bpQuat);
	bpQuat2 = bpQuat2 * bpTran2.getSqrtRotation();
	// Base pair positions wrt. lab frame
	off = (r11 + r21) * 0.5;
	off2 = (r12 + r22) * 0.5;
	off = bpQuat2.inv() * (off2 - off);
	(Vector3<fpType>(bpTran2) * 5.0).get(params + 3 * (4 * i + 2));
	off.get(params + 3 * (4 * i + 3));
}

/****************************************************************************/
template<class fpType>
void computeParams(int n, Vector3<fpType>* r1, Vector3<fpType>* r2,
		Quaternion<fpType> *d1, Quaternion<fpType> *d2, fpType *params) {
	Quaternion<fpType> bpTran;
	Quaternion<fpType> bpTran2;
	Quaternion<fpType> bpQuat;
	Quaternion<fpType> bpQuat2;
	Vector3<fpType> off;
	Vector3<fpType> off2;

	// Compute parameters
	bpTran = (d2[0].inv() * d1[0]);
	bpQuat = d2[0] * bpTran.getSqrtRotation();
	for (int i = 0; i < n - 1; ++i) {
		computeParamsHelper(i, params, //
				r1[i], r2[i], r1[i + 1], r2[i + 1], //
				d1[i], d2[i], d1[i + 1], d2[i + 1], //
				bpTran, bpQuat, //
				bpTran2, bpQuat2, //
				off, off2);
	}
	// Compute last set of intra
	off = (r1[n - 1] - r2[n - 1]);
	off = bpQuat.inv() * off;
	// Store last set of instra
	(Vector3<fpType>(bpTran) * 5.0).get(params + 3 * (4 * (n - 1)));
	off.get(params + 3 * (4 * (n - 1) + 1));
}

// Define the specializations that will be usable
template
void computeParams<double>(int n, Vector3<double>* r1, Vector3<double>* r2,
		Quaternion<double> *d1, Quaternion<double> *d2, double *params);
template
void computeParams<float>(int n, Vector3<float>* r1, Vector3<float>* r2,
		Quaternion<float> *d1, Quaternion<float> *d2, float *params);

/****************************************************************************/
template<class fpType>
void computeParamsArr(int n, fpType** r1, fpType** r2, fpType***d1, fpType***d2,
		fpType *params) {
	Quaternion<fpType> bpTran;
	Quaternion<fpType> bpTran2;
	Quaternion<fpType> bpQuat;
	Quaternion<fpType> bpQuat2;
	Quaternion<fpType> d1tmp;
	Quaternion<fpType> d2tmp;
	Quaternion<fpType> d1tmp2;
	Quaternion<fpType> d2tmp2;

	Vector3<fpType> off;
	Vector3<fpType> off2;
	Vector3<fpType> r1tmp;
	Vector3<fpType> r2tmp;
	Vector3<fpType> r1tmp2;
	Vector3<fpType> r2tmp2;

	// Compute parameters
	d1tmp = Matrix3<fpType>(d1[0][0], d1[0][1], d1[0][2]);
	d2tmp = Matrix3<fpType>(d2[0][0], d2[0][1], d2[0][2]);
	r1tmp = Vector3<fpType>(r1[0]);
	r2tmp = Vector3<fpType>(r2[0]);
	bpTran = (d2tmp.inv() * d1tmp);
	bpQuat = d2tmp * bpTran.getSqrtRotation();
	for (int i = 0; i < n - 1; ++i) {
		d1tmp2 = Matrix3<fpType>(d1[i + 1][0], d1[i + 1][1], d1[i + 1][2]);
		d2tmp2 = Matrix3<fpType>(d2[i + 1][0], d2[i + 1][1], d2[i + 1][2]);
		r1tmp2 = Vector3<fpType>(r1[i + 1]);
		r2tmp2 = Vector3<fpType>(r2[i + 1]);

		// Compute intra
		off = (r1tmp - r2tmp);
		off = bpQuat.inv() * off;

		// Store params
		(Vector3<fpType>(bpTran) * 5.0).get(params + 3 * (4 * i));
		off.get(params + 3 * (4 * i + 1));

		// Compute inter
		bpTran2 = bpTran;
		bpQuat2 = bpQuat;
		bpTran = (d2tmp2.inv() * d1tmp2);
		bpQuat = d2tmp2 * bpTran.getSqrtRotation();
		// Junction frame wrt. lab
		bpTran2 = (bpQuat2.inv() * bpQuat);
		bpQuat2 = bpQuat2 * bpTran2.getSqrtRotation();
		// Base pair positions wrt. lab frame
		off = (r1tmp + r2tmp) * 0.5;
		off2 = (r1tmp2 + r2tmp2) * 0.5;
		off = bpQuat2.inv() * (off2 - off);
		(Vector3<fpType>(bpTran2) * 5.0).get(params + 3 * (4 * i + 2));
		off.get(params + 3 * (4 * i + 3));

		d1tmp = d1tmp2;
		d2tmp = d2tmp2;
		r1tmp = r1tmp2;
		r2tmp = r2tmp2;
	}
	// Compute last set of intra
	off = (Vector3<fpType>(r1[n - 1]) - Vector3<fpType>(r2[n - 1]));
	off = bpQuat.inv() * off;
	// Store last set of instra
	(Vector3<fpType>(bpTran) * 5.0).get(params + 3 * (4 * (n - 1)));
	off.get(params + 3 * (4 * (n - 1) + 1));
}

// Define the specializations that will be usable
template
void computeParamsArr<double>(int n, double** r1, double** r2, double ***d1,
		double ***d2, double *params);
template
void computeParamsArr<float>(int n, float** r1, float** r2, float ***d1,
		float ***d2, float *params);

/****************************************************************************/
// A part of the helper functionality for reconstruction
template<class fpType>
void nextBasePair(int i, const fpType *params, //
		Vector3<fpType> &q, Quaternion<fpType> &g, bool baseModel) {
	Quaternion<fpType> L;
	Vector3<fpType> zeta;
	Quaternion<fpType> h;

	int firstIndex = 0;
	int secondIndex = 0;
	if (baseModel) {
		firstIndex = 12 * i + 6;
		secondIndex = 12 * i + 9;
	} else {
		firstIndex = 6 * i;
		secondIndex = 6 * i + 3;
	}

	L = Vector3<fpType>(params + firstIndex) * 0.2;
	h = g * L.getSqrtRotation();
	g = g * L;

	zeta = Vector3<fpType>(params + secondIndex);
	q = q + h * zeta;
}

/****************************************************************************/
// A part of the helper functionality for reconstruction
template<class fpType>
void nextBasePairM(int i, const fpType *params, //
		Vector3<fpType> &q, Matrix3<fpType> &g, bool baseModel) {
	Matrix3<fpType> sqrtL;
	Vector3<fpType> cay;
	Vector3<fpType> zeta;
	Matrix3<fpType> h;

	int firstIndex = 0;
	int secondIndex = 0;
	if (baseModel) {
		firstIndex = 12 * i + 6;
		secondIndex = 12 * i + 9;
	} else {
		firstIndex = 6 * i;
		secondIndex = 6 * i + 3;
	}

	cay = Vector3<fpType>(params + firstIndex) * 0.2;
	sqrtL = cay.getSqrtRotation();
	h = g * sqrtL;
	g = h * sqrtL;

	zeta = Vector3<fpType>(params + secondIndex);
	q = q + h * zeta;
}

/****************************************************************************/
// A helper function used for reconstruction
template<class fpType>
void reconstructHelper(int i,
		const fpType *params, //
		Vector3<fpType> &r1, Vector3<fpType> &r2, Quaternion<fpType> &d1,
		Quaternion<fpType> &d2, //
		Vector3<fpType> &q, Quaternion<fpType> &g, bool last = false) {
	Vector3<fpType> cay = Vector3<fpType>(params + 12 * i) * 0.2;
	Quaternion<fpType> sqrtLambda = cay.getSqrtRotation();

	Vector3<fpType> halfXiRot;

	// Compute parameters
	d1 = g * sqrtLambda;
	d2 = g * sqrtLambda.inv();
	halfXiRot = g * Vector3<fpType>(params + 12 * i + 3) * 0.5;
	r1 = q + halfXiRot;
	r2 = q - halfXiRot;

	// Prepare bp frame for next iteration
	if (!last) {
		nextBasePair(i, params, q, g, true);
	}
}

/****************************************************************************/
// A helper function used for reconstruction
template<class fpType>
void reconstructHelperM(int i,
		const fpType *params, //
		Vector3<fpType> &r1, Vector3<fpType> &r2, Matrix3<fpType> &d1,
		Matrix3<fpType> &d2, //
		Vector3<fpType> &q, Matrix3<fpType> &g, bool last = false) {
	Vector3<fpType> cay = Vector3<fpType>(params + 12 * i) * 0.2;
	Matrix3<fpType> sqrtLambda = cay.getSqrtRotation();

	Vector3<fpType> halfXiRot;

	// Compute parameters
	d1 = g * sqrtLambda;
	d2 = g * sqrtLambda.getTranspose();
	halfXiRot = g * Vector3<fpType>(params + 12 * i + 3) * 0.5;
	r1 = q + halfXiRot;
	r2 = q - halfXiRot;

	// Prepare bp frame for next iteration
	if (!last) {
		nextBasePairM(i, params, q, g, true);
	}
}

/****************************************************************************/
template<class fpType>
void reconstruct(int n, const fpType *params, Vector3<fpType>* r1,
		Vector3<fpType>* r2, Quaternion<fpType> *d1, Quaternion<fpType> *d2) {
	Quaternion<fpType> g;

	Vector3<fpType> q;

	int i = 0;
	for (; i < n - 1; ++i) {
		reconstructHelper(i, params, //
				r1[i], r2[i], d1[i], d2[i], //
				q, g);
	}
	// In the last step the base pair data for next iteration
	// shouldn't be computed (this would read from outside of memory)
	reconstructHelper(i, params, //
			r1[i], r2[i], d1[i], d2[i], //
			q, g, true);
}

// Define the specializations that will be usable
template
void reconstruct<double>(int n, const double *params, Vector3<double>* r1,
		Vector3<double>* r2, Quaternion<double> *d1, Quaternion<double> *d2);
template
void reconstruct<float>(int n, const float *params, Vector3<float>* r1,
		Vector3<float>* r2, Quaternion<float> *d1, Quaternion<float> *d2);

/****************************************************************************/
template<class fpType>
void reconstructM(int n, const fpType *params, Vector3<fpType>* r1,
		Vector3<fpType>* r2, Matrix3<fpType> *d1, Matrix3<fpType> *d2) {
	Matrix3<fpType> g(true);

	Vector3<fpType> q;

	int i = 0;
	for (; i < n - 1; ++i) {
		reconstructHelperM(i, params, //
				r1[i], r2[i], d1[i], d2[i], //
				q, g);
	}
	// In the last step the base pair data for next iteration
	// shouldn't be computed (this would read from outside of memory)
	reconstructHelperM(i, params, //
			r1[i], r2[i], d1[i], d2[i], //
			q, g, true);
}

// Define the specializations that will be usable
template
void reconstructM<double>(int n, const double *params, Vector3<double>* r1,
		Vector3<double>* r2, Matrix3<double> *d1, Matrix3<double> *d2);
template
void reconstructM<float>(int n, const float *params, Vector3<float>* r1,
		Vector3<float>* r2, Matrix3<float> *d1, Matrix3<float> *d2);

/****************************************************************************/
template<class fpType>
void reconstructVec(const vector<fpType> &params, vector<Vector3<fpType> > &r1,
		vector<Vector3<fpType> > &r2, vector<Quaternion<fpType> > &d1,
		vector<Quaternion<fpType> > &d2) {
	int n = (params.size() + 6) / 12;

	r1.resize(n);
	d1.resize(n);
	r2.resize(n);
	d2.resize(n);

	reconstruct(n, &(params[0]), &(r1[0]), &(r2[0]), &(d1[0]), &(d2[0]));
}

// Define the specializations that will be usable
template
void reconstructVec<double>(const vector<double> &params,
		vector<Vector3<double> > &r1, vector<Vector3<double> > &r2,
		vector<Quaternion<double> > &d1, vector<Quaternion<double> > &d2);
template
void reconstructVec<float>(const vector<float> &params,
		vector<Vector3<float> > &r1, vector<Vector3<float> > &r2,
		vector<Quaternion<float> > &d1, vector<Quaternion<float> > &d2);

/****************************************************************************/
template<class fpType>
void reconstructMVec(const vector<fpType> &params, vector<Vector3<fpType> > &r1,
		vector<Vector3<fpType> > &r2, vector<Matrix3<fpType> > &d1,
		vector<Matrix3<fpType> > &d2) {
	int n = (params.size() + 6) / 12;

	r1.resize(n);
	d1.resize(n);
	r2.resize(n);
	d2.resize(n);

	reconstructM(n, &(params[0]), &(r1[0]), &(r2[0]), &(d1[0]), &(d2[0]));
}

// Define the specializations that will be usable
template
void reconstructMVec<double>(const vector<double> &params,
		vector<Vector3<double> > &r1, vector<Vector3<double> > &r2,
		vector<Matrix3<double> > &d1, vector<Matrix3<double> > &d2);
template
void reconstructMVec<float>(const vector<float> &params,
		vector<Vector3<float> > &r1, vector<Vector3<float> > &r2,
		vector<Matrix3<float> > &d1, vector<Matrix3<float> > &d2);

/****************************************************************************/
template<class fpType>
void reconstructArr(int n, fpType *params, fpType** r1, fpType** r2,
		fpType***d1, fpType***d2) {
	Quaternion<fpType> g;
	Quaternion<fpType> d1tmp;
	Quaternion<fpType> d2tmp;

	Vector3<fpType> q;
	Vector3<fpType> r1tmp;
	Vector3<fpType> r2tmp;

	int i = 0;
	for (; i < n - 1; ++i) {
		// Compute parameters
		reconstructHelper(i, params, //
				r1tmp, r2tmp, d1tmp, d2tmp, //
				q, g);

		// Store the data
		((Matrix3<fpType> ) d1tmp).get(d1[i][0], d1[i][1], d1[i][2]);
		((Matrix3<fpType> ) d2tmp).get(d2[i][0], d2[i][1], d2[i][2]);
		r1tmp.get(r1[i]);
		r2tmp.get(r2[i]);
	}
	// In the last step the base pair data for next iteration
	// shouldn't be computed (this would reed from outside of memory)
	reconstructHelper(i, params, //
			r1tmp, r2tmp, d1tmp, d2tmp, //
			q, g, true);

	// Store the data for the last step
	((Matrix3<fpType> ) d1tmp).get(d1[i][0], d1[i][1], d1[i][2]);
	((Matrix3<fpType> ) d2tmp).get(d2[i][0], d2[i][1], d2[i][2]);
	r1tmp.get(r1[i]);
	r2tmp.get(r2[i]);
}

// Define the specializations that will be usable
template
void reconstructArr<double>(int n, double *params, double** r1, double** r2,
		double***d1, double***d2);
template
void reconstructArr<float>(int n, float *params, float** r1, float** r2,
		float***d1, float***d2);

/****************************************************************************/
template<class fpType>
void reconstructBP(int n, const fpType *params, algebra3d::Vector3<fpType>* r,
		algebra3d::Quaternion<fpType> *d, bool baseModel) {
	r[0] = Vector3<fpType>(0.0, 0.0, 0.0);
	d[0] = Quaternion<fpType>(0.0, 0.0, 0.0, 1.0);

	// Compute data for base pairs
	for (int i = 1; i < n; ++i) {
		r[i] = r[i - 1];
		d[i] = d[i - 1];
		nextBasePair(i - 1, params, r[i], d[i], baseModel);
	}
}

// Define the specializations that will be usable
template
void reconstructBP<double>(int n, const double *params,
		algebra3d::Vector3<double>* r, algebra3d::Quaternion<double> *d,
		bool baseModel);
template
void reconstructBP<float>(int n, const float *params,
		algebra3d::Vector3<float>* r, algebra3d::Quaternion<float> *d,
		bool baseModel);

/****************************************************************************/
template<class fpType>
void reconstructBPM(int n, const fpType *params, algebra3d::Vector3<fpType>* r,
		algebra3d::Matrix3<fpType> *d, bool baseModel) {
	r[0] = Vector3<fpType>(0.0, 0.0, 0.0);
	d[0] = Matrix3<fpType>(true);

	// Compute data for base pairs
	for (int i = 1; i < n; ++i) {
		r[i] = r[i - 1];
		d[i] = d[i - 1];
		nextBasePairM(i - 1, params, r[i], d[i], baseModel);
	}
}

// Define the specializations that will be usable
template
void reconstructBPM<double>(int n, const double *params,
		algebra3d::Vector3<double>* r, algebra3d::Matrix3<double> *d,
		bool baseModel);
template
void reconstructBPM<float>(int n, const float *params,
		algebra3d::Vector3<float>* r, algebra3d::Matrix3<float> *d,
		bool baseModel);

/****************************************************************************/
template<class fpType>
void reconstructBPVec(const std::vector<fpType> &params,
		std::vector<algebra3d::Vector3<fpType> > &r,
		std::vector<algebra3d::Quaternion<fpType> > &d, bool baseModel) {
	int n = (baseModel) ? params.size() / 12 + 1 : params.size() / 6 + 1;

	r.resize(n);
	d.resize(n);

	reconstructBP(n, &(params[0]), &(r[0]), &(d[0]), baseModel);
}

// Define the specializations that will be usable
template
void reconstructBPVec<double>(const vector<double> &params,
		vector<Vector3<double> > &r, vector<Quaternion<double> > &d,
		bool baseModel);
template
void reconstructBPVec<float>(const vector<float> &params,
		vector<Vector3<float> > &r, vector<Quaternion<float> > &d,
		bool baseModel);

/****************************************************************************/
template<class fpType>
void reconstructBPMVec(const std::vector<fpType> &params,
		std::vector<algebra3d::Vector3<fpType> > &r,
		std::vector<algebra3d::Matrix3<fpType> > &d, bool baseModel) {
	int n = (baseModel) ? params.size() / 12 + 1 : params.size() / 6 + 1;

	r.resize(n);
	d.resize(n);

	reconstructBPM(n, &(params[0]), &(r[0]), &(d[0]), baseModel);
}

// Define the specializations that will be usable
template
void reconstructBPMVec<double>(const vector<double> &params,
		vector<Vector3<double> > &r, vector<Matrix3<double> > &d,
		bool baseModel);
template
void reconstructBPMVec<float>(const vector<float> &params,
		vector<Vector3<float> > &r, vector<Matrix3<float> > &d,
		bool baseModel);

/****************************************************************************/
template<class fpType>
void reconstructBPArr(int n, fpType *params, fpType** r, fpType***d,
		bool baseModel) {
	Vector3<fpType> pos = Vector3<fpType>(0.0, 0.0, 0.0);
	Quaternion<fpType> rot = Quaternion<fpType>(0.0, 0.0, 0.0, 1.0);

	// Compute data for base pairs
	for (int i = 0; i < n; ++i) {
		pos.get(r[i]);
		((Matrix3<fpType> ) rot).get(d[i][0], d[i][1], d[i][2]);
		nextBasePair(i, params, pos, rot, baseModel);
	}
}

// Define the specializations that will be usable
template
void reconstructBPArr<double>(int n, double *params, double** r, double***d,
		bool baseModel);
template
void reconstructBPArr<float>(int n, float *params, float** r, float***d,
		bool baseModel);

/****************************************************************************/
template<class fpType>
void reconstructFML(int n, fpType *params, //
		Vector3<fpType> &rFirst, Quaternion<fpType> &dFirst, //
		Vector3<fpType> &rMid, Quaternion<fpType> &dMid, //
		Vector3<fpType> &rLast, Quaternion<fpType> &dLast, //
		bool baseModel) {
	Quaternion<fpType> Lambda;
	Quaternion<fpType> sqrtLambda;
	Quaternion<fpType> g;

	Vector3<fpType> halfXiRot;
	Vector3<fpType> q;
	Vector3<fpType> zeta;

	rFirst = Vector3<fpType>(0.0, 0.0, 0.0);
	dFirst = Quaternion<fpType>(0.0, 0.0, 0.0, 1.0);

	int mid = (n % 2 == 0) ? n / 2 : (n + 1) / 2;
	int i = 0;

	// Initialize the mid base frame data for iteration
	rMid = rFirst;
	dMid = dFirst;

	// Compute data for base pairs until the mid'th one
	for (; i < mid; ++i) {
		nextBasePair(i, params, rMid, dMid, baseModel);
	}

	// Initialize the last base frame data for iteration
	rLast = rMid;
	dLast = dMid;

	// Compute data for base pairs until the end
	for (; i < n - 1; ++i) {
		nextBasePair(i, params, rLast, dLast, baseModel);
	}
}

// Define the specializations that will be usable
template
void reconstructFML<double>(int n, double *params, //
		Vector3<double> &rFirst, Quaternion<double> &dFirst, //
		Vector3<double> &rMid, Quaternion<double> &dMid, //
		Vector3<double> &rLast, Quaternion<double> &dLast, //
		bool baseModel);
template
void reconstructFML<float>(int n, float *params, //
		Vector3<float> &rFirst, Quaternion<float> &dFirst, //
		Vector3<float> &rMid, Quaternion<float> &dMid, //
		Vector3<float> &rLast, Quaternion<float> &dLast, //
		bool baseModel);

/****************************************************************************/
template<class fpType>
void reconstructFMLArr(int n, fpType *params, //
		fpType *rFirst, //
		fpType **dFirst, //
		fpType *rMid, //
		fpType **dMid, //
		fpType *rLast, //
		fpType **dLast, //
		bool baseModel) {

	Vector3<fpType> rFirstV;
	Quaternion<fpType> dFirstQ;
	Vector3<fpType> rMidV;
	Quaternion<fpType> dMidQ;
	Vector3<fpType> rLastV;
	Quaternion<fpType> dLastQ;

	reconstructFML(n, params, rFirstV, dFirstQ, rMidV, dMidQ, rLastV, dLastQ,
			baseModel);

	rFirstV.get(rFirst);
	((Matrix3<fpType> ) dFirstQ).get(dFirst[0], dFirst[1], dFirst[2]);
	rMidV.get(rMid);
	((Matrix3<fpType> ) dMidQ).get(dMid[0], dMid[1], dMid[2]);
	rLastV.get(rLast);
	((Matrix3<fpType> ) dLastQ).get(dLast[0], dLast[1], dLast[2]);
}

// Define the specializations that will be usable
template
void reconstructFMLArr<double>(int n, double *params, //
		double *rFirst, //
		double **dFirst, //
		double *rMid, //
		double **dMid, //
		double *rLast, //
		double **dLast, //
		bool baseModel);
template
void reconstructFMLArr<float>(int n, float *params, //
		float *rFirst, //
		float **dFirst, //
		float *rMid, //
		float **dMid, //
		float *rLast, //
		float **dLast, //
		bool baseModel);

} // namespace

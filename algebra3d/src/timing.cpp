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
 * @file
 * A simple binary for performance tests
 */

#include<cstdlib>

#include<sys/time.h>

#include<iostream>
#include<iomanip>
#include<vector>

#include<algebra3d/global.h>
#include<algebra3d/Matrix3.h>
#include<algebra3d/Vector3.h>
#include<algebra3d/Quaternion.h>

using namespace std;
using namespace algebra3d;

// Make it possible to decide at compile time which floating point type to use
#ifdef FPTYPE_FLOAT
typedef float fpType;
#else
typedef double fpType;
#endif

/**
 * Returns the number of microseconds from Epoch. This is a simple way of
 * calculating elapsed time with microsecond resolution.
 *
 * @return The number of microseconds from Epoch.
 */
unsigned long long getTime() {
	struct timeval time;
	unsigned long long t;

	gettimeofday(&time, NULL);
	t = time.tv_sec;
	t *= 1000000;
	t += time.tv_usec;

	return t;
}

// A macro for the repetitive code of each timing run
#define RUN(CODE, LABEL) \
	minTime = 1LL << 60; \
	for (int t = 0; t < numTrials; ++t) { \
		runningTime = getTime(); \
		for (int r = 0; r < numReps; ++r) { \
			for (int i = 0; i < numTests; ++i) { \
				CODE \
			} \
		} \
		runningTime = getTime() - runningTime; \
		if (minTime > runningTime) { \
			minTime = runningTime; \
		} \
	} \
	cout << setw(colWidth) << LABEL << setw(timeWidth) \
			<< minTime / 1000000.0 << " s" << endl;

/**
 * The main function that runs the simple performance tests.
 *
 * @param argc The number of command line arguments
 * @param argv The command line arguments
 * @return 0 if execution is correct, a positive integer on errors.
 */
int main(int argc, char **argv) {
	int numTests;
	int numReps;
	int numTrials;
	int colWidth = 50;
	int timeWidth = 10;

	unsigned long long runningTime;
	unsigned long long minTime;
	Quaternion<fpType> quat;
	Quaternion<fpType> quatHalf;
	Vector3<fpType> quatAx;

	Vector3<fpType> v;
	Vector3<fpType> vInit(1, 1, 1);
	Vector3<fpType> cayHalf;

	Matrix3<fpType> eye(true);
	Matrix3<fpType> mat;

	vector<Vector3<fpType> > cays;
	vector<Quaternion<fpType> > quats;
	vector<Matrix3<fpType> > mats;

	// Check if all required arguments are given
	if (argc != 4) {
		cerr << "Usage:\n"
				"    " << argv[0]
				<< " num_tests num_reps num_trials\n"
						"Where:\n"
						"    num_tests is the number of randomly generated rotations to use\n"
						"    num_reps  is the number of repetitions each num_tests tests are run\n"
						"    num_trials  is the number of trials to perform to choose the fastest result\n"
				<< endl;
		return 1;
	}

	numTests = atoi(argv[1]);
	numReps = atoi(argv[2]);
	numTrials = atoi(argv[3]);

	if (sizeof(fpType) == sizeof(float)) {
		cout << "Single precission test" << endl;
	} else {
		cout << "Double precission test" << endl;
	}

	// Allocate memory
	cays.resize(numTests);
	quats.resize(numTests);
	mats.resize(numTests);

	cout << "Genetaring testing data (" << numTests << " random rotations)"
			<< endl;
	cout.flush();

	// Initialize random number generator for reproducibility
	srand(12345678);

	for (int i = 0; i < numTests; ++i) {
		quats[i].x = (rand() / (double) RAND_MAX);
		quats[i].y = (rand() / (double) RAND_MAX);
		quats[i].z = (rand() / (double) RAND_MAX);
		quats[i].w = (rand() / (double) RAND_MAX);
		quats[i].normalize();

		cays[i] = quats[i];
		mats[i] = quats[i];
	}

	cout << "\nTiming conversions between representations of rotations\n"
			<< numReps << " times for all " << numTests << " rotations "
			<< endl;

#define Q2M mat = quats[i];
	RUN(Q2M, "Quaternion to Matrix: ")

#define Q2V v = quats[i];
	RUN(Q2V, "Quaternion to Vector: ")

	cout << endl;

#define M2V v = mats[i];
	RUN(M2V, "Matrix to Vector: ")

#define M2Q quat = mats[i];
	RUN(M2Q, "Matrix to Quaternion: ")

	cout << endl;

#define V2Q quat = cays[i];
	RUN(V2Q, "Vector to Quaternion: ")

#define V2Q2M mat = cays[i];
	RUN(V2Q2M, "Vector to Quaternion to Matrix: ")

#define V2M \
		fpType tana2 = cays[i].dot(cays[i]); \
		Matrix3<fpType> m(cays[i].crossMatrix()); \
		mat = m * m; \
		mat *= 0.5; \
		mat += m; \
		mat *= (fpType) (4.0 / (4.0 + tana2)); \
		mat[0][0] += 1.0; \
		mat[1][1] += 1.0; \
		mat[2][2] += 1.0;
	RUN(V2M, "direct Vector to Matrix: ")

	cout << "\nTiming different methods of applying rotation to a vector\n"
			<< numReps << " times for all " << numTests << " rotations "
			<< endl;

	v = vInit;
#define Q2MXV v = quats[i] * v;
	RUN(Q2MXV, "Quaternion x vector (through a Matrix): ")

	v = vInit;
#define QXV \
		quat[0] = v[0]; \
		quat[1] = v[1]; \
		quat[2] = v[2]; \
		quat[3] = 0.0; \
		quat = quats[i] * quat * quats[i].inv(); \
		v[0] = quat[0]; \
		v[1] = quat[1]; \
		v[2] = quat[2];
	RUN(QXV, "Quaternion x vector (tripple product): ")

	v = vInit;
#define QXV2 \
		quatAx.x = quats[i].x; \
		quatAx.y = quats[i].y; \
		quatAx.z = quats[i].z; \
		v = (quats[i].w * quats[i].w + quatAx.dot(quatAx)) * v \
			+ 2.0 * (v.dot(quatAx)) * quatAx \
			+ 2.0 * quats[i].w * (quatAx.cross(v));
	RUN(QXV2, "Quaternion x vector (tripple product fast): ")

	v = vInit;
#define MXV v = mats[i] * v;
	RUN(MXV, "Matrix x vector: ")

	cout << "\nTiming different methods of multiplication\n" << numReps
			<< " times for all " << numTests << " rotations " << endl;

#define QXQ quats[i] = quats[i] * quats[(i + 1) % numTests];
	RUN(QXQ, "Quaternion mult: ")

#define MXM mats[i] = mats[i] * mats[(i + 1) % numTests];
	RUN(MXM, "Matrix mult: ")

	cout << "\nTiming different methods of computing half rotation\n" << numReps
			<< " times for all " << numTests << " rotations " << endl;

#define VHALF cayHalf = cays[i].getSqrtRotation();
	RUN(VHALF, "Cayley sqrt: ")

#define QHALF quatHalf = quats[i].getSqrtRotation();
	RUN(QHALF, "Quaterion sqrt: ")

#define QSCALE quatHalf = quats[i].getScaledAngleRotation(0.5);
	RUN(QSCALE, "Quaterion scale angle: ")

	return 0;
}


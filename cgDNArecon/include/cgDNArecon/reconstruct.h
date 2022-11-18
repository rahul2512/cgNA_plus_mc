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

/**
 * @file
 * Declarations of all the global functions.
 */

#ifndef _RECONSTRUCT_H_
#define _RECONSTRUCT_H_

#include <vector>

#include <algebra3d/global.h>

/* Forward declarations */
/**
 * @namespace algebra3d
 * @brief Dependency: Code for algebra of 3D transformations.
 *
 * This namespace is defined by the
 * <a href="http://lcvmwww.epfl.ch/algebra3d">algebra3d</a> project.
 */
namespace algebra3d {
/**
 * @brief A class to store vectors of 3 elements.
 *
 * This class is defined by the
 * <a href="http://lcvmwww.epfl.ch/algebra3d">algebra3d</a> project.
 */
template<class fpType> class Vector3;
/**
 * @brief A class to store quaternions.
 *
 * This class is defined is defined by the
 * <a href="http://lcvmwww.epfl.ch/algebra3d">algebra3d</a> project.
 */
template<class fpType> class Quaternion;
} /* namespace algebra3d */

/**
 * @namespace cgdna_recon
 * @brief A namespace of the whole library.
 */
namespace cgdna_recon {

/**
 * Computes rigid base model parameters from the provided data and stores them
 * in the provided array.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #computeParamsD is another name for #computeParams \<double\> for
 *  Vector3 s and algebra3d#Quaternion s.
 *  * #computeParamsF is another name for #computeParams \<float\> for
 *  Vector3 s and Quaternion s.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size).
 * @param[in] r1 An array of vectors specifying positions of bases of the main
 * strand.
 * @param[in] r2 An array of vectors specifying positions of bases of the
 * complimentary strand.
 * @param[in] d1 An array of quaternions specifying rotations of bases of the
 * main strand.
 * @param[in] d2 An array of quaternions specifying rotations of bases of the
 * complimentary strand.
 * @param[out] params 1-D array to store the computed parameters (should be at
 * least of length 12 * (n - 1) + 6).
 *
 * @see computeParamsD
 * @see computeParamsF
 */
template<class fpType>
void computeParams(int nbp, algebra3d::Vector3<fpType>* r1,
		algebra3d::Vector3<fpType>* r2, algebra3d::Quaternion<fpType> *d1,
		algebra3d::Quaternion<fpType> *d2, fpType *params);

/**
 * A name for the #computeParams function to use with double (equal to
 * computeParams \<double\>) for Vector3 s and Quaternion s.
 *
 * @see computeParams
 */
void (* const computeParamsD)(int nbp, algebra3d::Vector3<double>* r1,
		algebra3d::Vector3<double>* r2, algebra3d::Quaternion<double> *d1,
		algebra3d::Quaternion<double> *d2,
		double *params) = computeParams<double>;

/**
 * A name for the #computeParams function to use with float (equal to
 * computeParams \<float\>) for Vector3 s and Quaternion s.
 *
 * @see computeParams
 */
void (* const computeParamsF)(int nbp, algebra3d::Vector3<float>* r1,
		algebra3d::Vector3<float>* r2, algebra3d::Quaternion<float> *d1,
		algebra3d::Quaternion<float> *d2, float *params) = computeParams<float>;

/**
 * Computes rigid base model parameters from the provided data and stores them
 * in the provided array.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #computeParamsArrD is another name for #computeParamsArr \<double\> for
 *  arrays of doubles.
 *  * #computeParamsArrF is another name for #computeParamsArr \<float\> for
 *  arrays of floats.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size).
 * @param[in] r1 Double pointer (assumed to be [n][3]) specifying positions of
 * bases of the main strand.
 * @param[in] r2 Double pointer (assumed to be [n][3]) specifying positions of
 * bases of the complimentary strand.
 * @param[in] d1 Triple pointer (assumed to be [n][3][3]) specifying rotations
 * of bases of the main strand. The [3][3] arrays are column-major matrices:
 * [0][2] is the third component of the first frame vector.
 * @param[in] d2 Triple pointer (assumed to be [n][3][3]) specifying rotations
 * of bases of the complimentary strand. The [3][3] arrays are column-major
 * matrices: [0][2] is the third component of the first frame vector.
 * @param[out] params 1-D array to store the computed parameters (should be at
 * least of length 12 * (n - 1) + 6).
 *
 * @see computeParamsArrD
 * @see computeParamsArrF
 */
template<class fpType>
void computeParamsArr(int nbp, fpType** r1, fpType** r2, fpType***d1,
		fpType***d2, fpType *params);

/**
 * A name for the #computeParamsArr function to use with double (equal to
 * computeParamsArr \<double\>) for arrays of doubles.
 *
 * @see computeParamsArr
 */
void (* const computeParamsArrD)(int nbp, double** r1, double** r2, double***d1,
		double ***d2, double *params) = computeParamsArr<double>;

/**
 * A name for the #computeParamsArr function to use with float (equal to
 * computeParamsArr \<float\>) for arrays of floats.
 *
 * @see computeParamsArr
 */
void (* const computeParamsArrF)(int nbp, float** r1, float** r2, float***d1,
		float***d2, float *params) = computeParamsArr<float>;

/**
 * Reconstructs base frames from rigid base model parameters.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructD is another name for #reconstruct \<double\> for Vector3 s
 *  and Quaternion s.
 *  * #reconstructF is another name for #reconstruct \<float\> for Vector3 s
 *  and Quaternion s.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be at least
 * of length 12 * (n - 1) + 6).
 * @param[out] r1 An array to store the computed positions of bases of the
 * main strand.
 * @param[out] r2 An array to store vectors specifying the computed positions
 * of bases of the complimentary strand.
 * @param[out] d1 An array to store quaternions specifying the computed
 * rotations of bases of the main strand.
 * @param[out] d2 An array to store quaternions specifying the computed
 * rotations of bases of the complimentary strand.
 *
 * @see reconstructD
 * @see reconstructF
 */
template<class fpType>
void reconstruct(int nbp, const fpType *params, algebra3d::Vector3<fpType>* r1,
		algebra3d::Vector3<fpType>* r2, algebra3d::Quaternion<fpType> *d1,
		algebra3d::Quaternion<fpType> *d2);

/**
 * A name for the #reconstruct function to use with double (equal to
 * reconstruct \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstruct
 */
void (* const reconstructD)(int nbp, const double *params,
		algebra3d::Vector3<double>* r1, algebra3d::Vector3<double>* r2,
		algebra3d::Quaternion<double> *d1,
		algebra3d::Quaternion<double> *d2) = reconstruct<double>;

/**
 * A name for the #reconstruct function to use with float (equal to
 * reconstruct \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstruct
 */
void (* const reconstructF)(int nbp, const float *params,
		algebra3d::Vector3<float>* r1, algebra3d::Vector3<float>* r2,
		algebra3d::Quaternion<float> *d1,
		algebra3d::Quaternion<float> *d2) = reconstruct<float>;

/**
 * Reconstructs base frames from rigid base model parameters.
 *
 * Uses Matrix3 for composing rotations and rotating vectors, Vector3 (Cayley
 * vectors) for evaluating half rotations.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructD is another name for #reconstruct \<double\> for Vector3 s
 *  and Matrix3 s.
 *  * #reconstructF is another name for #reconstruct \<float\> for Vector3 s
 *  and Matrix3 s.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be at least
 * of length 12 * (n - 1) + 6).
 * @param[out] r1 An array to store the computed positions of bases of the
 * main strand.
 * @param[out] r2 An array to store vectors specifying the computed positions
 * of bases of the complimentary strand.
 * @param[out] d1 An array to store matrices specifying the computed
 * rotations of bases of the main strand.
 * @param[out] d2 An array to store matrices specifying the computed
 * rotations of bases of the complimentary strand.
 *
 * @see reconstructMD
 * @see reconstructMF
 */
template<class fpType>
void reconstructM(int nbp, const fpType *params, algebra3d::Vector3<fpType>* r1,
		algebra3d::Vector3<fpType>* r2, algebra3d::Matrix3<fpType> *d1,
		algebra3d::Matrix3<fpType> *d2);

/**
 * A name for the #reconstructM function to use with double (equal to
 * reconstructM \<double\>) for Vector3 s and Matrix3 s.
 *
 * @see reconstruct
 */
void (* const reconstructMD)(int nbp, const double *params,
		algebra3d::Vector3<double>* r1, algebra3d::Vector3<double>* r2,
		algebra3d::Matrix3<double> *d1,
		algebra3d::Matrix3<double> *d2) = reconstructM<double>;

/**
 * A name for the #reconstructM function to use with float (equal to
 * reconstructM \<float\>) for Vector3 s and Matrix3 s.
 *
 * @see reconstruct
 */
void (* const reconstructMF)(int nbp, const float *params,
		algebra3d::Vector3<float>* r1, algebra3d::Vector3<float>* r2,
		algebra3d::Matrix3<float> *d1,
		algebra3d::Matrix3<float> *d2) = reconstructM<float>;

/**
 * Reconstructs base frames from rigid base model parameters. The first
 * std#vector should be of size 12*N -6 so that the number of base pairs (N) is
 * computed using that size.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructVecD is another name for #reconstructVec \<double\> for
 *  Vector3 s and Quaternion s.
 *  * #reconstructVecF is another name for #reconstructVec \<float\> for
 *  Vector3 s and Quaternion s.
 *
 * @param[in] params std#vector with rigid base parameters (should be at least
 * of length 12 * (n - 1) + 6).
 * @param[out] r1 std#vector to store the computed positions of bases of the
 * main strand.
 * @param[out] r2 std#vector to store vectors specifying the computed positions
 * of bases of the complimentary strand.
 * @param[out] d1 std#vector to store quaternions specifying the computed
 * rotations of bases of the main strand.
 * @param[out] d2 std#vector to store quaternions specifying the computed
 * rotations of bases of the complimentary strand.
 *
 * @see reconstructVecD
 * @see reconstructVecF
 */
template<class fpType>
void reconstructVec(const std::vector<fpType> &params,
		std::vector<algebra3d::Vector3<fpType> > &r1,
		std::vector<algebra3d::Vector3<fpType> > &r2,
		std::vector<algebra3d::Quaternion<fpType> > &d1,
		std::vector<algebra3d::Quaternion<fpType> > &d2);

/**
 * A name for the #reconstructVec function to use with double (equal to
 * reconstructVec \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructVec
 */
void (* const reconstructVecD)(const std::vector<double> &params,
		std::vector<algebra3d::Vector3<double> > &r1,
		std::vector<algebra3d::Vector3<double> > &r2,
		std::vector<algebra3d::Quaternion<double> > &d1,
		std::vector<algebra3d::Quaternion<double> > &d2) = reconstructVec<double>;

/**
 * A name for the #reconstructVec function to use with float (equal to
 * reconstructVec \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructVec
 */
void (* const reconstructVecF)(const std::vector<float> &params,
		std::vector<algebra3d::Vector3<float> > &r1,
		std::vector<algebra3d::Vector3<float> > &r2,
		std::vector<algebra3d::Quaternion<float> > &d1,
		std::vector<algebra3d::Quaternion<float> > &d2) = reconstructVec<float>;

/**
 * Reconstructs base frames from rigid base model parameters. The first
 * std#vector should be of size 12*N -6 so that the number of base pairs (N) is
 * computed using that size.
 *
 * Uses Matrix3 for composing rotations and rotating vectors, Vector3 (Cayley
 * vectors) for evaluating half rotations.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructMVecD is another name for #reconstructVec \<double\> for
 *  Vector3 s and Matrix3 s.
 *  * #reconstructMVecF is another name for #reconstructVec \<float\> for
 *  Vector3 s and Matrix3 s.
 *
 * @param[in] params std#vector with rigid base parameters (should be at least
 * of length 12 * (n - 1) + 6).
 * @param[out] r1 std#vector to store the computed positions of bases of the
 * main strand.
 * @param[out] r2 std#vector to store vectors specifying the computed positions
 * of bases of the complimentary strand.
 * @param[out] d1 std#vector to store matrices specifying the computed
 * rotations of bases of the main strand.
 * @param[out] d2 std#vector to store matrices specifying the computed
 * rotations of bases of the complimentary strand.
 *
 * @see reconstructMVecD
 * @see reconstructMVecF
 */
template<class fpType>
void reconstructMVec(const std::vector<fpType> &params,
		std::vector<algebra3d::Vector3<fpType> > &r1,
		std::vector<algebra3d::Vector3<fpType> > &r2,
		std::vector<algebra3d::Matrix3<fpType> > &d1,
		std::vector<algebra3d::Matrix3<fpType> > &d2);

/**
 * A name for the #reconstructMVec function to use with double (equal to
 * reconstructMVec \<double\>) for Vector3 s and Matrix3 s.
 *
 * @see reconstructMVec
 */
void (* const reconstructMVecD)(const std::vector<double> &params,
		std::vector<algebra3d::Vector3<double> > &r1,
		std::vector<algebra3d::Vector3<double> > &r2,
		std::vector<algebra3d::Matrix3<double> > &d1,
		std::vector<algebra3d::Matrix3<double> > &d2) = reconstructMVec<double>;

/**
 * A name for the #reconstructMVec function to use with float (equal to
 * reconstructMVec \<float\>) for Vector3 s and Matrix3 s.
 *
 * @see reconstructMVec
 */
void (* const reconstructMVecF)(const std::vector<float> &params,
		std::vector<algebra3d::Vector3<float> > &r1,
		std::vector<algebra3d::Vector3<float> > &r2,
		std::vector<algebra3d::Matrix3<float> > &d1,
		std::vector<algebra3d::Matrix3<float> > &d2) = reconstructMVec<float>;

/**
 * Reconstructs base frames from rigid base model parameters.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors and for return values.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructArrD is another name for #reconstructArr \<double\> for
 *  arrays of doubles.
 *  * #reconstructArrF is another name for #reconstructArr \<float\> for arrays
 *  of floats.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be at least
 * of length 12 * (n - 1) + 6).
 * @param[out] r1 Double pointer (assumed to be [n][3]) to store positions of
 * bases of the main strand.
 * @param[out] r2 Double pointer (assumed to be [n][3]) to store positions of
 * bases of the complimentary strand.
 * @param[out] d1 Triple pointer (assumed to be [n][3][3]) to store rotations
 * of bases of the main strand. The [3][3] arrays are column-major matrices:
 * [0][2] is the third component of the first frame vector.
 * @param[out] d2 Triple pointer (assumed to be [n][3][3]) to store rotations
 * of bases of the complimentary strand. The [3][3] arrays are column-major
 * matrices: [0][2] is the third component of the first frame vector.
 *
 * @see reconstructArrD
 * @see reconstructArrF
 */
template<class fpType>
void reconstructArr(int nbp, fpType *params, fpType** r1, fpType** r2,
		fpType***d1, fpType***d2);

/**
 * A name for the #reconstructArr function to use with double (equal to
 * reconstructArr \<double\>) for arrays of doubles.
 *
 * @see reconstructArr
 */
void (* const reconstructArrD)(int nbp, double *params, double** r1,
		double** r2, double ***d1, double ***d2) = reconstructArr<double>;

/**
 * A name for the #reconstructArr function to use with float (equal to
 * reconstructArr \<float\>) for arrays of doubles.
 *
 * @see reconstructArr
 */
void (* const reconstructArrF)(int nbp, float *params, float** r1, float** r2,
		float ***d1, float ***d2) = reconstructArr<float>;

/**
 * Reconstructs base pair frames from rigid base or rigid base pair model
 * parameters.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructBPD is another name for #reconstructBP \<double\> for
 *  Vector3 s and Quaternion s.
 *  * #reconstructBPF is another name for #reconstructBP \<float\> for
 *  Vector3 s and Quaternion s.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be of length
 * 6 * n).
 * @param[out] r An array to store the computed positions of base pairs.
 * @param[out] d An array to store quaternions specifying the computed
 * rotations of base pairs.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructBPD
 * @see reconstructBPF
 */
template<class fpType>
void reconstructBP(int nbp, const fpType *params, algebra3d::Vector3<fpType>* r,
		algebra3d::Quaternion<fpType> *d, bool baseModel);

/**
 * A name for the #reconstructBP function to use with double (equal to
 * reconstructBP \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBP
 */
void (* const reconstructBPD)(int nbp, const double *params,
		algebra3d::Vector3<double>* r, algebra3d::Quaternion<double> *d,
		bool baseModel) = reconstructBP<double>;

/**
 * A name for the #reconstructBP function to use with float (equal to
 * reconstructBP \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBP
 */
void (* const reconstructBPF)(int nbp, const float *params,
		algebra3d::Vector3<float>* r, algebra3d::Quaternion<float> *d,
		bool baseModel) = reconstructBP<float>;

/**
 * Reconstructs base pair frames from rigid base or rigid base pair model
 * parameters.
 *
* Uses Matrix3 for composing rotations and rotating vectors, Vector3 (Cayley
* vectors) for evaluating half rotations.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructBPMD is another name for #reconstructBPM \<double\> for
 *  Vector3 s and Matrix3 s.
 *  * #reconstructBPMF is another name for #reconstructBPM \<float\> for
 *  Vector3 s and Matrix3 s.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be of length
 * 6 * n).
 * @param[out] r An array to store the computed positions of base pairs.
 * @param[out] d An array to store matrices specifying the computed
 * rotations of base pairs.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructBPMD
 * @see reconstructBPMF
 */
template<class fpType>
void reconstructBPM(int nbp, const fpType *params, algebra3d::Vector3<fpType>* r,
		algebra3d::Matrix3<fpType> *d, bool baseModel);

/**
 * A name for the #reconstructBPM function to use with double (equal to
 * reconstructBPM \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBPM
 */
void (* const reconstructBPMD)(int nbp, const double *params,
		algebra3d::Vector3<double>* r, algebra3d::Quaternion<double> *d,
		bool baseModel) = reconstructBP<double>;

/**
 * A name for the #reconstructBPM function to use with float (equal to
 * reconstructBPM \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBPM
 */
void (* const reconstructBPMF)(int nbp, const float *params,
		algebra3d::Vector3<float>* r, algebra3d::Quaternion<float> *d,
		bool baseModel) = reconstructBP<float>;

/**
 * Reconstructs base pair frames from rigid base or rigid base pair model
 * parameters. The first std#vector should be of size 6*N so that the number of
 * base pairs (N) is computed using that size.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructBPVecD is another name for #reconstructBPVec \<double\> for
 *  Vector3 s and Quaternion s.
 *  * #reconstructBPVecF is another name for #reconstructBPVec \<float\> for
 *  Vector3 s and Quaternion s.
 *
 * @param[in] params std#vector with rigid base pair parameters (should be of
 * length 6 * N.
 * @param[out] r std#vector to store the computed positions of base pairs.
 * @param[out] d std#vector to store quaternions specifying the computed
 * rotations of base pairs.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructBPVecD
 * @see reconstructBPVecF
 */
template<class fpType>
void reconstructBPVec(const std::vector<fpType> &params,
		std::vector<algebra3d::Vector3<fpType> > &r,
		std::vector<algebra3d::Quaternion<fpType> > &d, bool baseModel);

/**
 * A name for the #reconstructBPVec function to use with double (equal to
 * reconstructBPVec \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBPVec
 */
void (* const reconstructBPVecD)(const std::vector<double> &params,
		std::vector<algebra3d::Vector3<double> > &r,
		std::vector<algebra3d::Quaternion<double> > &d,
		bool baseModel) = reconstructBPVec<double>;

/**
 * A name for the #reconstructBPVec function to use with float (equal to
 * reconstructBPVec \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBPVec
 */
void (* const reconstructBPVecF)(const std::vector<float> &params,
		std::vector<algebra3d::Vector3<float> > &r,
		std::vector<algebra3d::Quaternion<float> > &d,
		bool baseModel) = reconstructBPVec<float>;

/**
 * Reconstructs base pair frames from rigid base or rigid base pair model
 * parameters. The first std#vector should be of size 6*N so that the number of
 * base pairs (N) is computed using that size.
 *
 * Uses Matrix3 for composing rotations and rotating vectors, Vector3 (Cayley
 * vectors) for evaluating half rotations.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructBPMVecD is another name for #reconstructBPVec \<double\> for
 *  Vector3 s and Quaternion s.
 *  * #reconstructBPMVecF is another name for #reconstructBPVec \<float\> for
 *  Vector3 s and Quaternion s.
 *
 * @param[in] params std#vector with rigid base pair parameters (should be of
 * length 6 * N.
 * @param[out] r std#vector to store the computed positions of base pairs.
 * @param[out] d std#vector to store quaternions specifying the computed
 * rotations of base pairs.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructBPMVecD
 * @see reconstructBPMVecF
 */
template<class fpType>
void reconstructBPMVec(const std::vector<fpType> &params,
		std::vector<algebra3d::Vector3<fpType> > &r,
		std::vector<algebra3d::Matrix3<fpType> > &d, bool baseModel);

/**
 * A name for the #reconstructBPMVec function to use with double (equal to
 * reconstructBPMVec \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBPMVec
 */
void (* const reconstructBPMVecD)(const std::vector<double> &params,
		std::vector<algebra3d::Vector3<double> > &r,
		std::vector<algebra3d::Matrix3<double> > &d,
		bool baseModel) = reconstructBPMVec<double>;

/**
 * A name for the #reconstructBPMVec function to use with float (equal to
 * reconstructBPMVec \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructBPMVec
 */
void (* const reconstructBPMVecF)(const std::vector<float> &params,
		std::vector<algebra3d::Vector3<float> > &r,
		std::vector<algebra3d::Matrix3<float> > &d,
		bool baseModel) = reconstructBPMVec<float>;

/**
 * Reconstructs base frames from rigid base or rigid base pair model
 * parameters.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors and return values.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructArrD is another name for #reconstructArr \<double\> for
 *  arrays of doubles.
 *  * #reconstructArrF is another name for #reconstructArr \<float\> for arrays
 *  of floats.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be at least
 * of length 6 * n.
 * @param[out] r Double pointer (assumed to be [n][3]) to store positions of
 * base pairs.
 * @param[out] d Triple pointer (assumed to be [n][3][3]) to store rotations
 * of base pairs. The [3][3] arrays store matrices column-major:
 * [0][2] is the third component of the first frame vector.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructBPArrD
 * @see reconstructBPArrF
 */
template<class fpType>
void reconstructBPArr(int nbp, fpType *params, fpType** r, fpType***d,
		bool baseModel);

/**
 * A name for the #reconstructBPArr function to use with double (equal to
 * reconstructBPArr \<double\>) for arrays of doubles.
 *
 * @see reconstructBPArr
 */
void (* const reconstructBPArrD)(int nbp, double *params, double** r,
		double ***d, bool baseModel) = reconstructBPArr<double>;

/**
 * A name for the #reconstructBPArr function to use with float (equal to
 * reconstructBPArr \<float\>) for arrays of doubles.
 *
 * @see reconstructBPArr
 */
void (* const reconstructBPArrF)(int nbp, float *params, float** r, float ***d,
		bool baseModel) = reconstructBPArr<float>;

/**
 * Reconstructs the first, middle and last base frames from rigid base or rigid
 * base pair model parameters.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors.
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructFMLD is another name for #reconstructFML \<double\> for
 *  Vector3 s and Quaternion s.
 *  * #reconstructFMLF is another name for #reconstructFML \<float\> for
 *  Vector3 s and Quaternion s.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size)
 * @param[in] params 1-D array with rigid base parameters (should be at least
 * of length 12 * n - 6 for rigid bases or 6 * n for rigid base pairs).
 * @param[out] rFirst The computed position of the first base pair of the main
 * strand.
 * @param[out] dFirst The computed rotation of the first base pair of the main
 * strand.
 * @param[out] rMid The computed position of the middle base pair of the main
 * strand.
 * @param dMid The computed rotation of the middle base pair of the main
 * strand.
 * @param[out] rLast The computed position of the last base pair of the main
 * strand.
 * @param[out] dLast The computed rotation of the last base pair of the main
 * strand.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructFMLD
 * @see reconstructFMLF
 */
template<class fpType>
void reconstructFML(int nbp, fpType *params, //
		algebra3d::Vector3<fpType> &rFirst, //
		algebra3d::Quaternion<fpType> &dFirst, //
		algebra3d::Vector3<fpType> &rMid, //
		algebra3d::Quaternion<fpType> &dMid, //
		algebra3d::Vector3<fpType> &rLast, //
		algebra3d::Quaternion<fpType> &dLast, //
		bool baseModel);

/**
 * A name for the #reconstructFML function to use with double (equal to
 * reconstructFML \<double\>) for Vector3 s and Quaternion s.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors. Returns the Matrix3 s
 *
 * @see reconstructFML
 */
void (* const reconstructFMLD)(int nbp,
		double *params, //
		algebra3d::Vector3<double> &rFirst,
		algebra3d::Quaternion<double> &d1irst, //
		algebra3d::Vector3<double> &rMid, algebra3d::Quaternion<double> &dMid, //
		algebra3d::Vector3<double> &rLast, algebra3d::Quaternion<double> &dLast, //
		bool baseModel) = reconstructFML<double>;

/**
 * A name for the #reconstructFML function to use with float (equal to
 * reconstructFML \<float\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructFML
 */
void (* const reconstructFMLF)(int nbp, float *params, //
		algebra3d::Vector3<float> &rFirst, algebra3d::Quaternion<float> &d1irst, //
		algebra3d::Vector3<float> &rMid, algebra3d::Quaternion<float> &dMid, //
		algebra3d::Vector3<float> &rLast, algebra3d::Quaternion<float> &dLast, //
		bool baseModel)= reconstructFML<float>;

/**
 * Reconstructs the first, middle and last base frames from rigid base model
 * parameters.
 *
 * Uses Quaternion s for composing rotations and evaluating half rotations.
 * Converts Quaternion s to Matrix3 s to rotate vectors. Returns the Matrix3 s
 *
 * Two function pointers have been added to make it easier to use.
 *  * #reconstructFMLArrD is another name for #reconstructFMLArr \<double\> for
 *  arrays of doubles.
 *  * #reconstructFMLArrF is another name for #reconstructFMLArr \<float\> for
 *  arrays of floats.
 *
 * @param[in] nbp Number of base pairs (all arrays are assumed to have
 * appropriate size).
 * @param[in] params 1-D array with rigid base parameters (should be at least
 * of length 12 * (n - 1) + 6).
 * @param[out] rFirst Pointer (assumed to be [3]) to store position of the
 * first base pair of the main strand.
 * @param[out] dFirst Double pointer (assumed to be [3][3]) to store rotation
 * of the first base pair of the main strand. The [3][3] arrays are
 * column-major matrices: [0][2] is the third component of the first frame
 * vector.
 * @param[out] rMid Pointer (assumed to be [3]) to store position of the
 * middle base pair of the main strand.
 * @param[out] dMid Double pointer (assumed to be [3][3]) to store rotation
 * of the middle base pair of the main strand. The [3][3] arrays are
 * column-major matrices: [0][2] is the third component of the first frame
 * vector.
 * @param[out] rLast Pointer (assumed to be [3]) to store position of the
 * last base pair of the main strand.
 * @param[out] dLast Double pointer (assumed to be [3][3]) to store rotation
 * of the last base pair of the main strand. The [3][3] arrays are
 * column-major matrices: [0][2] is the third component of the first frame
 * vector.
 * @param[in] baseModel A flag indicating which model are the parameters for.
 * `true` - base model (12*n - 6 parameters), `false` -
 * base pair model (6 * n parameters)
 *
 * @see reconstructFMLArrD
 * @see reconstructFMLArrF
 */
template<class fpType>
void reconstructFMLArr(int nbp, fpType *params, //
		fpType *rFirst, //
		fpType **dFirst, //
		fpType *rMid, //
		fpType **dMid, //
		fpType *rLast, //
		fpType **dLast, //
		bool baseModel);

/**
 * A name for the #reconstructFMLArr function to use with double (equal to
 * reconstructFMLArr \<double\>) for Vector3 s and Quaternion s.
 *
 * @see reconstructFMLArr
 */
void (* const reconstructFMLArrD)(int nbp, double *params, //
		double *rFirst, //
		double **dFirst, //
		double *rMid, //
		double **dMid, //
		double *rLast, //
		double **dLast, //
		bool baseModel) = reconstructFMLArr<double>;

/**
 * A name for the #reconstructFMLArr function to use with float (equal to
 * reconstructFMLArr \<float\*) for Vector3 s and Quaternion s.
 *
 * @see reconstructFMLArr
 */
void (* const reconstructFMLArrF)(int nbp, float *params, //
		float *rFirst, //
		float **dFirst, //
		float *rMid, //
		float **dMid, //
		float *rLast, //
		float **dLast, //
		bool baseModel) = reconstructFMLArr<float>;

} // namespace

#endif // _RECONSTRUCT_H_

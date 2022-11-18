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
 * Implementation of the Matrix class.
 *
 * @author Jarek Glowacki
 */

#include<algebra3d/global.h>
#include<algebra3d/Matrix3.h>
#include<algebra3d/Vector3.h>
#include<algebra3d/Quaternion.h>
#include<cstring>

#include<iostream>

namespace algebra3d {

template<class fpType>
const int Matrix3<fpType>::ADJUGATE_INDICES[3][3][4][2] = { //
		{ //
				{ { 1, 1 }, { 2, 2 }, { 2, 1 }, { 1, 2 } }, //
						{ { 2, 1 }, { 0, 2 }, { 0, 1 }, { 2, 2 } }, //
						{ { 0, 1 }, { 1, 2 }, { 1, 1 }, { 0, 2 } } //
				}, //
				{ //
				{ { 2, 0 }, { 1, 2 }, { 1, 0 }, { 2, 2 } }, //
						{ { 0, 0 }, { 2, 2 }, { 2, 0 }, { 0, 2 } }, //
						{ { 1, 0 }, { 0, 2 }, { 0, 0 }, { 1, 2 } } //
				}, //
				{ //
				{ { 1, 0 }, { 2, 1 }, { 2, 0 }, { 1, 1 } }, //
						{ { 2, 0 }, { 0, 1 }, { 0, 0 }, { 2, 1 } }, //
						{ { 0, 0 }, { 1, 1 }, { 1, 0 }, { 0, 1 } } //
				} //
		};

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::initialize(bool identity) {
	// Initialize as 0
	memset(data, 0, 9 * sizeof(fpType));
	// For identity set diagonal elements to 1.0
	for (int row = 0; row < 3; ++row) {
		data[row][row] = identity;
	}
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(bool identity) {
	initialize(identity);
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(const fpType (&data)[3][3]) {
	memcpy(this->data, data, 9 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(const fpType*data) {
	memcpy(this->data, data, 9 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(const Vector3<fpType> &col1,
		const Vector3<fpType> &col2, const Vector3<fpType> &col3) {
	col1.get(data[0]);
	col2.get(data[1]);
	col3.get(data[2]);
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(fpType m11, fpType m12, fpType m13, //
		fpType m21, fpType m22, fpType m23, //
		fpType m31, fpType m32, fpType m33) {
	data[0][0] = m11;
	data[0][1] = m21;
	data[0][2] = m31;

	data[1][0] = m12;
	data[1][1] = m22;
	data[1][2] = m32;

	data[2][0] = m31;
	data[2][1] = m32;
	data[2][2] = m33;
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(const Matrix3<fpType> &mat) {
	memcpy(this->data, mat.data, 9 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(const Vector3<fpType> &cvec) {
	*this = Quaternion<fpType>(cvec);
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(Vector3<fpType> axis, fpType angle) {
	initialize(true);

	axis.normalize();
	fpType sina = sin(angle);
	fpType cosa = cos(angle);

	// Rodrigues' rotation formula
	*this *= cosa;
	*this += ((fpType) 1.0 - cosa) * axis.outer(axis)
			+ sina * axis.crossMatrix();
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::Matrix3(const Quaternion<fpType> &quat) {
	fpType q1q1 = quat.x * quat.x;
	fpType q1q2 = quat.x * quat.y;
	fpType q1q3 = quat.x * quat.z;
	fpType q1q4 = quat.x * quat.w;

	fpType q2q2 = quat.y * quat.y;
	fpType q2q3 = quat.y * quat.z;
	fpType q2q4 = quat.y * quat.w;

	fpType q3q3 = quat.z * quat.z;
	fpType q3q4 = quat.z * quat.w;

	fpType q4q4 = quat.w * quat.w;

	// d1 (first column)
	data[0][0] = q1q1 - q2q2 - q3q3 + q4q4;
	data[0][1] = 2.0 * (q1q2 + q3q4);
	data[0][2] = 2.0 * (q1q3 - q2q4);

	// d2 (second column)
	data[1][0] = 2.0 * (q1q2 - q3q4);
	data[1][1] = -q1q1 + q2q2 - q3q3 + q4q4;
	data[1][2] = 2.0 * (q2q3 + q1q4);

	// d3 (first column)
	data[2][0] = 2.0 * (q1q3 + q2q4);
	data[2][1] = 2.0 * (q2q3 - q1q4);
	data[2][2] = -q1q1 - q2q2 + q3q3 + q4q4;
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType>::~Matrix3() {
	// Does nothing
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::set(const fpType data[3][3]) {
	memcpy(this->data, data, 9 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::set(const Vector3<fpType> &col1,
		const Vector3<fpType> &col2, const Vector3<fpType> &col3) {
	col1.get(data[0]);
	col2.get(data[1]);
	col3.get(data[2]);
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::get(fpType storage[3][3]) const {
	memcpy(storage, this->data, 9 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::get(fpType * const col1, fpType * const col2,
		fpType * const col3) const {
	memcpy(col1, data[0], 3 * sizeof(fpType));
	memcpy(col2, data[1], 3 * sizeof(fpType));
	memcpy(col3, data[2], 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::setColumn(int index, const Vector3<fpType> &col) {
	index = (index < 0) ? 0 : (index > 2) ? 2 : index;
	col.get(data[index]);
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::getColumn(int index, fpType * const &col) {
	index = (index < 0) ? 0 : (index > 2) ? 2 : index;
	memcpy(col, data[index], 3);
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Matrix3<fpType>::getColumn(int index) {
	index = (index < 0) ? 0 : (index > 2) ? 2 : index;
	return Vector3<fpType>(data[index]);
}

/****************************************************************************/
template<class fpType>
fpType Matrix3<fpType>::getTrace() const {
	return data[0][0] + data[1][1] + data[2][2];
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::transpose() {
	fpType tmp;
	for (int row = 0; row < 3; ++row) {
		for (int col = row + 1; col < 3; ++col) {
			tmp = data[col][row];
			data[col][row] = data[row][col];
			data[row][col] = tmp;
		}
	}
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType> Matrix3<fpType>::getTranspose() const {
	Matrix3<fpType> result(data);
	result.transpose();
	return result;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Matrix3<fpType>::getAxis() const {
	// Don't use inverse Rodrigues' formula not to fail for rotations by pi
	return Quaternion<fpType>(*this).getAxis();
}

/****************************************************************************/
template<class fpType>
fpType Matrix3<fpType>::getAngle() const {
	// The formula using trace is not used for accuracy reasons at pi
	return Quaternion<fpType>(*this).getAngle();
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::invert() {
	fpType de = det();

	if (algebra3d::compare<fpType>(de, 0.0) > 0) {
		fpType detInv = 1.0 / de;
		Vector3<fpType> x0(data[0]);
		Vector3<fpType> x1(data[1]);
		Vector3<fpType> x2(data[2]);

		Vector3<fpType> cross0 = detInv * x1.cross(x2);
		Vector3<fpType> cross1 = detInv * x2.cross(x0);
		Vector3<fpType> cross2 = detInv * x0.cross(x1);

		for (int i = 0; i < 3; ++i) {
			data[i][0] = cross0[i];
			data[i][1] = cross1[i];
			data[i][2] = cross2[i];
		}
//		Matrix3 tmp(data);
//		const int (*indeces)[2];
//		for (int row = 0; row < 3; ++row) {
//			for (int col = 0; col < 3; ++col) {
//				indeces = ADJUGATE_INDICES[col][row];
//				data[col][row] = detInv //
//				* ((tmp[indeces[0][0]][indeces[0][1]] //
//				* tmp[indeces[1][0]][indeces[1][1]]) //
//				- (tmp[indeces[2][0]][indeces[2][1]] //
//				* tmp[indeces[3][0]][indeces[3][1]]));
//			}
//		}
	} else {
		memset(data, 0, 9 * sizeof(fpType));
	}
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType> Matrix3<fpType>::getInverse() const {
	Matrix3<fpType> result(data);
	result.invert();
	return result;
}

/****************************************************************************/
template<class fpType>
fpType Matrix3<fpType>::det() {
	return data[0][0] * data[1][1] * data[2][2] //
	+ data[1][0] * data[2][1] * data[0][2] //
	+ data[2][0] * data[0][1] * data[1][2] //
	- data[0][0] * data[2][1] * data[1][2] //
	- data[1][0] * data[0][1] * data[2][2] //
	- data[2][0] * data[1][1] * data[0][2];
}

/****************************************************************************/
template<class fpType>
bool Matrix3<fpType>::isOrthogonal() {
	return Matrix3(true) == (*this) * getTranspose();
}

/****************************************************************************/
template<class fpType>
void Matrix3<fpType>::operator=(const Matrix3<fpType> &mat) {
	memcpy(data, mat.data, 9 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
const fpType (&Matrix3<fpType>::operator[](int i1) const)[3] {
			i1 = (i1 < 0) ? 0 : ((i1 > 2) ? 2 : i1);
			return data[i1];
		}

		/****************************************************************************/
		template<class fpType>
		fpType (&Matrix3<fpType>::operator[](int i1))[3] {
					i1 = (i1 < 0) ? 0 : ((i1 > 2) ? 2 : i1);
					return data[i1];
				}

				/****************************************************************************/
				template<class fpType>
				bool Matrix3<fpType>::compare(const Matrix3<fpType> &mat,
						const fpType &precision) const {
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							if (algebra3d::compare<fpType>(data[col][row],
									mat.data[col][row], precision) != 0) {
								return false;
							}
						}
					}
					return true;
				}

				/****************************************************************************/
				template<class fpType>
				bool Matrix3<fpType>::operator==(
						const Matrix3<fpType> &mat) const {
					return compare(mat);
				}

				/****************************************************************************/
				template<class fpType>
				bool Matrix3<fpType>::operator!=(
						const Matrix3<fpType> &mat) const {
					return !((*this) == mat);
				}

				/****************************************************************************/
				template<class fpType>
				void Matrix3<fpType>::operator+=(const Matrix3<fpType> &mat) {
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							data[col][row] += mat.data[col][row];
						}
					}
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> Matrix3<fpType>::operator+(
						const Matrix3<fpType> &mat) const {
					Matrix3<fpType> result(data);
					result += mat;
					return result;
				}

				/****************************************************************************/
				template<class fpType>
				void Matrix3<fpType>::operator-=(const Matrix3<fpType> &mat) {
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							data[col][row] -= mat.data[col][row];
						}
					}
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> Matrix3<fpType>::operator-() const {
					Matrix3 result(data);
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							result[col][row] = -result[col][row];
						}
					}
					return result;
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> Matrix3<fpType>::operator-(
						const Matrix3<fpType> &mat) const {
					Matrix3<fpType> result(data);
					result -= mat;
					return result;
				}

				/****************************************************************************/
				template<class fpType>
				Vector3<fpType> Matrix3<fpType>::operator*(
						const Vector3<fpType> &vec) const {
//					return Vector3<fpType>(data[0]) * vec.x
//							+ Vector3<fpType>(data[1]) * vec.y
//							+ Vector3<fpType>(data[2]) * vec.z;

					return Vector3<fpType>(
							data[0][0] * vec.x + data[1][0] * vec.y
									+ data[2][0] * vec.z,
							data[0][1] * vec.x + data[1][1] * vec.y
									+ data[2][1] * vec.z,
							data[0][2] * vec.x + data[1][2] * vec.y
									+ data[2][2] * vec.z);
				}

				/****************************************************************************/
				template<class fpType>
				void Matrix3<fpType>::operator*=(const Matrix3<fpType> &mat) {
					Matrix3<fpType> tmpm(data);
					fpType tmp;
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							tmp = 0.0;
							for (int j = 0; j < 3; ++j) {
								tmp += tmpm[j][row] * mat.data[col][j];
							}
							data[col][row] = tmp;
						}
					}
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> Matrix3<fpType>::operator*(
						const Matrix3<fpType> &mat) const {
					Matrix3<fpType> result(data);
					result *= mat;
					return result;
				}

				/****************************************************************************/
				template<class fpType>
				void Matrix3<fpType>::operator*=(const fpType &scale) {
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							data[col][row] *= scale;
						}
					}
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> Matrix3<fpType>::operator*(
						const fpType &scale) const {
					Matrix3 result(data);
					result *= scale;
					return result;
				}

				/****************************************************************************/
				template<class fpType>
				void Matrix3<fpType>::operator/=(const fpType &scale) {
					for (int row = 0; row < 3; ++row) {
						for (int col = 0; col < 3; ++col) {
							data[col][row] /= scale;
						}
					}
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> Matrix3<fpType>::operator/(
						const fpType &scale) const {
					Matrix3 result(data);
					result /= scale;
					return result;
				}

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> operator*(const double &scale,
						const Matrix3<fpType> &mat) {
					return mat * scale;
				}

				// Define the specializations that will be usable
				template
				Matrix3<double> operator*(const double &scale,
						const Matrix3<double> &mat);
				template
				Matrix3<float> operator*(const double &scale,
						const Matrix3<float> &mat);

				/****************************************************************************/
				template<class fpType>
				Matrix3<fpType> operator*(const float &scale,
						const Matrix3<fpType> &mat) {
					return mat * scale;
				}

				// Define the specializations that will be usable
				template
				Matrix3<double> operator*(const float &scale,
						const Matrix3<double> &mat);
				template
				Matrix3<float> operator*(const float &scale,
						const Matrix3<float> &mat);

				/****************************************************************************/
				template<class fpType>
				std::ostream& operator<<(std::ostream &stream,
						const Matrix3<fpType> &mat) {
					for (int row = 0; row < 2; ++row) {
						for (int col = 0; col < 3; ++col) {
							stream << mat[col][row] << " ";
						}
						stream << std::endl;
					}
					for (int col = 0; col < 3; ++col) {
						stream << mat[col][2] << " ";
					}
					return stream;
				}

				// Define the specializations that will be usable
				template
				std::ostream& operator<<<double>(std::ostream &stream,
						const Matrix3<double> &mat);
				template
				std::ostream& operator<<<float>(std::ostream &stream,
						const Matrix3<float> &mat);

				/****************************************************************************/
				template<class fpType>
				std::istream& operator>>(std::istream &stream,
						Matrix3<fpType> &mat) {
					for (int i = 0; i < 3; ++i) {
						for (int j = 0; j < 3; ++j) {
							stream >> mat.data[j][i];
						}
					}
					return stream;
				}

				// Define the specializations that will be usable
				template
				std::istream& operator>><double>(std::istream &stream,
						Matrix3<double> &mat);
				template
				std::istream& operator>><float>(std::istream &stream,
						Matrix3<float> &mat);

				/****************************************************************************/
				// Define the specializations that will be usable
				template class Matrix3<double> ;
				template class Matrix3<float> ;

				} /* namespace algebra3d */

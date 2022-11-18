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
 * Implementation of the Vector class.
 *
 * @author Jarek Glowacki
 */

#include<algebra3d/global.h>
#include<algebra3d/Vector3.h>
#include<algebra3d/Matrix3.h>
#include<algebra3d/Quaternion.h>
#include<cstring>

namespace algebra3d {

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3() :
		x(data[0]), y(data[1]), z(data[2]) {
	memset(data, 0, 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3(const fpType * const data) :
		x(this->data[0]), y(this->data[1]), z(this->data[2]) {
	memcpy(this->data, data, 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3(fpType x, fpType y, fpType z) :
		x(data[0]), y(data[1]), z(data[2]) {
	this->x = x;
	this->y = y;
	this->z = z;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3(const Vector3<fpType> &vec) :
		x(data[0]), y(data[1]), z(data[2]) {
	memcpy(this->data, vec.data, 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3(const Vector3<fpType> &axis, fpType angle) :
		x(data[0]), y(data[1]), z(data[2]) {
	*this = axis;
	normalize();
	*this *= 2.0 * tan(angle * 0.5);
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3(const Quaternion<fpType> &quat) :
		x(data[0]), y(data[1]), z(data[2]) {
	x = 2.0 * quat.x / quat.w;
	y = 2.0 * quat.y / quat.w;
	z = 2.0 * quat.z / quat.w;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::Vector3(const Matrix3<fpType> &mat) :
		x(data[0]), y(data[1]), z(data[2]) {
	fpType mult = mat.getTrace();
	// Watch out for the singular case: mult == -1.0
	mult = 2.0 / (mult + 1);
	Matrix3<fpType> crossMat = mat - mat.getTranspose();
	x = crossMat[1][2] * mult;
	y = crossMat[2][0] * mult;
	z = crossMat[0][1] * mult;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType>::~Vector3() {
	// Does nothing
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::set(const fpType data[3]) {
	memcpy(this->data, data, 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::get(fpType *storage) const {
	memcpy(storage, this->data, 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
fpType Vector3<fpType>::dot(const Vector3<fpType> &vec) const {
	return data[0] * vec.data[0] //
	+ data[1] * vec.data[1] //
	+ data[2] * vec.data[2];
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::cross(const Vector3<fpType> &vec) const {
	return Vector3<fpType>(y * vec.z - z * vec.y, //
	z * vec.x - x * vec.z, //
	x * vec.y - y * vec.x);
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType> Vector3<fpType>::outer(const Vector3<fpType> &vec) const {
	return Matrix3<fpType>(*this * vec.x, *this * vec.y, *this * vec.z);
}

/****************************************************************************/
template<class fpType>
Matrix3<fpType> Vector3<fpType>::crossMatrix() const {
	Matrix3<fpType> result;
	result[1][0] = -z;
	result[0][1] = z;

	result[2][0] = y;
	result[0][2] = -y;

	result[2][1] = -x;
	result[1][2] = x;
	return result;
}

/****************************************************************************/
template<class fpType>
fpType Vector3<fpType>::getNorm() const {
	return sqrt(this->dot(*this));
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::normalize() {
	fpType norm = getNorm();
	*this /= norm;
}

/****************************************************************************/
template<class fpType>
bool Vector3<fpType>::isZero() const {
	return (algebra3d::compare<fpType>(x, 0.0) == 0)
			&& (algebra3d::compare<fpType>(y, 0.0) == 0)
			&& (algebra3d::compare<fpType>(z, 0.0) == 0);
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::getUnit() const {
	fpType norm = getNorm();
	return *this / norm;
}

/****************************************************************************/
template<class fpType>
fpType Vector3<fpType>::getAngle() const {
	fpType length = getNorm();
	return 2.0 * atan(0.5 * length);
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::getSqrtRotation() const {
	double norm2 = this->dot(*this);
	// The formulae
	// cos(a/2) = 2cos^2(a/4) - 1
	// sin(a/2) = 2sin(a/4)cos(a/4)
	// can be inverted with uniqueness preserved as a is a the rotation angle so
	// cos(a/4) = sqrt((cos(a/2) + 1)/2)
	// sin(a/4) = sin(a/2)/(2cos(a/4))
	// wich implies
	// tan(a/4) = sin(a/2)/(cos^2(a/2) + 1)
	// using
	// cos(a/2) = 1/sqrt(tan^2(a/2) + 1)
	// sin(a/2) = tan(a/2)/sqrt(tan^2(a/2) + 1)
	// this implies the relation
	return *this * 2.0 / (2.0 + sqrt(norm2 + 4));
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::operator=(const Vector3<fpType> &vec) {
	memcpy(data, vec.data, 3 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
const fpType &Vector3<fpType>::operator[](int ind) const {
	ind = (ind < 0) ? 0 : ((ind > 2) ? 2 : ind);
	return data[ind];
}

/****************************************************************************/
template<class fpType>
fpType &Vector3<fpType>::operator[](int ind) {
	ind = (ind < 0) ? 0 : ((ind > 2) ? 2 : ind);
	return data[ind];
}

/****************************************************************************/
template<class fpType>
bool Vector3<fpType>::compare(const Vector3<fpType> &vec,
		const fpType &precision) const {
	return (algebra3d::compare(x, vec.x, precision) == 0.0)
			&& (algebra3d::compare(y, vec.y, precision) == 0.0)
			&& (algebra3d::compare(z, vec.z, precision) == 0.0);
}

/****************************************************************************/
template<class fpType>
bool Vector3<fpType>::operator==(const Vector3<fpType> &vec) const {
	return compare(vec);
}

/****************************************************************************/
template<class fpType>
bool Vector3<fpType>::operator!=(const Vector3<fpType> &vec) const {
	return !((*this) == vec);
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::operator+=(const Vector3<fpType> &vec) {
	x += vec.x;
	y += vec.y;
	z += vec.z;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator+(const Vector3<fpType> &vec) const {
	Vector3<fpType> result(data);
	result += vec;
	return result;
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::operator-=(const Vector3<fpType> &vec) {
	x -= vec.x;
	y -= vec.y;
	z -= vec.z;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator-() const {
	return Vector3<fpType>(-x, -y, -z);
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator-(const Vector3<fpType> &m) const {
	Vector3<fpType> result(data);
	result -= m;
	return result;
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::operator*=(const Matrix3<fpType> &mat) {
	x = this->dot(Vector3<fpType>(mat[0]));
	y = this->dot(Vector3<fpType>(mat[1]));
	z = this->dot(Vector3<fpType>(mat[2]));
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator*(const Matrix3<fpType> &mat) const {
	return Vector3<fpType>(this->dot(Vector3<fpType>(mat[0])), //
	this->dot(Vector3<fpType>(mat[1])), //
	this->dot(Vector3<fpType>(mat[2])));
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator*(const Vector3<fpType> &vec) const {
	return Matrix3<fpType>(*this) * vec;
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::operator*=(const fpType &s) {
	x *= s;
	y *= s;
	z *= s;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator*(const fpType &s) const {
	Vector3<fpType> result(data);
	result *= s;
	return result;
}

/****************************************************************************/
template<class fpType>
void Vector3<fpType>::operator/=(const fpType &s) {
	x /= s;
	y /= s;
	z /= s;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Vector3<fpType>::operator/(const fpType &s) const {
	Vector3<fpType> result(data);
	result /= s;
	return result;
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> operator*(const double &s, const Vector3<fpType> &vec) {
	return vec * s;
}

// Define the specializations that will be usable
template
Vector3<double> operator*<double>(const double &scale,
		const Vector3<double> &mat);
template
Vector3<float> operator*<float>(const double &scale, const Vector3<float> &mat);

/****************************************************************************/
template<class fpType>
Vector3<fpType> operator*(const float &s, const Vector3<fpType> &vec) {
	return vec * s;
}

// Define the specializations that will be usable
template
Vector3<double> operator*<double>(const float &scale,
		const Vector3<double> &mat);
template
Vector3<float> operator*<float>(const float &scale, const Vector3<float> &mat);

/****************************************************************************/
template<class fpType>
std::ostream& operator<<(std::ostream &stream, const Vector3<fpType> &vec) {
	return stream << vec.x << " " << vec.y << " " << vec.z;
}

// Define the specializations that will be usable
template
std::ostream& operator<<<double>(std::ostream &stream,
		const Vector3<double> &mat);
template
std::ostream& operator<<<float>(std::ostream &stream,
		const Vector3<float> &mat);

/****************************************************************************/
template<class fpType>
std::istream& operator>>(std::istream &stream, Vector3<fpType> &vec) {
	return stream >> vec.x >> vec.y >> vec.z;
}

// Define the specializations that will be usable
template
std::istream& operator>><double>(std::istream &stream, Vector3<double> &mat);
template
std::istream& operator>><float>(std::istream &stream, Vector3<float> &mat);

/****************************************************************************/
// Define the specializations that will be usable
template class Vector3<double> ;
template class Vector3<float> ;

} /* namespace algebra3d */

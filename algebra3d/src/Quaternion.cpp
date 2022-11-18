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
#include<algebra3d/Quaternion.h>
#include<algebra3d/Vector3.h>
#include<algebra3d/Matrix3.h>
#include<cstring>

namespace algebra3d {

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion() :
		x(data[0]), y(data[1]), z(data[2]), w(data[3]) {
	memset(data, 0, 3 * sizeof(fpType));
	w = 1.0;
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion(const fpType (&data)[4]) :
		x(this->data[0]), y(this->data[1]), z(this->data[2]), w(this->data[3]) {
	memcpy(this->data, data, 4 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion(fpType x, fpType y, fpType z, fpType w) :
		x(data[0]), y(data[1]), z(data[2]), w(data[3]) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->w = w;
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion(const Quaternion<fpType> &quat) :
		x(data[0]), y(data[1]), z(data[2]), w(data[3]) {
	memcpy(this->data, quat.data, 4 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion(const Vector3<fpType> &vec) :
		x(data[0]), y(data[1]), z(data[2]), w(data[3]) {
	// Cos of half angle a is always positive
	// Hence using sin(a)/cos(a) = tan(a) and sin(a)^2 + cos(a)^2 = 1
	w = sqrt(4.0 / (vec.dot(vec) + 4.0));

	x = 0.5 * vec.x * w;
	y = 0.5 * vec.y * w;
	z = 0.5 * vec.z * w;
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion(Vector3<fpType> axis, fpType angle) :
		x(data[0]), y(data[1]), z(data[2]), w(data[3]) {
	angle = normalizeAngle(angle) * 0.5;
	if ((algebra3d::compare<fpType>(axis.x, 0.0) == 0
			&& algebra3d::compare<fpType>(axis.y, 0.0) == 0
			&& algebra3d::compare<fpType>(axis.z, 0.0) == 0)
			|| (algebra3d::compare<fpType>(angle, (fpType) 0.0) == 0)) {
		x = y = z = 0.0;
		w = 1.0;
	} else {
		fpType sinhalf = sin(angle);
		// The angle was normalized and divided by 2
		// Hence the cosine is positive
		fpType coshalf = sqrt(1 - sinhalf * sinhalf);
		fpType compar = algebra3d::compare<fpType>(coshalf, (fpType) 0.0);
		axis.normalize();

		if (compar != 0) {
			x = compar * axis.x * sinhalf;
			y = compar * axis.y * sinhalf;
			z = compar * axis.z * sinhalf;
			w = compar * coshalf;
		} else {
			x = axis.x * sinhalf;
			y = axis.y * sinhalf;
			z = axis.z * sinhalf;
			w = coshalf;
		}
	}
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::Quaternion(const Matrix3<fpType> &mat) :
		x(data[0]), y(data[1]), z(data[2]), w(data[3]) {

	fpType mult;
	// Allways choose the version with positive w component
	int signW;

	// Compute 2 w*w
	w = 1.0 + mat[0][0] + mat[1][1] + mat[2][2];
	// Choose a component not close to 0
	if (w >= 0.5) {
		w = 0.5 * sqrt(w);
		mult = 0.25 / w;
		x = mult * (mat[1][2] - mat[2][1]);
		y = mult * (mat[2][0] - mat[0][2]);
		z = mult * (mat[0][1] - mat[1][0]);
	} else {
		// Compute 2 x*x
		x = 1.0 + mat[0][0] - mat[1][1] - mat[2][2];
		if (x >= 0.5) {
			x = 0.5 * sqrt(x);
			mult = 0.25 / x;
			w = mult * (mat[1][2] - mat[2][1]);
			signW = signStar<fpType>(w);
			x *= signW;
			y = signW * mult * (mat[1][0] + mat[0][1]);
			z = signW * mult * (mat[2][0] + mat[0][2]);
			w *= signW;

		} else {
			// Compute 2 y*y
			y = 1.0 - mat[0][0] + mat[1][1] - mat[2][2];
			if (y >= 0.5) {
				y = 0.5 * sqrt(y);
				mult = 0.25 / y;
				w = mult * (mat[2][0] - mat[0][2]);
				signW = signStar<fpType>(w);
				x = signW * mult * (mat[0][1] + mat[1][0]);
				y *= signW;
				z = signW * mult * (mat[2][1] + mat[1][2]);
				w *= signW;
			} else {
				// Compute 2 z*z
				z = 0.5 * sqrt(1.0 - mat[0][0] - mat[1][1] + mat[2][2]);
				mult = 0.25 / z;
				w = mult * (mat[0][1] - mat[1][0]);
				signW = signStar<fpType>(w);
				x = signW * mult * (mat[0][2] + mat[2][0]);
				y = signW * mult * (mat[1][2] + mat[2][1]);
				z *= signW;
				w *= signW;
			}
		}
	}
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType>::~Quaternion() {
	// Does nothing
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::set(const fpType *data) {
	memcpy(this->data, data, 4 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::set(const Matrix3<fpType> &mat) {
	*this = Quaternion<fpType>(mat);
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::get(fpType *storage) const {
	memcpy(storage, this->data, 4 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
fpType Quaternion<fpType>::getAngle() const {
	return normalizeAngle(2.0 * acos(sign<fpType>(w) * w));
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Quaternion<fpType>::getAxis() const {
	// If the w component represents -cos(angle/2) the sign of the others
	// needs to be flipped
	int signStarW = signStar<fpType>(w);
	return Vector3<fpType>(signStarW * x, signStarW * y, signStarW * z).getUnit();
}

/****************************************************************************/
template<class fpType>
fpType Quaternion<fpType>::getNorm() const {
	return sqrt(x * x + y * y + z * z + w * w);
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::normalize() {
	fpType normInv = 1.0 / getNorm();
	x *= normInv;
	y *= normInv;
	z *= normInv;
	w *= normInv;
}

/****************************************************************************/
template<class fpType>
bool Quaternion<fpType>::isUnit() const {
	return algebra3d::compare<fpType>(getNorm(), (fpType) 1.0) == 0;
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType> Quaternion<fpType>::getSqrtRotation() const {
	// Sing with the convention of 1 for 0
	int signStarW = signStar<fpType>(w);

	double scale = 1.0 / sqrt(2.0 * (signStarW * w + 1));
	return Quaternion<fpType>(scale * x, scale * y, scale * z,
			scale * (w + signStarW));
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType> Quaternion<fpType>::getSqrtRotationTowards(
		const Quaternion<fpType> &quat) const {
	return (this->inv() * quat).getSqrtRotation();
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::operator=(const Quaternion<fpType> &vec) {
	memcpy(data, vec.data, 4 * sizeof(fpType));
}

/****************************************************************************/
template<class fpType>
const fpType &Quaternion<fpType>::operator[](int ind) const {
	ind = (ind < 0) ? 0 : ((ind > 3) ? 3 : ind);
	return data[ind];
}

/****************************************************************************/
template<class fpType>
fpType &Quaternion<fpType>::operator[](int ind) {
	ind = (ind < 0) ? 0 : ((ind > 3) ? 3 : ind);
	return data[ind];
}

/****************************************************************************/
template<class fpType>
bool Quaternion<fpType>::compareRotations(const Quaternion<fpType> &quat,
		fpType precision) const {
	return (algebra3d::compare<fpType>(x, quat.x, precision) == 0
			&& algebra3d::compare<fpType>(y, quat.y, precision) == 0
			&& algebra3d::compare<fpType>(z, quat.z, precision) == 0
			&& algebra3d::compare<fpType>(w, quat.w, precision) == 0)
			|| (algebra3d::compare<fpType>(-x, quat.x, precision) == 0
					&& algebra3d::compare<fpType>(-y, quat.y, precision) == 0
					&& algebra3d::compare<fpType>(-z, quat.z, precision) == 0
					&& algebra3d::compare<fpType>(-w, quat.w, precision) == 0);
}

/****************************************************************************/
template<class fpType>
bool Quaternion<fpType>::compare(const Quaternion<fpType> &quat,
		const fpType &precision) const {
	return (algebra3d::compare(x, quat.x, precision) == 0)
			&& (algebra3d::compare(y, quat.y, precision) == 0)
			&& (algebra3d::compare(z, quat.z, precision) == 0)
			&& (algebra3d::compare(w, quat.w, precision) == 0);
}

/****************************************************************************/
template<class fpType>
bool Quaternion<fpType>::operator==(const Quaternion<fpType> &quat) const {
	return compare(quat);
}

/****************************************************************************/
template<class fpType>
bool Quaternion<fpType>::operator!=(const Quaternion<fpType> &quat) const {
	return !((*this) == quat);
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType> Quaternion<fpType>::inv() const {
	// Invert the axis
	// This way this times this inverse gives (0, 0, 0, 1)
	return Quaternion<fpType>(-x, -y, -z, w);
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::operator*=(const Quaternion<fpType> &quat) {
	fpType newValues[4];
	newValues[0] = w * quat.x + x * quat.w + y * quat.z - z * quat.y;
	newValues[1] = w * quat.y - x * quat.z + y * quat.w + z * quat.x;
	newValues[2] = w * quat.z + x * quat.y - y * quat.x + z * quat.w;
	newValues[3] = w * quat.w - x * quat.x - y * quat.y - z * quat.z;
	set(newValues);
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType> Quaternion<fpType>::operator*(
		const Quaternion<fpType> &quat) const {
	Quaternion<fpType> result(*this);
	result *= quat;
	return result;
}

/****************************************************************************/
template<class fpType>
void Quaternion<fpType>::scaleAngle(const fpType &mult) {
	if (algebra3d::compare<fpType>(x, 0.0) != 0
			|| algebra3d::compare<fpType>(y, 0.0) != 0
			|| algebra3d::compare<fpType>(z, 0.0) != 0) {
		*this = Quaternion<fpType>(getAxis(),
				normalizeAngle(getAngle() * mult));
	}
}

/****************************************************************************/
template<class fpType>
Quaternion<fpType> Quaternion<fpType>::getScaledAngleRotation(
		const fpType &mult) const {
	if (algebra3d::compare<fpType>(x, 0.0) == 0
			&& algebra3d::compare<fpType>(y, 0.0) == 0
			&& algebra3d::compare<fpType>(z, 0.0) == 0) {
		return Quaternion(0.0, 0.0, 0.0, w);
	}
	return Quaternion<fpType>(getAxis(), getAngle() * mult);
}

/****************************************************************************/
template<class fpType>
Vector3<fpType> Quaternion<fpType>::operator*(
		const Vector3<fpType> &vec) const {
	return Matrix3<fpType>(*this) * vec;
}

/****************************************************************************/
template<class fpType>
std::ostream& operator<<(std::ostream &stream, const Quaternion<fpType> &quat) {
	return stream << quat.x << " " << quat.y << " " << quat.z << " " << quat.w;
}

// Define the specializations that will be usable
template
std::ostream& operator<<<double>(std::ostream &stream,
		const Quaternion<double> &mat);
template
std::ostream& operator<<<float>(std::ostream &stream,
		const Quaternion<float> &mat);

/****************************************************************************/
template<class fpType>
std::istream& operator>>(std::istream &stream, Quaternion<fpType> &quat) {
	return stream >> quat.x >> quat.y >> quat.z >> quat.w;
}

// Define the specializations that will be usable
template
std::istream& operator>><double>(std::istream &stream, Quaternion<double> &mat);
template
std::istream& operator>><float>(std::istream &stream, Quaternion<float> &mat);

// Define the specializations that will be usable
template class Quaternion<double> ;
template class Quaternion<float> ;

} /* namespace algebra3d */

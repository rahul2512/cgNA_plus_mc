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
 * @file Quaternion.h
 *
 * Declaration of the Quaternion class.
 */

#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include<algebra3d/global.h>

namespace algebra3d {

/* Forward declarations */
template<class fpType>
class Matrix3;
template<class fpType>
class Vector3;

/**
 * @brief A class to store quaternions.
 *
 * The main goal is to represent rotations, however the functionality is not
 * limited to that.
 *
 * @author Jarek Glowacki
 */
template<class fpType>
class Quaternion {
	// Friendship declarations
	template<class fpT>
	friend std::istream& operator>>(std::istream &stream,
			Quaternion<fpT> &quat);

protected:
	/** The quaternion data. */
	fpType data[4];

public:
	/** Reference to the first element of #data. */
	fpType &x;
	/** Reference to the second element of #data. */
	fpType &y;
	/** Reference to the third element of #data. */
	fpType &z;
	/** Reference to the fourth element of #data. */
	fpType &w;

public:
	/**
	 * Creates an instance of a identity quaternion.
	 */
	Quaternion();

	/**
	 * Creates an instance using the provided data.
	 *
	 * @param[in] data Reference to the data to initialize the quaternion
	 * with.
	 */
	Quaternion(const fpType(&data)[4]);

	/**
	 * Initializes the vector with the provided values.
	 *
	 * @param[in] x The first element.
	 * @param[in] y The second element.
	 * @param[in] z The third element.
	 * @param[in] w The fourth element.
	 */
	Quaternion(fpType x, fpType y, fpType z, fpType w);

	/**
	 * Copies elements of the provided quaternion into this one.
	 *
	 * @param[in] quat The quaternion to copy.
	 */
	Quaternion(const Quaternion<fpType> &quat);

	/**
	 * Copies the provided rotation represented by a Cayley vector.
	 *
	 * @param[in] vec The Cayley vector to copy.
	 */
	Quaternion(const Vector3<fpType> &vec);

	/**
	 * Creates a quaternion for the provided rotation given by axis and angle.
	 *
	 * @param[in] axis The rotation axis (not necessary a unit vector).
	 * @param[in] angle The angle of rotation.
	 */
	Quaternion(Vector3<fpType> axis, fpType angle);

	/**
	 * Copies the provided rotation represented by a matrix.
	 *
	 * @param[in] mat The matrix to copy.
	 */
	Quaternion(const Matrix3<fpType> &mat);

	/**
	 * Does nothing. Virtual classes require virtual destructors.
	 */
	virtual ~Quaternion();

	/**
	 * Set the quaternion data to the values in the provided array. First 4
	 * positions from the given pointer are read.
	 *
	 * @param[in] data The array with values to set.
	 */
	virtual void set(const fpType *data);

	/**
	 * Set the quaternion data to represent the rotation given as a matrix.
	 *
	 * @param[in] mat The rotation matrix.
	 */
	virtual void set(const Matrix3<fpType> &mat);

	/**
	 * Copies the quaternion data to the provided array. First 4 positions
	 * from the given pointer are written).
	 *
	 * @param[out] storage The array to copy the data to.
	 */
	virtual void get(fpType *storage) const;

	/**
	 * Returns the angle of the rotation represented by this quaternion. If
	 * the quaternions is not a unit one the value is meaningless.
	 *
	 * @return The angle of the rotation represented by this quaternion.
	 */
	virtual fpType getAngle() const;

	/**
	 * Returns the axis of the rotation represented by this quaternion. If
	 * the quaternions is not a unit one the value is meaningless.
	 *
	 * @return The axis of the rotation represented by this quaternion.
	 */
	virtual Vector3<fpType> getAxis() const;

	/**
	 * Returns the L2 norm of the quaternion.
	 *
	 * @return The L2 norm of the quaternion.
	 */
	virtual fpType getNorm() const;

	/**
	 * Normalizes the quaternion (divides all elements by the norm).
	 */
	virtual void normalize();

	/**
	 * Returns `true` if this is a unit quaternion,
	 * `false` otherwise.
	 *
	 * @return `true` if this is a unit quaternion,
	 * `false` otherwise.
	 */
	virtual bool isUnit() const;

	/**
	 * Returns the quaternion representation of half of the rotation
	 * represented by this instance.
	 *
	 * @return the quaternion representation of half of the rotation
	 * represented by this instance.
	 */
	virtual Quaternion<fpType> getSqrtRotation() const;

	/**
	 * Computes the rotation that is half the transition between this one and
	 * the provided one. That is to say that:
	 *
	 *     Quaternion q1, q2;
	 *     ...
	 *     Quaternion trans = q1.getSqrtRotationTowards(q2);
	 *     // The following is true:
	 *     q1 * trans * trans == q2;
	 *
	 * Returns it as a Quaternion.
	 *
	 * @param[in] quat The other rotation.
	 * @return the Quaternion that represents the square root of transition.
	 */
	virtual Quaternion<fpType> getSqrtRotationTowards(
			const Quaternion<fpType> &quat) const;

	/**
	 * Copies elements of the provided quaternion into this one.
	 *
	 * @param[in] quat The quaternion to copy.
	 */
	virtual void operator=(const Quaternion<fpType> &quat);

	/**
	 * Returns a reference to the element of the quaternion with the provided
	 * index. If the index is incorrect the closest correct one is used is
	 * returned. Version for constant instances.
	 *
	 * @param[in] ind Index of the element to return.
	 * @return Reference to the selected element.
	 */
	virtual const fpType &operator[](int ind) const;

	/**
	 * Returns a reference to the element of the quaternion with the provided
	 * index. If the index is incorrect the closest correct one is used is
	 * returned.
	 *
	 * @param[in] ind Index of the element to return.
	 * @return Reference to the selected element.
	 */
	virtual fpType &operator[](int ind);

	/**
	 * Checks if the rotation represented by this quaternion is identical with
	 * the one represented by the argument.
	 *
	 * @param[in] quat The quaternion to compare this with.
	 * @param[in] precision The precision to use; by default the global
	 * PRECISION is used.
	 *
	 * @return `true` if the rotations are identical,
	 * `false` otherwise.
	 */
	virtual bool compareRotations(const Quaternion<fpType> &quat,
			fpType precision = PRECISION) const;

	/**
	 * Checks if all elements of the provided quaternion are the same as in
	 * this one up to the provided precision. The comparison is simply checking
	 * if the difference between numbers is smaller or equal the precision.
	 *
	 * @param[in] quat The quaternion to compare this with.
	 * @param[in] precision The precision to use; by default the global
	 * PRECISION is used.
	 *
	 * @return `true` if all entries are equal up to the precision,
	 * `false` otherwise.
	 */
	virtual bool compare(const Quaternion<fpType> &quat,
			const fpType &precision = PRECISION) const;

	/**
	 * Checks if all elements of the provided quaternion are the same as this.
	 *
	 * @param[in] quat The quaternion to compare this with.
	 *
	 * @return `true` if all entries are equal,
	 * `false` otherwise.
	 */
	virtual bool operator==(const Quaternion<fpType> &quat) const;

	/**
	 * Comparing operator. Checks if any element of this quaternion differs
	 * from respective element of the parameter.
	 *
	 * @param[in] quat The quaternion vector to compare this with.
	 * @return `false` if all entries are equal,
	 * `true` otherwise.
	 */
	virtual bool operator!=(const Quaternion<fpType> &quat) const;

	/**
	 * Returns a quaternion that represents inverse of this rotation, i.e.
	 * same angle, negative axis the axis of this rotation. For a non-unit
	 * quaternion may result in an error.
	 *
	 * @return A quaternion that represents a rotation that is inverse of this
	 * one.
	 */
	virtual Quaternion<fpType> inv() const;

	/**
	 * Post-multiplies this quaternion by the provided one.
	 *
	 * @param[in] quat The quaternion to post-multiply this by.
	 */
	virtual void operator*=(const Quaternion<fpType> &quat);

	/**
	 * Returns the product of this quaternion and the provided one.
	 *
	 * @param[in] quat The matrix to multiply this vector.
	 * @return The product of this quaternion and the provided one.
	 */
	virtual Quaternion<fpType> operator*(const Quaternion<fpType> &quat) const;

	/**
	 * Changes the angle of rotation represented by this quaternion by
	 * multiplying it by the provided number.
	 *
	 * @param mult The value to multiply the angle by.
	 */
	virtual void scaleAngle(const fpType &mult);

	/**
	 * Returns a quaternion representing a rotation around the same axis by an
	 * angle that is a given multiply of this rotation angle.
	 *
	 * @param mult The value to multiply the angle by.
	 */
	virtual Quaternion<fpType> getScaledAngleRotation(const fpType &mult) const;

	/**
	 * Returns the effect of rotating the vector according to this quaternion.
	 * If this is not a unit quaternion this might result in an error.
	 *
	 * @param[in] vec The vector to transform according to this quaternion.
	 * @return The effect of rotating the vector according to this quaternion.
	 */
	virtual Vector3<fpType> operator*(const Vector3<fpType> &vec) const;
};

/**
 * Outputs the provided quaternion to the provided stream.
 *
 * @param[in] stream The stream to use.
 * @param[in] quat The quaternion to output.
 * @return The stream.
 */
template<class fpType>
std::ostream& operator<<(std::ostream &stream, const Quaternion<fpType> &quat);

/**
 * Reads data for the provided quaternion from the provided stream.
 *
 * @param[in] stream The stream to use.
 * @param[in] quat The quaternion to fill in.
 * @return The stream.
 */
template<class fpType>
std::istream& operator>>(std::istream &stream, Quaternion<fpType> &quat);

/**
 * @typedef #algebra3d#Quaternion \double\> #algebra3d#QuaternionD
 *
 * An abbreviated name for #algebra3d#Quaternion \<double\>
 */
typedef Quaternion<double> QuaternionD;

/**
 * @typedef #algebra3d#Quaternion \<double\> #algebra3d#QuaternionF
 *
 * An abbreviated name for #algebra3d#Quaternion \<float\>
 */
typedef Quaternion<float> QuaternionF;

}  /* namespace algebra3d */

#endif // _QUATERNION_H_

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
 * @file Vector3.h
 *
 * Declaration of the Vector class.
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include<algebra3d/global.h>

namespace algebra3d {

/* Forward declarations */
template<class fpType>
class Matrix3;
template<class fpType>
class Quaternion;

/**
 * @brief A class to store vectors of 3 elements.
 *
 * @author Jarek Glowacki
 */
template<class fpType>
class Vector3 {
	// Friendship declarations
	template<class fpT>
	friend std::istream& operator>>(std::istream &stream, Vector3<fpT> &vec);

protected:
	/** The vector data. */
	fpType data[3];

public:
	/** Reference to the first element of #data. */
	fpType &x;
	/** Reference to the second element of #data. */
	fpType &y;
	/** Reference to the third element of #data. */
	fpType &z;

public:
	/**
	 * Creates an instance of a 0 vector.
	 */
	Vector3();

	/**
	 * Creates an instance using the provided data.
	 *
	 * @param[in] data Reference to the data to initialize the vector with.
	 */
	Vector3(const fpType * const data);

	/**
	 * Initializes the vector with the provided values.
	 *
	 * @param[in] x First coordinate.
	 * @param[in] y Second coordinate.
	 * @param[in] z Third coordinate.
	 */
	Vector3(fpType x, fpType y, fpType z);

	/**
	 * Copies elements of the provided vector into this one.
	 *
	 * @param[in] vec The vector to copy.
	 */
	Vector3(const Vector3<fpType> &vec);

	/**
	 * Creates a Cayley vector for the provided rotation given by axis and
	 * angle.
	 *
	 * @param[in] axis The rotation axis (not necessary a unit vector).
	 * @param[in] angle The angle of rotation.
	 */
	Vector3(const Vector3<fpType> &axis, fpType angle);

	/**
	 * Creates a Cayley vector for the provided rotation given by a
	 * quaternion.
	 *
	 * @param[in] quat The quaternion representing the rotation.
	 */
	Vector3(const Quaternion<fpType> &quat);

	/**
	 * Creates a Cayley vector for the provided rotation given by a matrix.
	 * @note Watch out for the singular case with trace of the provided matrix
	 * equal to -1 where the norm of the Cayley vector tends to infinity.
	 *
	 * @param[in] mat The matrix of the rotation.
	 */
	Vector3(const Matrix3<fpType> &mat);

	/**
	 * Does nothing. Virtual classes require virtual destructors.
	 */
	virtual ~Vector3();

	/**
	 * Set the vector data to the values in the provided array. First 3
	 * positions from the given pointer are read.
	 *
	 * @param[in] data The array with values to set.
	 */
	virtual void set(const fpType *data);

	/**
	 * Copies the vector data to the provided array. First 3 positions from the
	 * given pointer are written).
	 *
	 * @param[out] storage The array to copy the data to.
	 */
	virtual void get(fpType *storage) const;

	/**
	 * Returns the dot product of this vector and the provided one.
	 *
	 * @param[in] vec The vector to multiply this by.
	 * @return The dot product of this vector and the provided one.
	 */
	virtual fpType dot(const Vector3<fpType> &vec) const;

	/**
	 * Returns the cross product of this vector with the provided one.
	 *
	 * @param[in] vec The other vector of the cross product.
	 * @return The cross product of this vector with the provided one.
	 */
	virtual Vector3<fpType> cross(const Vector3<fpType> &vec) const;

	/**
	 * Returns the outer product of this vector and the provided one.
	 *
	 * @param[in] vec The vector to multiply this by.
	 * @return The outer product of this vector and the provided one.
	 */
	virtual Matrix3<fpType> outer(const Vector3<fpType> &vec) const;

	/**
	 * Returns the cross product matrix of this vector, i.e. a cross product
	 * operator for the vector.
	 *
	 * @return The cross product matrix of this vector
	 */
	virtual Matrix3<fpType> crossMatrix() const;

	/**
	 * Returns the L2 norm of the vector (`sqrt(*this * *this)`).
	 *
	 * @return The L2 norm of the vector (`sqrt(*this * *this)`).
	 */
	virtual fpType getNorm() const;

	/**
	 * Normalizes the vector to unit length. Divides the vector by
	 * `sqrt(*this * *this)`.
	 */
	virtual void normalize();

	/**
	 * Returns `true` if all components are 0 (up to the precision), `false`
	 * otherwise.
	 *
	 * @return `true` if all components are 0 (up to the precision), `false`
	 * otherwise.
	 */
	bool isZero() const;

	/**
	 * Returns a unit vector in the direction of this one. For infinite and
	 * zero vectors their copies are returned.
	 *
	 * @return A unit vector in the direction of this one.
	 */
	Vector3<fpType> getUnit() const;

	/**
	 * Returns the angle of rotation represented by this Cayley vector.
	 *
	 * @return The angle of rotation represented by this Cayley vector.
	 */
	fpType getAngle() const;

	/**
	 * Returns the Cayley vector representation of half of the rotation
	 * represented by this instance.
	 *
	 * @return the Cayley vector representation of half of the rotation
	 * represented by this instance.
	 */
	Vector3<fpType> getSqrtRotation() const;

	/**
	 * Copies elements of the provided vector into this one.
	 *
	 * @param[in] vec The vector to copy.
	 */
	virtual void operator=(const Vector3<fpType> &vec);

	/**
	 * Returns a reference to the element of the vector with the provided
	 * index. If the index is incorrect the closest correct one is used is
	 * returned. Version for constant instances.
	 *
	 * @param[in] ind Index of the element to return.
	 * @return Reference to the selected element.
	 */
	virtual const fpType &operator[](int ind) const;

	/**
	 * Returns a reference to the element of the vector with the provided
	 * index. If the index is incorrect the closest correct one is used is
	 * returned.
	 *
	 * @param[in] ind Index of the element to return.
	 * @return Reference to the selected element.
	 */
	virtual fpType &operator[](int ind);

	/**
	 * Checks if all elements of the provided vector are the same as in this
	 * one up to the provided precision. The comparison is simply checking
	 * if the difference between numbers is smaller or equal the precision.
	 *
	 * @param[in] vec The vector to compare this with.
	 * @param[in] precision The precision to use; by default the global
	 * PRECISION is used.
	 *
	 * @return `true` if all entries are equal up to the precision,
	 * `false` otherwise.
	 */
	virtual bool compare(const Vector3<fpType> &vec, const fpType &precision =
			PRECISION) const;

	/**
	 * Comparing operator. Checks if all elements of the provided vector are
	 * the same as this.
	 *
	 * @param[in] vec The vector to compare this with.
	 *
	 * @return `true` if all entries are equal,
	 * `false` otherwise.
	 */
	virtual bool operator==(const Vector3<fpType> &vec) const;

	/**
	 * Comparing operator. Checks if any element of this vector differs from
	 * respective element of the parameter.
	 *
	 * @param[in] vec The vector to compare this with.
	 * @return `false` if all entries are equal,
	 * `true` otherwise.
	 */
	virtual bool operator!=(const Vector3<fpType> &vec) const;

	/**
	 * Adds the provided vector to this one.
	 *
	 * @param[in] vec The vector to add.
	 */
	virtual void operator+=(const Vector3<fpType> &vec);

	/**
	 * Returns the sum of this vector and the provided one.
	 *
	 * @param[in] vec The other argument of the sum.
	 * @return The sum of this vector and the provided one.
	 */
	virtual Vector3<fpType> operator+(const Vector3<fpType> &vec) const;

	/**
	 * Returns a copy of this vector multiplied by -1.
	 *
	 * @return A copy of this vector multiplied by -1.
	 */
	virtual Vector3<fpType> operator-() const;

	/**
	 * Subtracts the provided vector from this one.
	 *
	 * @param[in] vec The vector to subtract.
	 */
	virtual void operator-=(const Vector3<fpType> &vec);

	/**
	 * Returns the result of subtracting the provided vector from this one.
	 *
	 * @param[in] vec The matrix to subtract.
	 * @return The result of subtracting the provided vector from this one.
	 */
	virtual Vector3<fpType> operator-(const Vector3<fpType> &vec) const;

	/**
	 * Post-multiplies this vector by the provided matrix (v * m - vector
	 * treated as a 1x3).
	 *
	 * @param[in] mat The matrix to multiply this vector by.
	 * @return The product of this vector and the provided matrix
	 * (post-multiplication).
	 */
	void operator*=(const Matrix3<fpType> &mat);

	/**
	 * Returns the product of this vector and the provided matrix
	 * (v * m - vector treated as 1x3).
	 *
	 * @param[in] mat The matrix to multiply this vector by.
	 * @return The product of this vector and the provided matrix
	 * (post-multiplication).
	 */
	virtual Vector3<fpType> operator*(const Matrix3<fpType> &mat) const;

	/**
	 * Returns the effect of rotating the provided vector according to this
	 * Cayley vector.
	 *
	 * @param[in] vec The vector to transform according to this Cayley vector.
	 * @return The effect of rotating the provided vector according to this
	 * Cayley vector.
	 */
	virtual Vector3<fpType> operator*(const Vector3<fpType> &vec) const;

	/**
	 * Scales a this vector by the provided number.
	 *
	 * @param[in] s The scaling factor.
	 */
	virtual void operator*=(const fpType &s);

	/**
	 * Scales a copy of this vector by the provided number and returns it.
	 *
	 * @param[in] s The scaling factor.
	 * @return The scaled vector.
	 */
	virtual Vector3<fpType> operator*(const fpType &s) const;

	/**
	 * Scales a this matrix by the inverse of the provided number.
	 *
	 * @param[in] s The inverse scaling factor.
	 */
	virtual void operator/=(const fpType &s);

	/**
	 * Scales a copy of this vector by the inverse of the provided number
	 * and returns it.
	 *
	 * @param[in] s The inverse scaling factor.
	 * @return The scaled matrix.
	 */
	virtual Vector3<fpType> operator/(const fpType &s) const;
};

/**
 * Scales the provided vector by the provided double number.
 *
 * @param[in] scale The scaling factor.
 * @param[in] vec Vector to scale.
 * @return The scaled matrix.
 */
template<class fpType>
Vector3<fpType> operator*(const double &scale, const Vector3<fpType> &vec);

/**
 * Scales the provided vector by the provided float number.
 *
 * @param[in] scale The scaling factor.
 * @param[in] vec Vector to scale.
 * @return The scaled matrix.
 */
template<class fpType>
Vector3<fpType> operator*(const float &scale, const Vector3<fpType> &vec);

/**
 * Outputs the provided vector to the provided stream.
 *
 * @param[in] stream The stream to use.
 * @param[in] vec The vector to output.
 * @return The stream.
 */
template<class fpType>
std::ostream& operator<<(std::ostream &stream, const Vector3<fpType> &vec);

/**
 * Reads data for the provided vector from the provided stream.
 *
 * @param[in] stream The stream to use.
 * @param[in] vec The vector to fill in.
 * @return The stream.
 */
template<class fpType>
std::istream& operator>>(std::istream &stream, Vector3<fpType> &vec);

/**
 * @typedef #algebra3d#Vector3 \<double\> #algebra3d#Vector3D
 *
 * An abbreviated name for #algebra3d#Vector3 \<double\>
 */
typedef Vector3<double> Vector3D;

/**
 * @typedef #algebra3d#Vector3 \<double\> #algebra3d#Vector3F
 *
 * An abbreviated name for #algebra3d#Vector3 \<float\>
 */
typedef Vector3<float> Vector3F;

} /* namespace algebra3d */

#endif // _VECTOR_H_

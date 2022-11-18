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
 * @file Matrix3.h
 * Declaration of the Matrix class.
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include<algebra3d/global.h>

namespace algebra3d {

/* Forward declarations */
template<class fpType>
class Vector3;
template<class fpType>
class Quaternion;

/**
 * @brief A class to store 3x3 rotation matrices.
 *
 * Column first storage convention is used here. It is not limited to SO(3)
 * matrices.
 *
 * @author Jarek Glowacki
 */
template<class fpType>
class Matrix3 {
	// Friendship declarations
	template<class fpT>
	friend std::istream& operator>>(std::istream &stream, Matrix3<fpT> &mat);

protected:
	/**
	 * Indices for computing adjugate matrix. The indices are:
	 * - column of the final matrix
	 * - row of the final matrix
	 * - element of the 2 by 2 det (first the ones that come positive, than the
	 * ones that come negative)
	 * - a pair of column, row indices
	 */
	static const int ADJUGATE_INDICES[3][3][4][2];

	/** The matrix data. */
	fpType data[3][3];

protected:
	/**
	 * Initializes the instance of a matrix. Used by constructors. The
	 * parameter defines if identity (`true`) or zero matrix
	 * (`false`) is requested.
	 *
	 * @param[in] identity identity `true` or zero matrix
	 * `false` is requested; `false` by default.
	 */
	virtual void initialize(bool identity = false);

public:
	/**
	 * Creates an instance of a matrix. The parameter defines if identity
	 * (`true`) or zero matrix (`false`) is requested.
	 *
	 * @param[in] identity identity `true` or zero matrix
	 * `false` is requested; `false` by default.
	 */
	Matrix3(bool identity = false);

	/**
	 * Creates an instance using the provided data. Data is stored by
	 * columns (Matrix#data[2][0] is third column first row)
	 *
	 * @param[in] data Reference to the data to set the matrix to.
	 */
	Matrix3(const fpType(&data)[3][3]);

	/**
	 * Creates an instance using the provided data. First 9 positions. are
	 * read: first 3 - first column, second 3 - second column, third 3 - third
	 * column.
	 *
	 * @param[in] data Reference to the data to set the matrix to.
	 */
	Matrix3(const fpType *data);

	/**
	 * Creates an instance using the provided vectors as columns. Data is
	 * stored by columns, so Matrix#data[0] becomes the fist vector,
	 * Matrix#data[1] the second and so on.
	 *
	 * @param[in] col1 The first column of the matrix.
	 * @param[in] col2 The second column of the matrix.
	 * @param[in] col3 The third column of the matrix.
	 */
	Matrix3(const Vector3<fpType> &col1, const Vector3<fpType> &col2,
			const Vector3<fpType> &col3);

	/**
	 * Creates an instance using the provided values.
	 *
	 * @param[in] m11 The element to becomeMatrix#data[0][0].
	 * @param[in] m12 The element to becomeMatrix#data[1][0].
	 * @param[in] m13 The element to becomeMatrix#data[2][0].
	 * @param[in] m21 The element to becomeMatrix#data[0][1].
	 * @param[in] m22 The element to becomeMatrix#data[1][1].
	 * @param[in] m23 The element to becomeMatrix#data[2][1].
	 * @param[in] m31 The element to becomeMatrix#data[0][2].
	 * @param[in] m32 The element to becomeMatrix#data[1][2].
	 * @param[in] m33 The element to becomeMatrix#data[2][2].
	 */
	Matrix3(fpType m11, fpType m12, fpType m13, //
			fpType m21, fpType m22, fpType m23, //
			fpType m31, fpType m32, fpType m33);

	/**
	 * Copies elements of the provided matrix into this one.
	 *
	 * @param[in] mat The matrix to copy.
	 */
	Matrix3(const Matrix3<fpType> &mat);

	/**
	 * Creates a rotation matrix for the provided Cayley vector.
	 *
	 * @param[in] cvec The Cayley vector describing the rotation.
	 * @note Passes though a Quaternion for efficiency.
	 */
	Matrix3(const Vector3<fpType> &cvec);

	/**
	 * Creates a rotation matrix for the provided axis (not necessarily unit
	 * length) and angle (in radians).
	 *
	 * @param[in] axis The rotation axis (of any length).
	 * @param[in] angle The rotation angle in radians.
	 */
	Matrix3(Vector3<fpType> axis, fpType angle);

	/**
	 * Creates a rotation matrix for the provided quaternion.
	 *
	 * @param[in] quat The rotation axis (of any length).
	 */
	Matrix3(const Quaternion<fpType> &quat);

	/**
	 * Does nothing. Virtual classes require virtual destructors.
	 */
	virtual ~Matrix3();

	/**
	 * Sets the matrix data to the values in the provided array.
	 *
	 * @param[in] data The array with values to set.
	 */
	virtual void set(const fpType data[3][3]);

	/**
	 * Sets the matrix data to the data provided as column vectors.
	 *
	 * @param[in] col1 The first column of the matrix.
	 * @param[in] col2 The second column of the matrix.
	 * @param[in] col3 The third column of the matrix.
	 */
	virtual void set(const Vector3<fpType> &col1, const Vector3<fpType> &col2,
			const Vector3<fpType> &col3);

	/**
	 * Copies the matrix data to the provided array.
	 *
	 * @param[out] storage The array to copy the data to.
	 */
	virtual void get(fpType storage[3][3]) const;

	/**
	 * Copies the matrix data to the provided arrays. Each array is assumed to
	 * be at least of size [3].
	 *
	 * @param[out] col1 Storage for the first column of the matrix.
	 * @param[out] col2 Storage for the second column of the matrix.
	 * @param[out] col3 Storage for the third column of the matrix.
	 */
	virtual void get(fpType * const col1, fpType * const col2,
			fpType * const col3) const;

	/**
	 * Sets the column with the provided index to the provided vector.
	 *
	 * @param[in] index Index of the column to set.
	 * @param[in] col The column data.
	 */
	virtual void setColumn(int index, const Vector3<fpType> &col);

	/**
	 * Copies data of the column with the provided index to the provided array.
	 *
	 * @param[in] index Index of the column to copy.
	 * @param[out] col The storage place to copy to.
	 */
	virtual void getColumn(int index, fpType * const &col);

	/**
	 * Returns the requested column as a vector.
	 *
	 * @param[in] index Index of the column to copy.
	 */
	virtual Vector3<fpType> getColumn(int index);

	/**
	 * Returns the trace of the matrix.
	 *
	 * @return The trace of the matrix.
	 */
	fpType getTrace() const;

	/**
	 * Transposes this matrix.
	 */
	virtual void transpose();

	/**
	 * Returns the transpose of this matrix.
	 *
	 * @return The transpose of this matrix.
	 */
	virtual Matrix3 getTranspose() const;

	/**
	 * Returns the unit axis vector of the rotation this matrix represents.
	 * For identity and zero vector is returned. For rotation by PI
	 * infinity vector. If the matrix is not orthogonal this may result
	 * in an error.
	 *
	 * @return The axis of the rotation this matrix represents. For identity
	 * and zero vector is returned. For rotation by PI - infinity vector.
	 */
	Vector3<fpType> getAxis() const;

	/**
	 * Returns the angle of the rotation this matrix represents. If the matrix
	 * is not orthogonal this may result in an error.
	 *
	 * @return The angle of the rotation this matrix represents.
	 */
	fpType getAngle() const;

	/**
	 * Inverts this matrix. If the determinant is close to 0.0
	 * the matrix is set to 0.
	 */
	virtual void invert();

	/**
	 * Returns the inverse of this matrix. If the determinant is close to 0
	 * zero matrix is returned.
	 *
	 * @return The inverse of this matrix.
	 */
	virtual Matrix3 getInverse() const;

	/**
	 * Returns the determinant of this matrix.
	 *
	 * @return The determinant of this matrix.
	 */
	virtual fpType det();

	/**
	 * Returns `true` if the matrix is orthogonal,
	 * `false` otherwise.
	 *
	 * @return `true` if the matrix is orthogonal,
	 * `false` otherwise.
	 */
	virtual bool isOrthogonal();

	/**
	 * Copies elements of the provided matrix into this one.
	 *
	 * @param[in] mat The matrix to copy.
	 */
	virtual void operator=(const Matrix3 &mat);

	/**
	 * Returns a reference to the column with the provided index. If the index
	 * is incorrect the closest correct one is used is returned. Version for
	 * constant instances.
	 *
	 * @param[in] i1 Index of the column of the element to return.
	 * @return Reference to the selected column or `NULL` for
	 * incorrect index.
	 */
	virtual const fpType (&operator[](int i1) const)[3];

	/**
	 * Returns a reference to the column with the provided index. If
	 * the index is incorrect `NULL` is returned.
	 *
	 * @param[in] i1 Index of the column of the element to return.
	 * @return Reference to the selected column or `NULL` for
	 * incorrect index.
	 */
virtual	fpType (&operator[](int i1))[3];

	/**
	 * Checks if all elements of the provided matrix are the same as in this
	 * one up to the provided precision. The comparison is simply checking if
	 * the difference between numbers is smaller or equal the precision.
	 *
	 * @param[in] mat The matrix to compare this with.
	 * @param[in] precision The precision to use; by default the global
	 * PRECISION is used.
	 *
	 * @return `true` if all entries are equal up to the precision,
	 * `false` otherwise.
	 */
	virtual bool compare(const Matrix3<fpType> &mat, const fpType &precision =
			PRECISION) const;

	/**
	 * Comparing operator. Checks if all elements of the provided matrix are
	 * the same as this.
	 *
	 * @param[in] mat The matrix to compare this with.
	 *
	 * @return `true` if all entries are equal,
	 * `false` otherwise.
	 */
	virtual bool operator==(const Matrix3 &mat) const;

	/**
	 * Checks if any element of this matrix differs from
	 * respective element of the parameter.
	 *
	 * @param[in] mat The matrix to compare this with.
	 *
	 * @return `false` if all entries are equal,
	 * `true` otherwise.
	 */
	virtual bool operator!=(const Matrix3 &mat) const;

	/**
	 * Adds the provided matrix to this one.
	 *
	 * @param[in] mat The matrix to add.
	 */
	virtual void operator+=(const Matrix3 &mat);

	/**
	 * Returns the sum of this matrix and the provided one.
	 *
	 * @param[in] mat The other argument of the sum.
	 * @return The sum of this matrix and the provided one.
	 */
	virtual Matrix3 operator+(const Matrix3 &mat) const;

	/**
	 * Returns a copy of this matrix multiplied by -1.
	 *
	 * @return A copy of this matrix multiplied by -1.
	 */
	virtual Matrix3 operator-() const;

	/**
	 * Subtracts the provided matrix from this one.
	 *
	 * @param[in] mat The matrix to subtract.
	 */
	virtual void operator-=(const Matrix3 &mat);

	/**
	 * Returns the result of subtracting the provided matrix from this one.
	 *
	 * @param[in] mat The matrix to subtract.
	 * @return The result of subtracting the provided matrix from this one.
	 */
	virtual Matrix3 operator-(const Matrix3 &mat) const;

	/**
	 * Post-multiplies this matrix by the provided one.
	 *
	 * @param[in] mat The matrix to post-multiply this by.
	 */
	virtual void operator*=(const Matrix3 &mat);

	/**
	 * Post-multiplies a copy of this matrix by the provided one and returns
	 * the result.
	 *
	 * @param[in] mat The matrix to post-multiply this by.
	 * @return The result of the post-multiplication.
	 */
	virtual Matrix3 operator*(const Matrix3 &mat) const;

	/**
	 * Returns the product of this matrix and the provided vector.
	 * (m * v - vector treated as 3x1).
	 *
	 * @param[in] vec The vector to pre-multiply by this matrix.
	 * @return The result of the post-multiplication.
	 */
	virtual Vector3<fpType> operator*(const Vector3<fpType> &vec) const;

	/**
	 * Scales a this matrix by the provided number.
	 *
	 * @param[in] scale The scaling factor.
	 */
	virtual void operator*=(const fpType &scale);

	/**
	 * Scales a copy of this matrix by the provided number and returns it.
	 *
	 * @param[in] scale The scaling factor.
	 * @return The scaled matrix.
	 */
	virtual Matrix3 operator*(const fpType &scale) const;

	/**
	 * Scales a this matrix by the inverse of the provided number.
	 *
	 * @param[in] scale The inverse scaling factor.
	 */
	virtual void operator/=(const fpType &scale);

	/**
	 * Scales a copy of this matrix by the inverse of the provided number
	 * and returns it.
	 *
	 * @param[in] scale The inverse scaling factor.
	 * @return The scaled matrix.
	 */
	virtual Matrix3 operator/(const fpType &scale) const;
};

/**
 * Scales the provided matrix by the provided double number.
 *
 * @param[in] scale The scaling factor.
 * @param[in] mat Matrix to scale.
 * @return The scaled matrix.
 */
template<class fpType>
Matrix3<fpType> operator*(const double &scale, const Matrix3<fpType> &mat);

/**
 * Scales the provided matrix by the provided float number.
 *
 * @param[in] scale The scaling factor.
 * @param[in] mat Matrix to scale.
 * @return The scaled matrix.
 */
template<class fpType>
Matrix3<fpType> operator*(const float &scale, const Matrix3<fpType> &mat);

/**
 * Outputs the provided matrix to the provided stream..
 *
 * @param[in] stream The stream to use.
 * @param[in] mat The matrix to output.
 * @return The stream.
 */
template<class fpType>
std::ostream& operator<<(std::ostream &stream, const Matrix3<fpType> &mat);

/**
 * Reads data for the provided matrix from the provided stream.
 *
 * @param[in] stream The stream to use.
 * @param[in] mat The matrix to fill in.
 * @return The stream.
 */
template<class fpType>
std::istream& operator>>(std::istream &stream, Matrix3<fpType> &mat);

/**
 * @typedef #algebra3d#Matrix3 \<double\> #algebra3d#Matrix3D
 *
 * An abbreviated name for #algebra3d#Matrix3 \<double\>
 */
typedef Matrix3<double> Matrix3D;

/**
 * @typedef #algebra3d#Matrix3 \<float\> #algebra3d#Matrix3F
 *
 * An abbreviated name for #algebra3d#Matrix3 \<float\>
 */
typedef Matrix3<float> Matrix3F;

} /* namespace algebra3d */

#endif /* _MATRIX_H_ */

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.complex;

import java.io.Serializable;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.Precision;

/**
 * This class implements <a href="http://mathworld.wolfram.com/Quaternion.html">
 * quaternions</a> (Hamilton's hypercomplex numbers).
 * <p>It is a reimplementation of {@link Quaternion} using {@link CalculusFieldElement} elements.</p>
 * <p>Instance of this class are guaranteed to be immutable.</p>
 *
 * @author Alberto Foss&agrave;
 */
public final class FieldQuaternion<T extends CalculusFieldElement<T>> implements Serializable {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20230528L;

    /** First component (scalar part). */
    private final T q0;

    /** Second component (first vector part). */
    private final T q1;

    /** Third component (second vector part). */
    private final T q2;

    /** Fourth component (third vector part). */
    private final T q3;

    /**
     * Builds a quaternion from its components.
     *
     * @param a Scalar component.
     * @param b First vector component.
     * @param c Second vector component.
     * @param d Third vector component.
     */
    public FieldQuaternion(final T a, final T b, final T c, final T d) {
        this.q0 = a;
        this.q1 = b;
        this.q2 = c;
        this.q3 = d;
    }

    /**
     * Builds a quaternion from scalar and vector parts.
     *
     * @param scalar Scalar part of the quaternion.
     * @param v Components of the vector part of the quaternion.
     *
     * @throws MathIllegalArgumentException if the array length is not 3.
     */
    public FieldQuaternion(final T scalar, final T[] v)
            throws MathIllegalArgumentException {
        MathUtils.checkDimension(v.length, 3);
        this.q0 = scalar;
        this.q1 = v[0];
        this.q2 = v[1];
        this.q3 = v[2];
    }

    /**
     * Builds a pure quaternion from a vector (assuming that the scalar
     * part is zero).
     *
     * @param v Components of the vector part of the pure quaternion.
     */
    public FieldQuaternion(final T[] v) {
        this(v[0].getField().getZero(), v);
    }

    /**
     * Build a {@link FieldQuaternion} from a {@link Quaternion}.
     * @param field field for the components
     * @param q quaternion
     */
    public FieldQuaternion(final Field<T> field, final Quaternion q) {
        this(field.getZero().add(q.getQ0()), field.getZero().add(q.getQ1()), field.getZero().add(q.getQ2()), field.getZero().add(q.getQ3()));
    }

    /**
     * Returns the identity quaternion.
     * @param field field for the components
     * @return identity quaternion
     * @param <T> the type of field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> getIdentity(final Field<T> field) {
        return new FieldQuaternion<>(field, Quaternion.IDENTITY);
    }

    /**
     * Returns the zero quaternion.
     * @param field field for the components
     * @return zero quaternion
     * @param <T> the type of field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> getZero(final Field<T> field) {
        return new FieldQuaternion<>(field, Quaternion.ZERO);
    }

    /**
     * Returns the i quaternion.
     * @param field field for the components
     * @return i quaternion
     * @param <T> the type of field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> getI(final Field<T> field) {
        return new FieldQuaternion<>(field, Quaternion.I);
    }

    /**
     * Returns the j quaternion.
     * @param field field for the components
     * @return j quaternion
     * @param <T> the type of field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> getJ(final Field<T> field) {
        return new FieldQuaternion<>(field, Quaternion.J);
    }

    /**
     * Returns the k quaternion.
     * @param field field for the components
     * @return k quaternion
     * @param <T> the type of field elements
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> getK(final Field<T> field) {
        return new FieldQuaternion<>(field, Quaternion.K);
    }

    /**
     * Returns the conjugate quaternion of the instance.
     *
     * @return the conjugate quaternion
     */
    public FieldQuaternion<T> getConjugate() {
        return new FieldQuaternion<>(q0, q1.negate(), q2.negate(), q3.negate());
    }

    /**
     * Returns the Hamilton product of two quaternions.
     *
     * @param q1 First quaternion.
     * @param q2 Second quaternion.
     * @return the product {@code q1} and {@code q2}, in that order.
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> multiply(final FieldQuaternion<T> q1,
                                                                                  final FieldQuaternion<T> q2) {

        // Components of the first quaternion.
        final T q1a = q1.getQ0();
        final T q1b = q1.getQ1();
        final T q1c = q1.getQ2();
        final T q1d = q1.getQ3();

        // Components of the second quaternion.
        final T q2a = q2.getQ0();
        final T q2b = q2.getQ1();
        final T q2c = q2.getQ2();
        final T q2d = q2.getQ3();

        // Components of the product.
        final T w = q1a.multiply(q2a).subtract(q1b.multiply(q2b)).subtract(q1c.multiply(q2c)).subtract(q1d.multiply(q2d));
        final T x = q1a.multiply(q2b).add(q1b.multiply(q2a)).add(q1c.multiply(q2d)).subtract(q1d.multiply(q2c));
        final T y = q1a.multiply(q2c).subtract(q1b.multiply(q2d)).add(q1c.multiply(q2a)).add(q1d.multiply(q2b));
        final T z = q1a.multiply(q2d).add(q1b.multiply(q2c)).subtract(q1c.multiply(q2b)).add(q1d.multiply(q2a));

        return new FieldQuaternion<>(w, x, y, z);
    }

    /**
     * Returns the Hamilton product of the instance by a quaternion.
     *
     * @param q Quaternion.
     * @return the product of this instance with {@code q}, in that order.
     */
    public FieldQuaternion<T> multiply(final FieldQuaternion<T> q) {
        return multiply(this, q);
    }

    /**
     * Computes the sum of two quaternions.
     *
     * @param q1 Quaternion.
     * @param q2 Quaternion.
     * @return the sum of {@code q1} and {@code q2}.
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> add(final FieldQuaternion<T> q1,
                                                                             final FieldQuaternion<T> q2) {
        return new FieldQuaternion<>(q1.getQ0().add(q2.getQ0()),
                q1.getQ1().add(q2.getQ1()),
                q1.getQ2().add(q2.getQ2()),
                q1.getQ3().add(q2.getQ3()));
    }

    /**
     * Computes the sum of the instance and another quaternion.
     *
     * @param q Quaternion.
     * @return the sum of this instance and {@code q}
     */
    public FieldQuaternion<T> add(final FieldQuaternion<T> q) {
        return add(this, q);
    }

    /**
     * Subtracts two quaternions.
     *
     * @param q1 First Quaternion.
     * @param q2 Second quaternion.
     * @return the difference between {@code q1} and {@code q2}.
     */
    public static <T extends CalculusFieldElement<T>> FieldQuaternion<T> subtract(final FieldQuaternion<T> q1,
                                                                                  final FieldQuaternion<T> q2) {
        return new FieldQuaternion<>(q1.getQ0().subtract(q2.getQ0()),
                q1.getQ1().subtract(q2.getQ1()),
                q1.getQ2().subtract(q2.getQ2()),
                q1.getQ3().subtract(q2.getQ3()));
    }

    /**
     * Subtracts a quaternion from the instance.
     *
     * @param q Quaternion.
     * @return the difference between this instance and {@code q}.
     */
    public FieldQuaternion<T> subtract(final FieldQuaternion<T> q) {
        return subtract(this, q);
    }

    /**
     * Computes the dot-product of two quaternions.
     *
     * @param q1 Quaternion.
     * @param q2 Quaternion.
     * @return the dot product of {@code q1} and {@code q2}.
     */
    public static <T extends CalculusFieldElement<T>> T dotProduct(final FieldQuaternion<T> q1,
                                                                   final FieldQuaternion<T> q2) {
        return q1.getQ0().multiply(q2.getQ0())
                .add(q1.getQ1().multiply(q2.getQ1()))
                .add(q1.getQ2().multiply(q2.getQ2()))
                .add(q1.getQ3().multiply(q2.getQ3()));
    }

    /**
     * Computes the dot-product of the instance by a quaternion.
     *
     * @param q Quaternion.
     * @return the dot product of this instance and {@code q}.
     */
    public T dotProduct(final FieldQuaternion<T> q) {
        return dotProduct(this, q);
    }

    /**
     * Computes the squared norm of the quaternion.
     * @return the squared norm
     */
    public T getNormSq() {
        return q0.multiply(q0).add(q1.multiply(q1)).add(q2.multiply(q2)).add(q3.multiply(q3));
    }

    /**
     * Computes the norm of the quaternion.
     *
     * @return the norm.
     */
    public T getNorm() {
        return getNormSq().sqrt();
    }

    /**
     * Computes the normalized quaternion (the versor of the instance).
     * The norm of the quaternion must not be zero.
     *
     * @return a normalized quaternion.
     * @throws MathIllegalArgumentException if the norm of the quaternion is zero.
     */
    public FieldQuaternion<T> normalize() {

        final T norm = getNorm();
        if (norm.getReal() < Precision.SAFE_MIN) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NORM, norm);
        }

        final T inv = norm.reciprocal();
        return new FieldQuaternion<>(q0.multiply(inv), q1.multiply(inv), q2.multiply(inv), q3.multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    @SuppressWarnings("unchecked")
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof FieldQuaternion) {
            // FIXME how to properly do this class cast?
            final FieldQuaternion<T> q = (FieldQuaternion<T>) other;
            return q0 == q.getQ0() &&
                    q1 == q.getQ1() &&
                    q2 == q.getQ2() &&
                    q3 == q.getQ3();
        }

        return false;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        // "Effective Java" (second edition, p. 47).
        int result = 17;
        for (double comp : new double[] { q0.getReal(), q1.getReal(), q2.getReal(), q3.getReal() }) {
            final int c = MathUtils.hash(comp);
            result = 31 * result + c;
        }
        return result;
    }

    /**
     * Checks whether this instance is equal to another quaternion
     * within a given tolerance.
     *
     * @param q Quaternion with which to compare the current quaternion.
     * @param eps Tolerance.
     * @return {@code true} if each of the components are equal
     * within the allowed absolute error.
     */
    public boolean equals(final FieldQuaternion<T> q,
                          final double eps) {
        // FIXME CalculusFieldElements should implement equals(Object, eps)
        return Precision.equals(q0.getReal(), q.getQ0().getReal(), eps) &&
                Precision.equals(q1.getReal(), q.getQ1().getReal(), eps) &&
                Precision.equals(q2.getReal(), q.getQ2().getReal(), eps) &&
                Precision.equals(q3.getReal(), q.getQ3().getReal(), eps);
    }

    /**
     * Checks whether the instance is a unit quaternion within a given
     * tolerance.
     *
     * @param eps Tolerance (absolute error).
     * @return {@code true} if the norm is 1 within the given tolerance,
     * {@code false} otherwise
     */
    public boolean isUnitQuaternion(double eps) {
        return Precision.equals(getNorm().getReal(), 1d, eps);
    }

    /**
     * Checks whether the instance is a pure quaternion within a given
     * tolerance.
     *
     * @param eps Tolerance (absolute error).
     * @return {@code true} if the scalar part of the quaternion is zero.
     */
    public boolean isPureQuaternion(double eps) {
        return FastMath.abs(getQ0().getReal()) <= eps;
    }

    /**
     * Returns the polar form of the quaternion.
     *
     * @return the unit quaternion with positive scalar part.
     */
    public FieldQuaternion<T> getPositivePolarForm() {
        if (getQ0().getReal() < 0) {
            final FieldQuaternion<T> unitQ = normalize();
            // The quaternion of rotation (normalized quaternion) q and -q
            // are equivalent (i.e. represent the same rotation).
            return new FieldQuaternion<>(unitQ.getQ0().negate(),
                    unitQ.getQ1().negate(),
                    unitQ.getQ2().negate(),
                    unitQ.getQ3().negate());
        } else {
            return this.normalize();
        }
    }

    /**
     * Returns the inverse of this instance.
     * The norm of the quaternion must not be zero.
     *
     * @return the inverse.
     * @throws MathIllegalArgumentException if the norm (squared) of the quaternion is zero.
     */
    public FieldQuaternion<T> getInverse() {

        final T squareNorm = getNormSq();
        if (squareNorm.getReal() < Precision.SAFE_MIN) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NORM, squareNorm);
        }

        final T inv = squareNorm.reciprocal();
        return new FieldQuaternion<>(q0.multiply(inv),
                q1.multiply(inv).negate(),
                q2.multiply(inv).negate(),
                q3.multiply(inv).negate());
    }

    /**
     * Gets the first component of the quaternion (scalar part).
     *
     * @return the scalar part.
     */
    public T getQ0() {
        return q0;
    }

    /**
     * Gets the second component of the quaternion (first component
     * of the vector part).
     *
     * @return the first component of the vector part.
     */
    public T getQ1() {
        return q1;
    }

    /**
     * Gets the third component of the quaternion (second component
     * of the vector part).
     *
     * @return the second component of the vector part.
     */
    public T getQ2() {
        return q2;
    }

    /**
     * Gets the fourth component of the quaternion (third component
     * of the vector part).
     *
     * @return the third component of the vector part.
     */
    public T getQ3() {
        return q3;
    }

    /**
     * Gets the scalar part of the quaternion.
     *
     * @return the scalar part.
     * @see #getQ0()
     */
    public T getScalarPart() {
        return getQ0();
    }

    /**
     * Gets the three components of the vector part of the quaternion.
     *
     * @return the vector part.
     * @see #getQ1()
     * @see #getQ2()
     * @see #getQ3()
     */
    public T[] getVectorPart() {
        final T[] v = MathArrays.buildArray(q1.getField(), 3);
        v[0] = q1;
        v[1] = q2;
        v[2] = q3;
        return v;
    }

    /**
     * Returns a quaternion built from the real part of the instance components.
     * @return quaternion
     */
    public Quaternion toQuaternion() {
        return new Quaternion(q0.getReal(), q1.getReal(), q2.getReal(), q3.getReal());
    }

    /**
     * Multiplies the instance by a scalar.
     *
     * @param alpha Scalar factor.
     * @return a scaled quaternion.
     */
    public FieldQuaternion<T> multiply(final double alpha) {
        return new FieldQuaternion<>(q0.multiply(alpha), q1.multiply(alpha), q2.multiply(alpha), q3.multiply(alpha));
    }

    /**
     * Multiplies the instance by a scalar.
     * @param alpha scalar factor
     * @return scaled quaternion
     */
    public FieldQuaternion<T> multiply(final T alpha) {
        return new FieldQuaternion<>(q0.multiply(alpha), q1.multiply(alpha), q2.multiply(alpha), q3.multiply(alpha));
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        // TODO this might need to be adapted
        final String sp = " ";
        return "[" +
                q0 + sp +
                q1 + sp +
                q2 + sp +
                q3 +
                ']';
    }

}

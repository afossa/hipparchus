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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import java.util.Random;

public class FieldQuaternionTest {

    /** Epsilon for double comparison. */
    private static final double EPS = Math.ulp(1d);

    /** Epsilon for double comparison. */
    private static final double COMPARISON_EPS = 1e-14;

    @Test
    public final void testAccessors1() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(q0, q1, q2, q3));

        Assert.assertEquals(q0, q.getQ0().getReal(), 0);
        Assert.assertEquals(q1, q.getQ1().getReal(), 0);
        Assert.assertEquals(q2, q.getQ2().getReal(), 0);
        Assert.assertEquals(q3, q.getQ3().getReal(), 0);
    }

    @Test
    public final void testAccessors2() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(q0, q1, q2, q3));

        final Binary64 sP = q.getScalarPart();
        final Binary64[] vP = q.getVectorPart();

        Assert.assertEquals(q0, sP.getReal(), 0);
        Assert.assertEquals(q1, vP[0].getReal(), 0);
        Assert.assertEquals(q2, vP[1].getReal(), 0);
        Assert.assertEquals(q3, vP[2].getReal(), 0);
    }

    @Test
    public final void testAccessors3() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;

        final Binary64[] v = {new Binary64(q1), new Binary64(q2), new Binary64(q3)};
        final FieldQuaternion<Binary64> q = new FieldQuaternion<>(new Binary64(q0), v);

        final Binary64 sP = q.getScalarPart();
        final Binary64[] vP = q.getVectorPart();

        Assert.assertEquals(q0, sP.getReal(), 0);
        Assert.assertEquals(q1, vP[0].getReal(), 0);
        Assert.assertEquals(q2, vP[1].getReal(), 0);
        Assert.assertEquals(q3, vP[2].getReal(), 0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testWrongDimension() {
        new FieldQuaternion<>(new Binary64[]{new Binary64(1), new Binary64(2)});
    }

    @Test
    public void testGetIdentity() {
        final FieldQuaternion<Binary64> id = FieldQuaternion.getIdentity(Binary64Field.getInstance());
        Assert.assertEquals(1.0, id.getQ0().getReal(), 0);
        Assert.assertEquals(0.0, id.getQ1().getReal(), 0);
        Assert.assertEquals(0.0, id.getQ2().getReal(), 0);
        Assert.assertEquals(0.0, id.getQ3().getReal(), 0);
    }

    @Test
    public void testGetZero() {
        final FieldQuaternion<Binary64> zero = FieldQuaternion.getZero(Binary64Field.getInstance());
        Assert.assertEquals(0.0, zero.getQ0().getReal(), 0);
        Assert.assertEquals(0.0, zero.getQ1().getReal(), 0);
        Assert.assertEquals(0.0, zero.getQ2().getReal(), 0);
        Assert.assertEquals(0.0, zero.getQ3().getReal(), 0);
    }

    @Test
    public void testGetI() {
        final FieldQuaternion<Binary64> qi = FieldQuaternion.getI(Binary64Field.getInstance());
        Assert.assertEquals(0.0, qi.getQ0().getReal(), 0);
        Assert.assertEquals(1.0, qi.getQ1().getReal(), 0);
        Assert.assertEquals(0.0, qi.getQ2().getReal(), 0);
        Assert.assertEquals(0.0, qi.getQ3().getReal(), 0);
    }

    @Test
    public void testGetJ() {
        final FieldQuaternion<Binary64> qj = FieldQuaternion.getJ(Binary64Field.getInstance());
        Assert.assertEquals(0.0, qj.getQ0().getReal(), 0);
        Assert.assertEquals(0.0, qj.getQ1().getReal(), 0);
        Assert.assertEquals(1.0, qj.getQ2().getReal(), 0);
        Assert.assertEquals(0.0, qj.getQ3().getReal(), 0);
    }

    @Test
    public void testGetK() {
        final FieldQuaternion<Binary64> qk = FieldQuaternion.getK(Binary64Field.getInstance());
        Assert.assertEquals(0.0, qk.getQ0().getReal(), 0);
        Assert.assertEquals(0.0, qk.getQ1().getReal(), 0);
        Assert.assertEquals(0.0, qk.getQ2().getReal(), 0);
        Assert.assertEquals(1.0, qk.getQ3().getReal(), 0);
    }

    @Test
    public final void testConjugate() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(q0, q1, q2, q3));

        final FieldQuaternion<Binary64> qConjugate = q.getConjugate();

        Assert.assertEquals(q0, qConjugate.getQ0().getReal(), 0);
        Assert.assertEquals(-q1, qConjugate.getQ1().getReal(), 0);
        Assert.assertEquals(-q2, qConjugate.getQ2().getReal(), 0);
        Assert.assertEquals(-q3, qConjugate.getQ3().getReal(), 0);
    }

    @Test
    public final void testProductQuaternionQuaternion() {

        // Case : analytic test case
        final FieldQuaternion<Binary64> qA =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1, 0.5, -3, 4));
        final FieldQuaternion<Binary64> qB =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(6, 2, 1, -9));
        final FieldQuaternion<Binary64> qResult = FieldQuaternion.multiply(qA, qB);

        Assert.assertEquals(44, qResult.getQ0().getReal(), EPS);
        Assert.assertEquals(28, qResult.getQ1().getReal(), EPS);
        Assert.assertEquals(-4.5, qResult.getQ2().getReal(), EPS);
        Assert.assertEquals(21.5, qResult.getQ3().getReal(), EPS);

        // Conjugate of the product of two quaternions and product of their conjugates :
        // Conj(qA * qB) = Conj(qB) * Conj(qA)

        final FieldQuaternion<Binary64> conjugateOfProduct = qB.getConjugate().multiply(qA.getConjugate());
        final FieldQuaternion<Binary64> productOfConjugate = (qA.multiply(qB)).getConjugate();

        Assert.assertEquals(conjugateOfProduct.getQ0().getReal(), productOfConjugate.getQ0().getReal(), EPS);
        Assert.assertEquals(conjugateOfProduct.getQ1().getReal(), productOfConjugate.getQ1().getReal(), EPS);
        Assert.assertEquals(conjugateOfProduct.getQ2().getReal(), productOfConjugate.getQ2().getReal(), EPS);
        Assert.assertEquals(conjugateOfProduct.getQ3().getReal(), productOfConjugate.getQ3().getReal(), EPS);
    }

    @Test
    public final void testProductQuaternionVector() {

        // Case : Product between a vector and a quaternion : QxV
        final FieldQuaternion<Binary64> quaternion =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(4, 7, -1, 2));

        final Binary64[] vector = {new Binary64(2.0), new Binary64(1.0), new Binary64(3.0)};
        final FieldQuaternion<Binary64> qResultQxV =
                FieldQuaternion.multiply(quaternion, new FieldQuaternion<>(vector));

        Assert.assertEquals(-19, qResultQxV.getQ0().getReal(), EPS);
        Assert.assertEquals(3, qResultQxV.getQ1().getReal(), EPS);
        Assert.assertEquals(-13, qResultQxV.getQ2().getReal(), EPS);
        Assert.assertEquals(21, qResultQxV.getQ3().getReal(), EPS);

        // comparison with the result given by the formula :
        // qResult = (- vectorQ . vector) + (scalarQ * vector + vectorQ ^ vector)

        final Binary64[] vectorQ = quaternion.getVectorPart();

        final Binary64 scalarPartRefQxV = quaternion.getQ0().linearCombination(vectorQ, vector).negate();
        Assert.assertEquals(scalarPartRefQxV.getReal(), qResultQxV.getScalarPart().getReal(), EPS);

        // Case : Product between a vector and a quaternion : VxQ

        final FieldQuaternion<Binary64> qResultVxQ = FieldQuaternion.multiply(new FieldQuaternion<>(vector), quaternion);

        Assert.assertEquals(-19, qResultVxQ.getQ0().getReal(), EPS);
        Assert.assertEquals(13, qResultVxQ.getQ1().getReal(), EPS);
        Assert.assertEquals(21, qResultVxQ.getQ2().getReal(), EPS);
        Assert.assertEquals(3, qResultVxQ.getQ3().getReal(), EPS);

        // comparison with the result given by the formula :
        // qResult = (- vector . vectorQ) + (scalarQ * vector + vector ^ vectorQ)

        final Binary64 scalarPartRefVxQ = quaternion.getQ0().linearCombination(vectorQ, vector).negate();
        Assert.assertEquals(scalarPartRefVxQ.getReal(), qResultVxQ.getScalarPart().getReal(), EPS);
    }

    @Test
    public final void testDotProductQuaternionQuaternion() {
        // expected output
        final double expected = -6.;
        // inputs
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1, 2, 2, 1));
        final FieldQuaternion<Binary64> q2 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(3, -2, -1, -3));

        final Binary64 actual1 = FieldQuaternion.dotProduct(q1, q2);
        final Binary64 actual2 = q1.dotProduct(q2);

        Assert.assertEquals(expected, actual1.getReal(), EPS);
        Assert.assertEquals(expected, actual2.getReal(), EPS);
    }

    @Test
    public final void testScalarMultiplyDouble() {
        // expected outputs
        final double w = 1.6;
        final double x = -4.8;
        final double y = 11.20;
        final double z = 2.56;
        // inputs
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(0.5, -1.5, 3.5, 0.8));
        final double a = 3.2;

        final FieldQuaternion<Binary64> q = q1.multiply(a);

        Assert.assertEquals(w, q.getQ0().getReal(), COMPARISON_EPS);
        Assert.assertEquals(x, q.getQ1().getReal(), COMPARISON_EPS);
        Assert.assertEquals(y, q.getQ2().getReal(), COMPARISON_EPS);
        Assert.assertEquals(z, q.getQ3().getReal(), COMPARISON_EPS);
    }

    @Test
    public final void testScalarMultiplyField() {

        // expected outputs
        final double w = 1.6;
        final double x = -4.8;
        final double y = 11.20;
        final double z = 2.56;

        // inputs
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(0.5, -1.5, 3.5, 0.8));
        final Binary64 a = new Binary64(3.2);

        final FieldQuaternion<Binary64> q = q1.multiply(a);

        Assert.assertEquals(w, q.getQ0().getReal(), COMPARISON_EPS);
        Assert.assertEquals(x, q.getQ1().getReal(), COMPARISON_EPS);
        Assert.assertEquals(y, q.getQ2().getReal(), COMPARISON_EPS);
        Assert.assertEquals(z, q.getQ3().getReal(), COMPARISON_EPS);
    }

    @Test
    public final void testAddQuaternionQuaternion() {
        // expected outputs
        final double w = 4;
        final double x = -1;
        final double y = 2;
        final double z = -4;
        // inputs
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1., 2., -2., -1.));
        final FieldQuaternion<Binary64> q2 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(3., -3., 4., -3.));

        final FieldQuaternion<Binary64> qa = FieldQuaternion.add(q1, q2);
        final FieldQuaternion<Binary64> qb = q1.add(q2);

        Assert.assertEquals(w, qa.getQ0().getReal(), EPS);
        Assert.assertEquals(x, qa.getQ1().getReal(), EPS);
        Assert.assertEquals(y, qa.getQ2().getReal(), EPS);
        Assert.assertEquals(z, qa.getQ3().getReal(), EPS);

        Assert.assertEquals(w, qb.getQ0().getReal(), EPS);
        Assert.assertEquals(x, qb.getQ1().getReal(), EPS);
        Assert.assertEquals(y, qb.getQ2().getReal(), EPS);
        Assert.assertEquals(z, qb.getQ3().getReal(), EPS);
    }

    @Test
    public final void testSubtractQuaternionQuaternion() {
        // expected outputs
        final double w = -2.;
        final double x = 5.;
        final double y = -6.;
        final double z = 2.;
        // inputs
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1., 2., -2., -1.));
        final FieldQuaternion<Binary64> q2 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(3., -3., 4., -3.));

        final FieldQuaternion<Binary64> qa = FieldQuaternion.subtract(q1, q2);
        final FieldQuaternion<Binary64> qb = q1.subtract(q2);

        Assert.assertEquals(w, qa.getQ0().getReal(), EPS);
        Assert.assertEquals(x, qa.getQ1().getReal(), EPS);
        Assert.assertEquals(y, qa.getQ2().getReal(), EPS);
        Assert.assertEquals(z, qa.getQ3().getReal(), EPS);

        Assert.assertEquals(w, qb.getQ0().getReal(), EPS);
        Assert.assertEquals(x, qb.getQ1().getReal(), EPS);
        Assert.assertEquals(y, qb.getQ2().getReal(), EPS);
        Assert.assertEquals(z, qb.getQ3().getReal(), EPS);
    }

    @Test
    public final void testNorm() {

        final double q0 = 2;
        final double q1 = 1;
        final double q2 = -4;
        final double q3 = 3;
        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(q0, q1, q2, q3));
        final Binary64 norm = q.getNorm();

        Assert.assertEquals(FastMath.sqrt(30), norm.getReal(), 0);

        final Binary64 normSquareRef = FieldQuaternion.multiply(q, q.getConjugate()).getScalarPart();
        Assert.assertEquals(FastMath.sqrt(normSquareRef).getReal(), norm.getReal(), 0);
    }

    @Test
    public final void testNormalize() {

        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(2, 1, -4, -2));

        final FieldQuaternion<Binary64> versor = q.normalize();

        Assert.assertEquals(2.0 / 5.0, versor.getQ0().getReal(), 0);
        Assert.assertEquals(1.0 / 5.0, versor.getQ1().getReal(), 0);
        Assert.assertEquals(-4.0 / 5.0, versor.getQ2().getReal(), 0);
        Assert.assertEquals(-2.0 / 5.0, versor.getQ3().getReal(), 0);

        Assert.assertEquals(1, versor.getNorm().getReal(), 0);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public final void testNormalizeFail() {
        final FieldQuaternion<Binary64> zeroQ =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(0, 0, 0, 0));
        zeroQ.normalize();
    }

    @Test
    public final void testObjectEquals() {
        final Binary64 one = Binary64.ONE;

        final FieldQuaternion<Binary64> q1 = new FieldQuaternion<>(one, one, one, one);
        Assert.assertEquals(q1, q1);

        final FieldQuaternion<Binary64> q2 = new FieldQuaternion<>(one, one, one, one);
        Assert.assertEquals(q2, q1);

        final FieldQuaternion<Binary64> q3 = new FieldQuaternion<>(one, one.add(FastMath.nextUp(1.0)), one, one);
        Assert.assertNotEquals(q3, q1);
    }

    @Test
    public final void testQuaternionEquals() {
        final double inc = 1e-5;
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(2, 1, -4, -2));
        final FieldQuaternion<Binary64> q2 =
                new FieldQuaternion<>(q1.getQ0().add(inc), q1.getQ1(), q1.getQ2(), q1.getQ3());
        final FieldQuaternion<Binary64> q3 =
                new FieldQuaternion<>(q1.getQ0(), q1.getQ1().add(inc), q1.getQ2(), q1.getQ3());
        final FieldQuaternion<Binary64> q4 =
                new FieldQuaternion<>(q1.getQ0(), q1.getQ1(), q1.getQ2().add(inc), q1.getQ3());
        final FieldQuaternion<Binary64> q5 =
                new FieldQuaternion<>(q1.getQ0(), q1.getQ1(), q1.getQ2(), q1.getQ3().add(inc));

        Assert.assertFalse(q1.equals(q2, 0.9 * inc));
        Assert.assertFalse(q1.equals(q3, 0.9 * inc));
        Assert.assertFalse(q1.equals(q4, 0.9 * inc));
        Assert.assertFalse(q1.equals(q5, 0.9 * inc));

        Assert.assertTrue(q1.equals(q2, 1.1 * inc));
        Assert.assertTrue(q1.equals(q3, 1.1 * inc));
        Assert.assertTrue(q1.equals(q4, 1.1 * inc));
        Assert.assertTrue(q1.equals(q5, 1.1 * inc));
    }

    @Test
    public final void testQuaternionEquals2() {
        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1, 4, 2, 3));
        final double gap = 1e-5;
        final FieldQuaternion<Binary64> q2 = new FieldQuaternion<>(
                q1.getQ0().add(gap), q1.getQ1().add(gap), q1.getQ2().add(gap), q1.getQ3().add(gap));

        Assert.assertTrue(q1.equals(q2, 10 * gap));
        Assert.assertFalse(q1.equals(q2, gap));
        Assert.assertFalse(q1.equals(q2, gap / 10));
    }

    @Test
    public final void testIsUnitQuaternion() {
        final Random r = new Random(48);
        final int numberOfTrials = 1000;
        for (int i = 0; i < numberOfTrials; i++) {
            final FieldQuaternion<Binary64> q1 = new FieldQuaternion<>(Binary64Field.getInstance(),
                    new Quaternion(r.nextDouble(), r.nextDouble(), r.nextDouble(), r.nextDouble()));
            final FieldQuaternion<Binary64> q2 = q1.normalize();
            Assert.assertTrue(q2.isUnitQuaternion(COMPARISON_EPS));
        }

        final FieldQuaternion<Binary64> q = new FieldQuaternion<>(Binary64Field.getInstance(),
                new Quaternion(1, 1, 1, 1));
        Assert.assertFalse(q.isUnitQuaternion(COMPARISON_EPS));
    }

    @Test
    public final void testIsPureQuaternion() {
        final FieldQuaternion<Binary64> q1 = new FieldQuaternion<>(Binary64Field.getInstance(),
                new Quaternion(0, 5, 4, 8));
        Assert.assertTrue(q1.isPureQuaternion(EPS));

        final FieldQuaternion<Binary64> q2 = new FieldQuaternion<>(Binary64Field.getInstance(),
                new Quaternion(0 - EPS, 5, 4, 8));
        Assert.assertTrue(q2.isPureQuaternion(EPS));

        final FieldQuaternion<Binary64> q3 = new FieldQuaternion<>(Binary64Field.getInstance(),
                new Quaternion(0 - 1.1 * EPS, 5, 4, 8));
        Assert.assertFalse(q3.isPureQuaternion(EPS));

        final Random r = new Random(48);
        final Binary64[] v = {new Binary64(r.nextDouble()), new Binary64(r.nextDouble()), new Binary64(r.nextDouble())};
        final FieldQuaternion<Binary64> q4 = new FieldQuaternion<>(v);
        Assert.assertTrue(q4.isPureQuaternion(0));

        final FieldQuaternion<Binary64> q5 = new FieldQuaternion<>(Binary64.ZERO, v);
        Assert.assertTrue(q5.isPureQuaternion(0));
    }

    @Test
    public final void testGetInverse() {
        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1.5, 4, 2, -2.5));

        final FieldQuaternion<Binary64> inverseQ = q.getInverse();
        Assert.assertEquals(1.5 / 28.5, inverseQ.getQ0().getReal(), 0);
        Assert.assertEquals(-4.0 / 28.5, inverseQ.getQ1().getReal(), 0);
        Assert.assertEquals(-2.0 / 28.5, inverseQ.getQ2().getReal(), 0);
        Assert.assertEquals(2.5 / 28.5, inverseQ.getQ3().getReal(), 0);

        final FieldQuaternion<Binary64> product = FieldQuaternion.multiply(inverseQ, q);
        Assert.assertEquals(1, product.getQ0().getReal(), EPS);
        Assert.assertEquals(0, product.getQ1().getReal(), EPS);
        Assert.assertEquals(0, product.getQ2().getReal(), EPS);
        Assert.assertEquals(0, product.getQ3().getReal(), EPS);

        final FieldQuaternion<Binary64> qNul = FieldQuaternion.getZero(Binary64Field.getInstance());
        try {
            final FieldQuaternion<Binary64> inverseQNul = qNul.getInverse();
            Assert.fail("expecting MathIllegalArgumentException but got : " + inverseQNul);
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testGetPositivePolarForm() {

        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1.5, 4, 2, -2.5));
        final FieldQuaternion<Binary64> p1 = q1.getPositivePolarForm();
        final FieldQuaternion<Binary64> e1 = q1.normalize();

        Assert.assertEquals(e1.getQ0().getReal(), p1.getQ0().getReal(), 0);
        Assert.assertEquals(e1.getQ1().getReal(), p1.getQ1().getReal(), 0);
        Assert.assertEquals(e1.getQ2().getReal(), p1.getQ2().getReal(), 0);
        Assert.assertEquals(e1.getQ3().getReal(), p1.getQ3().getReal(), 0);

        final FieldQuaternion<Binary64> q2 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(-1.5, 4, 2, -2.5));
        final FieldQuaternion<Binary64> p2 = q2.getPositivePolarForm();
        final FieldQuaternion<Binary64> e2 = q2.normalize().multiply(-1.0);

        Assert.assertEquals(e2.getQ0().getReal(), p2.getQ0().getReal(), 0);
        Assert.assertEquals(e2.getQ1().getReal(), p2.getQ1().getReal(), 0);
        Assert.assertEquals(e2.getQ2().getReal(), p2.getQ2().getReal(), 0);
        Assert.assertEquals(e2.getQ3().getReal(), p2.getQ3().getReal(), 0);
    }

    @Test
    public void testToQuaternion() {

        final FieldQuaternion<Binary64> q1 =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1.5, 4, 2, -2.5));
        final Quaternion q2 = q1.toQuaternion();

        Assert.assertEquals(q1.getQ0().getReal(), q2.getQ0(), 0);
        Assert.assertEquals(q1.getQ1().getReal(), q2.getQ1(), 0);
        Assert.assertEquals(q1.getQ2().getReal(), q2.getQ2(), 0);
        Assert.assertEquals(q1.getQ3().getReal(), q2.getQ3(), 0);
    }

    @Test
    public final void testToString() {
        final FieldQuaternion<Binary64> q =
                new FieldQuaternion<>(Binary64Field.getInstance(), new Quaternion(1, 2, 3, 4));
        Assert.assertEquals("[1.0 2.0 3.0 4.0]", q.toString());
    }
}

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
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

package org.hipparchus.ode.nonstiff;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.interpolators.RungeKuttaFieldStateInterpolator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/**
 * This class implements the common part of all embedded Runge-Kutta
 * integrators for Ordinary Differential Equations.
 *
 * <p>These methods are embedded explicit Runge-Kutta methods with two
 * sets of coefficients allowing to estimate the error, their Butcher
 * arrays are as follows :</p>
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 *       |  b'1  b'2 ...   b's-1 b's
 * </pre>
 *
 * <p>In fact, we rather use the array defined by ej = bj - b'j to
 * compute directly the error rather than computing two estimates and
 * then comparing them.</p>
 *
 * <p>Some methods are qualified as <i>fsal</i> (first same as last)
 * methods. This means the last evaluation of the derivatives in one
 * step is the same as the first in the next step. Then, this
 * evaluation can be reused from one step to the next one and the cost
 * of such a method is really s-1 evaluations despite the method still
 * has s stages. This behaviour is true only for successful steps, if
 * the step is rejected after the error estimation phase, no
 * evaluation is saved. For an <i>fsal</i> method, we have cs = 1 and
 * asi = bi for all i.</p>
 *
 * @param <T> the type of the field elements
 */

public abstract class EmbeddedRungeKuttaFieldIntegrator<T extends CalculusFieldElement<T>>
    extends AdaptiveStepsizeFieldIntegrator<T>
    implements FieldExplicitRungeKuttaIntegrator<T> {

    /** Index of the pre-computed derivative for <i>fsal</i> methods. */
    private final int fsal;

    /** Time steps from Butcher array (without the first zero). */
    private final T[] c;

    /** Internal weights from Butcher array (without the first empty row). */
    private final T[][] a;

    /** External weights for the high order method from Butcher array. */
    private final T[] b;

    /** Time steps from Butcher array (without the first zero). */
    private double[] realC = new double[0];

    /** Internal weights from Butcher array (without the first empty row). Real version, optional. */
    private double[][] realA = new double[0][];

    /** External weights for the high order method from Butcher array. Real version, optional. */
    private double[] realB = new double[0];

    /** Stepsize control exponent. */
    private final double exp;

    /** Safety factor for stepsize control. */
    private T safety;

    /** Minimal reduction factor for stepsize control. */
    private T minReduction;

    /** Maximal growth factor for stepsize control. */
    private T maxGrowth;

    /** Flag setting whether coefficients in Butcher array are interpreted as Field or real numbers. */
    private boolean usingFieldCoefficients;

    /** Build a Runge-Kutta integrator with the given Butcher array.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param fsal index of the pre-computed derivative for <i>fsal</i> methods
     * or -1 if method is not <i>fsal</i>
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     */
    protected EmbeddedRungeKuttaFieldIntegrator(final Field<T> field, final String name, final int fsal,
                                                final double minStep, final double maxStep,
                                                final double scalAbsoluteTolerance,
                                                final double scalRelativeTolerance) {

        super(field, name, minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);

        this.fsal = fsal;
        this.c    = getC();
        this.a    = getA();
        this.b    = getB();

        exp = -1.0 / getOrder();

        // set the default values of the algorithm control parameters
        setSafety(field.getZero().add(0.9));
        setMinReduction(field.getZero().add(0.2));
        setMaxGrowth(field.getZero().add(10.0));

    }

    /** Build a Runge-Kutta integrator with the given Butcher array.
     * @param field field to which the time and state vector elements belong
     * @param name name of the method
     * @param fsal index of the pre-computed derivative for <i>fsal</i> methods
     * or -1 if method is not <i>fsal</i>
     * @param minStep minimal step (must be positive even for backward
     * integration), the last step can be smaller than this
     * @param maxStep maximal step (must be positive even for backward
     * integration)
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     */
    protected EmbeddedRungeKuttaFieldIntegrator(final Field<T> field, final String name, final int fsal,
                                                final double   minStep, final double maxStep,
                                                final double[] vecAbsoluteTolerance,
                                                final double[] vecRelativeTolerance) {

        super(field, name, minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
        this.usingFieldCoefficients = false;

        this.fsal = fsal;
        this.c    = getC();
        this.a    = getA();
        this.b    = getB();

        exp = -1.0 / getOrder();

        // set the default values of the algorithm control parameters
        setSafety(field.getZero().add(0.9));
        setMinReduction(field.getZero().add(0.2));
        setMaxGrowth(field.getZero().add(10.0));

    }

    /** Create an interpolator.
     * @param forward integration direction indicator
     * @param yDotK slopes at the intermediate points
     * @param globalPreviousState start of the global step
     * @param globalCurrentState end of the global step
     * @param mapper equations mapper for the all equations
     * @return external weights for the high order method from Butcher array
     */
    protected abstract RungeKuttaFieldStateInterpolator<T> createInterpolator(boolean forward, T[][] yDotK,
                                                                              FieldODEStateAndDerivative<T> globalPreviousState,
                                                                              FieldODEStateAndDerivative<T> globalCurrentState,
                                                                              FieldEquationsMapper<T> mapper);

    /** Get the order of the method.
     * @return order of the method
     */
    public abstract int getOrder();

    /** Get the safety factor for stepsize control.
     * @return safety factor
     */
    public T getSafety() {
        return safety;
    }

    /** Set the safety factor for stepsize control.
     * @param safety safety factor
     */
    public void setSafety(final T safety) {
        this.safety = safety;
    }

    /**
     * Setter for the flag between real or Field coefficients in the Butcher array.
     *
     * @param usingFieldCoefficients new value for flag
     */
    public void setUsingFieldCoefficients(boolean usingFieldCoefficients) {
        this.usingFieldCoefficients = usingFieldCoefficients;
    }

    /** {@inheritDoc} */
    @Override
    public boolean isUsingFieldCoefficients() {
        return usingFieldCoefficients;
    }

    /** {@inheritDoc} */
    @Override
    public int getNumberOfStages() {
        return b.length;
    }

    /** {@inheritDoc} */
    @Override
    protected FieldODEStateAndDerivative<T> initIntegration(FieldExpandableODE<T> eqn, FieldODEState<T> s0, T t) {
        if (!isUsingFieldCoefficients()) {
            realA = getRealA();
            realB = getRealB();
            realC = getRealC();
        }
        return super.initIntegration(eqn, s0, t);
    }

    /** {@inheritDoc} */
    @Override
    public FieldODEStateAndDerivative<T> integrate(final FieldExpandableODE<T> equations,
                                                   final FieldODEState<T> initialState, final T finalTime)
        throws MathIllegalArgumentException, MathIllegalStateException {

        sanityChecks(initialState, finalTime);
        setStepStart(initIntegration(equations, initialState, finalTime));
        final boolean forward = finalTime.subtract(initialState.getTime()).getReal() > 0;

        // create some internal working arrays
        final int   stages = getNumberOfStages();
        final T[][] yDotK  = MathArrays.buildArray(getField(), stages, -1);
        T[]   yTmp   = MathArrays.buildArray(getField(), equations.getMapper().getTotalDimension());

        // set up integration control objects
        T  hNew           = getField().getZero();
        boolean firstTime = true;

        // main integration loop
        setIsLastStep(false);
        do {

            // iterate over step size, ensuring local normalized error is smaller than 1
            double error = 10.0;
            while (error >= 1.0) {

                // first stage
                final T[] y = getStepStart().getCompleteState();
                yDotK[0] = getStepStart().getCompleteDerivative();

                if (firstTime) {
                    final StepsizeHelper helper = getStepSizeHelper();
                    final T[] scale = MathArrays.buildArray(getField(), helper.getMainSetDimension());
                    for (int i = 0; i < scale.length; ++i) {
                        scale[i] = helper.getTolerance(i, y[i].abs());
                    }
                    hNew = getField().getZero().add(initializeStep(forward, getOrder(), scale, getStepStart()));
                    firstTime = false;
                }

                setStepSize(hNew);
                if (forward) {
                    if (getStepStart().getTime().add(getStepSize()).subtract(finalTime).getReal() >= 0) {
                        setStepSize(finalTime.subtract(getStepStart().getTime()));
                    }
                } else {
                    if (getStepStart().getTime().add(getStepSize()).subtract(finalTime).getReal() <= 0) {
                        setStepSize(finalTime.subtract(getStepStart().getTime()));
                    }
                }

                // next stages
                if (isUsingFieldCoefficients()) {
                    FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(getEquations(),
                            getStepStart().getTime(), y, getStepSize(), a, c, yDotK);
                    yTmp = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y, yDotK, getStepSize(), b);
                } else {
                    FieldExplicitRungeKuttaIntegrator.applyInternalButcherWeights(getEquations(),
                            getStepStart().getTime(), y, getStepSize(), realA, realC, yDotK);
                    yTmp = FieldExplicitRungeKuttaIntegrator.applyExternalButcherWeights(y, yDotK, getStepSize(), realB);
                }

                incrementEvaluations(stages - 1);

                // estimate the error at the end of the step
                error = estimateError(yDotK, y, yTmp, getStepSize());
                if (error >= 1.0) {
                    // reject the step and attempt to reduce error by stepsize control
                    final T factor = MathUtils.min(maxGrowth,
                                                   MathUtils.max(minReduction, safety.multiply(FastMath.pow(error, exp))));
                    hNew = getStepSizeHelper().filterStep(getStepSize().multiply(factor), forward, false);
                }

            }
            final T   stepEnd = getStepStart().getTime().add(getStepSize());
            final T[] yDotTmp = (fsal >= 0) ? yDotK[fsal] : computeDerivatives(stepEnd, yTmp);
            final FieldODEStateAndDerivative<T> stateTmp = equations.getMapper().mapStateAndDerivative(stepEnd, yTmp, yDotTmp);

            // local error is small enough: accept the step, trigger events and step handlers
            setStepStart(acceptStep(createInterpolator(forward, yDotK, getStepStart(), stateTmp, equations.getMapper()),
                                    finalTime));

            if (!isLastStep()) {

                // stepsize control for next step
                final T factor = MathUtils.min(maxGrowth,
                                               MathUtils.max(minReduction, safety.multiply(FastMath.pow(error, exp))));
                final T  scaledH    = getStepSize().multiply(factor);
                final T  nextT      = getStepStart().getTime().add(scaledH);
                final boolean nextIsLast = forward ?
                                           nextT.subtract(finalTime).getReal() >= 0 :
                                           nextT.subtract(finalTime).getReal() <= 0;
                hNew = getStepSizeHelper().filterStep(scaledH, forward, nextIsLast);

                final T  filteredNextT      = getStepStart().getTime().add(hNew);
                final boolean filteredNextIsLast = forward ?
                                                   filteredNextT.subtract(finalTime).getReal() >= 0 :
                                                   filteredNextT.subtract(finalTime).getReal() <= 0;
                if (filteredNextIsLast) {
                    hNew = finalTime.subtract(getStepStart().getTime());
                }

            }

        } while (!isLastStep());

        final FieldODEStateAndDerivative<T> finalState = getStepStart();
        resetInternalState();
        return finalState;

    }

    /** Get the minimal reduction factor for stepsize control.
     * @return minimal reduction factor
     */
    public T getMinReduction() {
        return minReduction;
    }

    /** Set the minimal reduction factor for stepsize control.
     * @param minReduction minimal reduction factor
     */
    public void setMinReduction(final T minReduction) {
        this.minReduction = minReduction;
    }

    /** Get the maximal growth factor for stepsize control.
     * @return maximal growth factor
     */
    public T getMaxGrowth() {
        return maxGrowth;
    }

    /** Set the maximal growth factor for stepsize control.
     * @param maxGrowth maximal growth factor
     */
    public void setMaxGrowth(final T maxGrowth) {
        this.maxGrowth = maxGrowth;
    }

    /** Compute the error ratio.
     * @param yDotK derivatives computed during the first stages
     * @param y0 estimate of the step at the start of the step
     * @param y1 estimate of the step at the end of the step
     * @param h  current step
     * @return error ratio, greater than 1 if step should be rejected
     */
    protected abstract double estimateError(T[][] yDotK, T[] y0, T[] y1, T h);

}

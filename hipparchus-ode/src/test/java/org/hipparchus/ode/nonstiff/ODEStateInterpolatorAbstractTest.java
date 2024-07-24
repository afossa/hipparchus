/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.ode.nonstiff;


import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.sampling.AbstractODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public abstract class ODEStateInterpolatorAbstractTest {

    @Test
    public abstract void interpolationAtBounds();

    protected void doInterpolationAtBounds(double epsilon) {
        ODEStateInterpolator interpolator = setUpInterpolator(new SinCos(), 0.0, new double[] { 0.0, 1.0 }, 0.125);

        Assertions.assertEquals(0.0, interpolator.getPreviousState().getTime(), 1.0e-15);
        for (int i = 0; i < 2; ++i) {
            Assertions.assertEquals(interpolator.getPreviousState().getPrimaryState()[i],
                                interpolator.getInterpolatedState(interpolator.getPreviousState().getTime()).getPrimaryState()[i],
                                epsilon);
        }
        Assertions.assertEquals(0.125, interpolator.getCurrentState().getTime(), 1.0e-15);
        for (int i = 0; i < 2; ++i) {
            Assertions.assertEquals(interpolator.getCurrentState().getPrimaryState()[i],
                                interpolator.getInterpolatedState(interpolator.getCurrentState().getTime()).getPrimaryState()[i],
                                epsilon);
        }
        Assertions.assertFalse(interpolator.isPreviousStateInterpolated());
        Assertions.assertFalse(interpolator.isCurrentStateInterpolated());
    }

    @Test
    public abstract void interpolationInside();

    protected void doInterpolationInside(double epsilonSin, double epsilonCos) {

        SinCos sinCos = new SinCos();
        ODEStateInterpolator interpolator = setUpInterpolator(sinCos, 0.0, new double[] { 0.0, 1.0 }, 0.0125);

        int n = 100;
        double maxErrorSin = 0;
        double maxErrorCos = 0;
        for (int i = 0; i <= n; ++i) {
            double t =     ((n - i) * interpolator.getPreviousState().getTime() +
                                 i  * interpolator.getCurrentState().getTime()) / n;
            final double[] interpolated = interpolator.getInterpolatedState(t).getPrimaryState();
            double[] reference = sinCos.theoreticalState(t);
            maxErrorSin = FastMath.max(maxErrorSin, FastMath.abs(interpolated[0] - reference[0]));
            maxErrorCos = FastMath.max(maxErrorCos, FastMath.abs(interpolated[1] - reference[1]));
        }

        Assertions.assertEquals(0.0, maxErrorSin, epsilonSin);
        Assertions.assertEquals(0.0, maxErrorCos, epsilonCos);

        Assertions.assertFalse(interpolator.isPreviousStateInterpolated());
        Assertions.assertFalse(interpolator.isCurrentStateInterpolated());
    }

    @Test
    public abstract void restrictPrevious();

    protected void doRestrictPrevious(double epsilon, double epsilonDot) {

        AbstractODEStateInterpolator original   = setUpInterpolator(new SinCos(), 0.0, new double[] { 0.0, 1.0 }, 0.125);

        Assertions.assertFalse(original.isPreviousStateInterpolated());
        Assertions.assertFalse(original.isCurrentStateInterpolated());

        AbstractODEStateInterpolator restricted = original.restrictStep(original.getInterpolatedState(1.0 / 32),
                                                                        original.getCurrentState());

        Assertions.assertSame(original.getPreviousState(),       original.getGlobalPreviousState());
        Assertions.assertSame(original.getCurrentState(),        original.getGlobalCurrentState());
        Assertions.assertSame(original.getGlobalPreviousState(), restricted.getGlobalPreviousState());
        Assertions.assertSame(original.getGlobalCurrentState(),  restricted.getGlobalCurrentState());
        Assertions.assertNotSame(restricted.getPreviousState(),  restricted.getGlobalPreviousState());
        Assertions.assertSame(restricted.getCurrentState(),      restricted.getGlobalCurrentState());
        Assertions.assertEquals(1.0 / 32, restricted.getPreviousState().getTime(), 1.0e-15);
        Assertions.assertTrue(restricted.isPreviousStateInterpolated());
        Assertions.assertFalse(restricted.isCurrentStateInterpolated());

        checkRestricted(original, restricted, epsilon, epsilonDot);

    }

    @Test
    public abstract void restrictCurrent();

    protected void doRestrictCurrent(double epsilon, double epsilonDot) {

        AbstractODEStateInterpolator original   = setUpInterpolator(new SinCos(), 0.0, new double[] { 0.0, 1.0 }, 0.125);

        Assertions.assertFalse(original.isPreviousStateInterpolated());
        Assertions.assertFalse(original.isCurrentStateInterpolated());

        AbstractODEStateInterpolator restricted = original.restrictStep(original.getPreviousState(),
                                                                        original.getInterpolatedState(3.0 / 32));

        Assertions.assertSame(original.getPreviousState(),       original.getGlobalPreviousState());
        Assertions.assertSame(original.getCurrentState(),        original.getGlobalCurrentState());
        Assertions.assertSame(original.getGlobalPreviousState(), restricted.getGlobalPreviousState());
        Assertions.assertSame(original.getGlobalCurrentState(),  restricted.getGlobalCurrentState());
        Assertions.assertSame(restricted.getPreviousState(),     restricted.getGlobalPreviousState());
        Assertions.assertNotSame(restricted.getCurrentState(),   restricted.getGlobalCurrentState());
        Assertions.assertEquals(3.0 / 32, restricted.getCurrentState().getTime(), 1.0e-15);
        Assertions.assertFalse(restricted.isPreviousStateInterpolated());
        Assertions.assertTrue(restricted.isCurrentStateInterpolated());

        checkRestricted(original, restricted, epsilon, epsilonDot);

    }

    @Test
    public abstract void restrictBothEnds();

    protected void doRestrictBothEnds(double epsilon, double epsilonDot) {

        AbstractODEStateInterpolator original   = setUpInterpolator(new SinCos(), 0.0, new double[] { 0.0, 1.0 }, 0.125);

        Assertions.assertFalse(original.isPreviousStateInterpolated());
        Assertions.assertFalse(original.isCurrentStateInterpolated());

        AbstractODEStateInterpolator restricted = original.restrictStep(original.getInterpolatedState(1.0 / 32),
                                                                        original.getInterpolatedState(3.0 / 32));

        Assertions.assertSame(original.getPreviousState(),       original.getGlobalPreviousState());
        Assertions.assertSame(original.getCurrentState(),        original.getGlobalCurrentState());
        Assertions.assertSame(original.getGlobalPreviousState(), restricted.getGlobalPreviousState());
        Assertions.assertSame(original.getGlobalCurrentState(),  restricted.getGlobalCurrentState());
        Assertions.assertNotSame(restricted.getPreviousState(),  restricted.getGlobalPreviousState());
        Assertions.assertNotSame(restricted.getCurrentState(),   restricted.getGlobalCurrentState());
        Assertions.assertEquals(1.0 / 32, restricted.getPreviousState().getTime(), 1.0e-15);
        Assertions.assertEquals(3.0 / 32, restricted.getCurrentState().getTime(), 1.0e-15);
        Assertions.assertTrue(restricted.isPreviousStateInterpolated());
        Assertions.assertTrue(restricted.isCurrentStateInterpolated());

        checkRestricted(original, restricted, epsilon, epsilonDot);

    }

    @Test
    public abstract void degenerateInterpolation();

    protected void doDegenerateInterpolation() {
        AbstractODEStateInterpolator interpolator = setUpInterpolator(new SinCos(), 0.0, new double[] { 0.0, 1.0 }, 0.0);
        ODEStateAndDerivative interpolatedState = interpolator.getInterpolatedState(0.0);
        Assertions.assertEquals(0.0, interpolatedState.getTime(), 0.0);
        Assertions.assertEquals(0.0, interpolatedState.getPrimaryState()[0], 0.0);
        Assertions.assertEquals(1.0, interpolatedState.getPrimaryState()[1], 0.0);
        Assertions.assertEquals(1.0, interpolatedState.getPrimaryDerivative()[0], 0.0);
        Assertions.assertEquals(0.0, interpolatedState.getPrimaryDerivative()[1], 0.0);
    }

    private void checkRestricted(AbstractODEStateInterpolator original, AbstractODEStateInterpolator restricted,
                                 double epsilon, double epsilonDot) {
        for (double t = restricted.getPreviousState().getTime();
             t <= restricted.getCurrentState().getTime();
             t += 1.0 / 256) {
            ODEStateAndDerivative originalInterpolated   = original.getInterpolatedState(t);
            ODEStateAndDerivative restrictedInterpolated = restricted.getInterpolatedState(t);
            Assertions.assertEquals(t, originalInterpolated.getTime(), 1.0e-15);
            Assertions.assertEquals(t, restrictedInterpolated.getTime(), 1.0e-15);
            Assertions.assertEquals(originalInterpolated.getPrimaryState()[0],
                                restrictedInterpolated.getPrimaryState()[0],
                                epsilon);
            Assertions.assertEquals(originalInterpolated.getPrimaryState()[1],
                                restrictedInterpolated.getPrimaryState()[1],
                                epsilon);
            Assertions.assertEquals(originalInterpolated.getPrimaryDerivative()[0],
                                restrictedInterpolated.getPrimaryDerivative()[0],
                                epsilonDot);
            Assertions.assertEquals(originalInterpolated.getPrimaryDerivative()[1],
                                restrictedInterpolated.getPrimaryDerivative()[1],
                                epsilonDot);
        }

    }

    public interface ReferenceODE extends OrdinaryDifferentialEquation {
        double[]              theoreticalState(double t);
        DerivativeStructure[] theoreticalState(DerivativeStructure t);
    }

    protected abstract AbstractODEStateInterpolator setUpInterpolator(final ReferenceODE eqn,
                                                                      final double t0, final double[] y0,
                                                                      final double t1);

    private static class SinCos implements ReferenceODE {
        public int getDimension() {
            return 2;
        }
        public double[] computeDerivatives(final double t, final double[] y) {
            return new double[] {
                y[1], -y[0]
            };
        }
        public double[] theoreticalState(final double t) {
            return new double[] {
                FastMath.sin(t), FastMath.cos(t)
            };
        }
        public DerivativeStructure[] theoreticalState(final DerivativeStructure t) {
            return new DerivativeStructure[] {
                t.sin(), t.cos()
            };
        }
    }

}

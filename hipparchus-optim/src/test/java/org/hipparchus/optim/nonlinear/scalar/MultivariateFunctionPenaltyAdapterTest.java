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
package org.hipparchus.optim.nonlinear.scalar;

import org.hipparchus.analysis.MultivariateFunction;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.MaxEval;
import org.hipparchus.optim.PointValuePair;
import org.hipparchus.optim.SimplePointChecker;
import org.hipparchus.optim.nonlinear.scalar.noderiv.AbstractSimplex;
import org.hipparchus.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.hipparchus.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class MultivariateFunctionPenaltyAdapterTest {
    @Test
    public void testStartSimplexInsideRange() {
        final BiQuadratic biQuadratic = new BiQuadratic(2.0, 2.5, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionPenaltyAdapter wrapped
              = new MultivariateFunctionPenaltyAdapter(biQuadratic,
                                                       biQuadratic.getLower(),
                                                       biQuadratic.getUpper(),
                                                       1000.0, new double[] { 100.0, 100.0 });

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[] { 1.0, 0.5 });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(300),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { 1.5, 2.25 }));

        Assertions.assertEquals(biQuadratic.getBoundedXOptimum(), optimum.getPoint()[0], 2e-7);
        Assertions.assertEquals(biQuadratic.getBoundedYOptimum(), optimum.getPoint()[1], 2e-7);
    }

    @Test
    public void testStartSimplexOutsideRange() {
        final BiQuadratic biQuadratic = new BiQuadratic(2.0, 2.5, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionPenaltyAdapter wrapped
              = new MultivariateFunctionPenaltyAdapter(biQuadratic,
                                                       biQuadratic.getLower(),
                                                       biQuadratic.getUpper(),
                                                       1000.0, new double[] { 100.0, 100.0 });

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[] { 1.0, 0.5 });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(300),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -1.5, 4.0 }));

        Assertions.assertEquals(biQuadratic.getBoundedXOptimum(), optimum.getPoint()[0], 2e-7);
        Assertions.assertEquals(biQuadratic.getBoundedYOptimum(), optimum.getPoint()[1], 2e-7);
    }

    @Test
    public void testOptimumOutsideRange() {
        final BiQuadratic biQuadratic = new BiQuadratic(4.0, 0.0, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionPenaltyAdapter wrapped
            =  new MultivariateFunctionPenaltyAdapter(biQuadratic,
                                                      biQuadratic.getLower(),
                                                      biQuadratic.getUpper(),
                                                      1000.0, new double[] { 100.0, 100.0 });

        SimplexOptimizer optimizer = new SimplexOptimizer(new SimplePointChecker<PointValuePair>(1.0e-11, 1.0e-20));
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[] { 1.0, 0.5 });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(600),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -1.5, 4.0 }));

        Assertions.assertEquals(biQuadratic.getBoundedXOptimum(), optimum.getPoint()[0], 2e-7);
        Assertions.assertEquals(biQuadratic.getBoundedYOptimum(), optimum.getPoint()[1], 2e-7);
    }

    @Test
    public void testUnbounded() {
        final BiQuadratic biQuadratic = new BiQuadratic(4.0, 0.0,
                                                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                                                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        final MultivariateFunctionPenaltyAdapter wrapped
            = new MultivariateFunctionPenaltyAdapter(biQuadratic,
                                                     biQuadratic.getLower(),
                                                     biQuadratic.getUpper(),
                                                     1000.0, new double[] { 100.0, 100.0 });

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[] { 1.0, 0.5 });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(300),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -1.5, 4.0 }));

        Assertions.assertEquals(biQuadratic.getBoundedXOptimum(), optimum.getPoint()[0], 2e-7);
        Assertions.assertEquals(biQuadratic.getBoundedYOptimum(), optimum.getPoint()[1], 2e-7);
    }

    @Test
    public void testHalfBounded() {
        final BiQuadratic biQuadratic = new BiQuadratic(4.0, 4.0,
                                                        1.0, Double.POSITIVE_INFINITY,
                                                        Double.NEGATIVE_INFINITY, 3.0);
        final MultivariateFunctionPenaltyAdapter wrapped
              = new MultivariateFunctionPenaltyAdapter(biQuadratic,
                                                       biQuadratic.getLower(),
                                                       biQuadratic.getUpper(),
                                                       1000.0, new double[] { 100.0, 100.0 });

        SimplexOptimizer optimizer = new SimplexOptimizer(new SimplePointChecker<PointValuePair>(1.0e-10, 1.0e-20));
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[] { 1.0, 0.5 });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(400),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -1.5, 4.0 }));

        Assertions.assertEquals(biQuadratic.getBoundedXOptimum(), optimum.getPoint()[0], 2e-7);
        Assertions.assertEquals(biQuadratic.getBoundedYOptimum(), optimum.getPoint()[1], 2e-7);
    }

    private static class BiQuadratic implements MultivariateFunction {

        private final double xOptimum;
        private final double yOptimum;

        private final double xMin;
        private final double xMax;
        private final double yMin;
        private final double yMax;

        public BiQuadratic(final double xOptimum, final double yOptimum,
                           final double xMin, final double xMax,
                           final double yMin, final double yMax) {
            this.xOptimum = xOptimum;
            this.yOptimum = yOptimum;
            this.xMin     = xMin;
            this.xMax     = xMax;
            this.yMin     = yMin;
            this.yMax     = yMax;
        }

        public double value(double[] point) {
            // the function should never be called with out of range points
            Assertions.assertTrue(point[0] >= xMin);
            Assertions.assertTrue(point[0] <= xMax);
            Assertions.assertTrue(point[1] >= yMin);
            Assertions.assertTrue(point[1] <= yMax);

            final double dx = point[0] - xOptimum;
            final double dy = point[1] - yOptimum;
            return dx * dx + dy * dy;

        }

        public double[] getLower() {
            return new double[] { xMin, yMin };
        }

        public double[] getUpper() {
            return new double[] { xMax, yMax };
        }

        public double getBoundedXOptimum() {
            return (xOptimum < xMin) ? xMin : ((xOptimum > xMax) ? xMax : xOptimum);
        }

        public double getBoundedYOptimum() {
            return (yOptimum < yMin) ? yMin : ((yOptimum > yMax) ? yMax : yOptimum);
        }

    }

}

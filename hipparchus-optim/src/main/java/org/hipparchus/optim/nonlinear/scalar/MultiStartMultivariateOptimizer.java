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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.optim.BaseMultiStartMultivariateOptimizer;
import org.hipparchus.optim.PointValuePair;
import org.hipparchus.random.RandomVectorGenerator;

/**
 * Multi-start optimizer.
 * This class wraps an optimizer in order to use it several times in
 * turn with different starting points (trying to avoid being trapped
 * in a local extremum when looking for a global one).
 *
 */
public class MultiStartMultivariateOptimizer
    extends BaseMultiStartMultivariateOptimizer<PointValuePair> {
    /** Underlying optimizer. */
    private final MultivariateOptimizer optimizer;
    /** Found optima. */
    private final List<PointValuePair> optima;

    /**
     * Create a multi-start optimizer from a single-start optimizer.
     *
     * @param optimizer Single-start optimizer to wrap.
     * @param starts Number of starts to perform.
     * If {@code starts == 1}, the result will be same as if {@code optimizer}
     * is called directly.
     * @param generator Random vector generator to use for restarts.
     * @throws NullArgumentException if {@code optimizer} or {@code generator}
     * is {@code null}.
     * @throws MathIllegalArgumentException if {@code starts < 1}.
     */
    public MultiStartMultivariateOptimizer(final MultivariateOptimizer optimizer,
                                           final int starts,
                                           final RandomVectorGenerator generator)
        throws MathIllegalArgumentException, NullArgumentException {
        super(optimizer, starts, generator);
        this.optimizer = optimizer;
        this.optima   = new ArrayList<>();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public PointValuePair[] getOptima() {
        optima.sort(getPairComparator());
        return optima.toArray(new PointValuePair[0]);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void store(PointValuePair optimum) {
        optima.add(optimum);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void clear() {
        optima.clear();
    }

    /**
     * @return a comparator for sorting the optima.
     */
    private Comparator<PointValuePair> getPairComparator() {
        return new Comparator<PointValuePair>() {
            /** {@inheritDoc} */
            @Override
            public int compare(final PointValuePair o1,
                               final PointValuePair o2) {
                if (o1 == null) {
                    return (o2 == null) ? 0 : 1;
                } else if (o2 == null) {
                    return -1;
                }
                final double v1 = o1.getValue();
                final double v2 = o2.getValue();
                return (optimizer.getGoalType() == GoalType.MINIMIZE) ?
                    Double.compare(v1, v2) : Double.compare(v2, v1);
            }
        };
    }
}

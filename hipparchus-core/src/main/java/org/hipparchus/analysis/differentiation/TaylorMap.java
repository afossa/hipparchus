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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.MathUtils;

/** Container for a Taylor map.
 * <p>
 * A Taylor map is a set of n {@link DerivativeStructure}
 * \((f_1, f_2, \ldots, f_n)\) depending on m parameters \((p_1, p_2, \ldots, p_m)\),
 * with positive n and m.
 * </p>
 * @since 2.2
 */
public class TaylorMap {

    /** Evaluation point. */
    private final double[] point;

    /** Mapping functions. */
    private final DerivativeStructure[] functions;

    /** Simple constructor.
     * <p>
     * The number of number of parameters and derivation orders of all
     * functions must match.
     * </p>
     * @param point point at which map is evaluated
     * @param functions functions composing the map (must contain at least one element)
     */
    public TaylorMap(final double[] point, final DerivativeStructure[] functions) {
        if (point == null || point.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE,
                                                   point == null ? 0 : point.length);
        }
        if (functions == null || functions.length == 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE,
                                                   functions == null ? 0 : functions.length);
        }
        this.point     = point.clone();
        this.functions = functions.clone();
        final DSFactory factory0 = functions[0].getFactory();
        MathUtils.checkDimension(point.length, factory0.getCompiler().getFreeParameters());
        for (int i = 1; i < functions.length; ++i) {
            factory0.checkCompatibility(functions[i].getFactory());
        }
    }

    /** Constructor for identity map.
     * <p>
     * The identity is considered to be evaluated at origin.
     * </p>
     * @param parameters number of free parameters
     * @param order derivation order
     * @param nbFunctions number of functions
     */
    public TaylorMap(final int parameters, final int order, final int nbFunctions) {
        this(parameters, nbFunctions);
        final DSFactory factory = new DSFactory(parameters, order);
        for (int i = 0; i < nbFunctions; ++i) {
            functions[i] = factory.variable(i, 0.0);
        }
    }

    /** Build an empty map evaluated at origin.
     * @param parameters number of free parameters
     * @param nbFunctions number of functions
     */
    private TaylorMap(final int parameters, final int nbFunctions) {
        this.point     = new double[parameters];
        this.functions = new DerivativeStructure[nbFunctions];
    }

    /** Get the number of parameters of the map.
     * @return number of parameters of the map
     */
    public int getNbParameters() {
        return point.length;
    }

    /** Get the number of functions of the map.
     * @return number of functions of the map
     */
    public int getNbFunctions() {
        return functions.length;
    }

    /** Get the point at which map is evaluated.
     * @return point at which map is evaluated
     */
    public double[] getPoint() {
        return point.clone();
    }

    /** Get a function from the map.
     * @param i index of the function (must be between 0 included and {@link #getNbFunctions()} excluded
     * @return function at index i
     */
    public DerivativeStructure getFunction(final int i) {
        return functions[i];
    }

    /** Subtract two maps.
     * @param map map to subtract from instance
     * @return this - map
     */
    private TaylorMap subtract(final TaylorMap map) {
        final TaylorMap result = new TaylorMap(point.length, functions.length);
        for (int i = 0; i < result.functions.length; ++i) {
            result.functions[i] = functions[i].subtract(map.functions[i]);
        }
        return result;
    }

    /** Evaluate Taylor expansion of the map at some offset.
     * @param deltaP parameters offsets \((\Delta p_1, \Delta p_2, \ldots, \Delta p_n)\)
     * @return value of the Taylor expansion at \((p_1 + \Delta p_1, p_2 + \Delta p_2, \ldots, p_n + \Delta p_n)\)
     */
    public double[] value(final double... deltaP) {
        final double[] value = new double[functions.length];
        for (int i = 0; i < functions.length; ++i) {
            value[i] = functions[i].taylor(deltaP);
        }
        return value;
    }

    /** Compose the instance with another Taylor map as \(\mathrm{this} \circ \mathrm{other}\).
     * @param other map with which instance must be composed
     * @return composed map \(\mathrm{this} \circ \mathrm{other}\)
     */
    public TaylorMap compose(final TaylorMap other) {

        // safety check
        MathUtils.checkDimension(getNbParameters(), other.getNbFunctions());

        final DerivativeStructure[] composed = new DerivativeStructure[functions.length];
        for (int i = 0; i < functions.length; ++i) {
            composed[i] = functions[i].rebase(other.functions);
        }

        return new TaylorMap(other.point, composed);

    }

    /** Invert the instance.
     * <p>
     * Consider {@link #value(double[]) Taylor expansion} of the map with
     * small parameters offsets \((\Delta p_1, \Delta p_2, \ldots, \Delta p_n)\)
     * which leads to evaluation offsets \((f_1 + df_1, f_2 + df_2, \ldots, f_n + df_n)\).
     * The map inversion defines a Taylor map that computes \((\Delta p_1,
     * \Delta p_2, \ldots, \Delta p_n)\) from \((df_1, df_2, \ldots, df_n)\).
     * </p>
     * <p>
     * The map must be square to be invertible (i.e. the number of functions and the
     * number of parameters in the functions must match)
     * </p>
     * @param decomposer matrix decomposer to user for inverting the linear part
     * @return inverted map
     * @see <a href="https://doi.org/10.1016/S1076-5670(08)70228-3">chapter
     * 2 of Advances in Imaging and Electron Physics, vol 108
     * by Martin Berz</a>
     */
    public TaylorMap invert(final MatrixDecomposer decomposer) {

        final DSFactory  factory  = functions[0].getFactory();
        final DSCompiler compiler = factory.getCompiler();
        final int        n        = functions.length;

        // safety check
        MathUtils.checkDimension(n, functions[0].getFreeParameters());

        // set up an indirection array between linear terms and complete derivatives arrays
        final int[] indirection    = new int[n];
        int linearIndex = 0;
        for (int k = 1; linearIndex < n; ++k) {
            if (compiler.getPartialDerivativeOrdersSum(k) == 1) {
                indirection[linearIndex++] = k;
            }
        }

        // separate linear and non-linear terms
        final RealMatrix linear      = MatrixUtils.createRealMatrix(n, n);
        final TaylorMap  nonLinearTM = new TaylorMap(n, n);
        for (int i = 0; i < n; ++i) {
            nonLinearTM.functions[i] = factory.build(functions[i].getAllDerivatives());
            nonLinearTM.functions[i].setDerivativeComponent(0, 0.0);
            for (int j = 0; j < n; ++j) {
                final int k = indirection[j];
                linear.setEntry(i, j, functions[i].getDerivativeComponent(k));
                nonLinearTM.functions[i].setDerivativeComponent(k, 0.0);
            }
        }

        // invert the linear part
        final RealMatrix linearInvert = decomposer.decompose(linear).getInverse();

        // convert the invert of linear part back to a Taylor map
        final TaylorMap  linearInvertTM = new TaylorMap(n, n);
        for (int i = 0; i < n; ++i) {
            linearInvertTM.functions[i] = new DerivativeStructure(factory);
            for (int j = 0; j < n; ++j) {
                linearInvertTM.functions[i].setDerivativeComponent(indirection[j], linearInvert.getEntry(i, j));
            }
        }

        // perform fixed-point evaluation of the inverse
        // adding one derivation order at each iteration
        final TaylorMap identity = new TaylorMap(n, compiler.getOrder(), n);
        TaylorMap invertTM = linearInvertTM;
        for (int k = 1; k < compiler.getOrder(); ++k) {
            invertTM = linearInvertTM.compose(identity.subtract(nonLinearTM.compose(invertTM)));
        }

        // set the constants
        for (int i = 0; i < n; ++i) {
            invertTM.point[i] = functions[i].getValue();
            invertTM.functions[i].setDerivativeComponent(0, point[i]);
        }

        return invertTM;

    }

}

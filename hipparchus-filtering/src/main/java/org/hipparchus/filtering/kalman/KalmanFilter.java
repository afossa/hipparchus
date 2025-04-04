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

package org.hipparchus.filtering.kalman;

import org.hipparchus.exception.MathRuntimeException;

/**
 * Interface representing a Kalman filter.
 * @param <T> the type of the measurements
 * @since 1.3
 */
public interface KalmanFilter<T extends Measurement> extends KalmanEstimate {

    /** Perform one estimation step.
     * @param measurement single measurement to handle
     * @return estimated state after measurement has been considered
     * @exception MathRuntimeException if estimation fails
     */
    ProcessEstimate estimationStep(T measurement)
        throws MathRuntimeException;

    /** Set the filter observer callback.
     * @param observer the observer
     */
    void setObserver(KalmanObserver observer);
}

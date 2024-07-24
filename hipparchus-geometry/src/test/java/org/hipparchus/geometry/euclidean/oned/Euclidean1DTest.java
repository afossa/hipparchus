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
package org.hipparchus.geometry.euclidean.oned;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.geometry.Space;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

public class Euclidean1DTest {

    @Test
    public void testDimension() {
        Assertions.assertEquals(1, Euclidean1D.getInstance().getDimension());
    }

    @Test
    public void testSubSpace() {
        assertThrows(Euclidean1D.NoSubSpaceException.class, () -> {
            Euclidean1D.getInstance().getSubSpace();
        });
    }

    @Test
    public void testSerialization() {
        Space e1 = Euclidean1D.getInstance();
        Space deserialized = (Space) UnitTestUtils.serializeAndRecover(e1);
        Assertions.assertTrue(e1 == deserialized);
    }

}

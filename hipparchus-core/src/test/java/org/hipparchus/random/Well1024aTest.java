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
package org.hipparchus.random;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class Well1024aTest extends RandomGeneratorAbstractTest {

    @Override
    protected RandomGenerator makeGenerator() {
        return new Well1024a(1001);
    }

    @Test
    public void testReferenceCode() {
        Well1024a mt = new Well1024a(new int[] {
            740849862,  1202665156,  -199039369,  -259008301,  -291878969, -1164428990, -1565918811,   491009864,
          -1883086670,  1383450241,  1244617256,   689006653, -1576746370, -1307940314,  1421489086,  1742094000,
           -595495729,  1047766204,  1875773301, -1637793284,  1379017098,   262792705,   191880010,  -251000180,
          -1753047622,  -972355720,    90626881,  1644693418,  1503365577,   439653419,  1806361562,  1268823869
       });
        int[] refInt = {
           -1478749726,  -1645579484,  -2075363835,  -2063444174,  -1834148336,  -1769045872,    -40711346,   1717441026,
            2130656771,    783441285,    570433609,   1560023451,    653233971,   1368672434,    -72036215,   1071111800,
             933776492,     26114960,     49888778,   1808107515,   1092989296,    754848337,   1336994364,  -1987450448,
            -691190146,  -1803870839,   1110716866,   1173269113,   -391000050,   2014216908,    180756301,   -382891013,
           -1908154585,   1580737629,   1080267957,   -125532248,   2094530239,   2132964485,   -438596348,   -760299445,
            1058181869,   2050816800,  -1534429037,    -62552782,    824524142,   -818590371,  -1857695907,   -684762866,
            -156556543,   -902759995,   -880795194,  -1387351132,  -1263017515,    448006597,    201038266,   1929826313,
            -455367306,    672963027,   2000073013,  -1546842042,    446341090,   1001696686,   -779919012,   -347722602,
           -1342821677,   1639571150,   -835315755,   1505585376,    367004975,  -2035864404,  -1786623553,   1249724913,
             182435312,   1444514513,   1815333708,   1333772382,    299664001,   -284691169,   2034403374,   1423310887,
           -1319051884,   1557286441,   -445198266,   -251809030,   1602786123,    944036382,  -1020529634,    258344235,
             685254367,   1838964943,   -156674528,   -979736602,   -538312836,    234643178,    211152102,   -635498640,
           -1036733933,  -1347589147,   -565609042,  -1358714165,    508618483,   -786364693,   2071450261,   1206956772,
            -678931458,    167690617,    144698821,   1719720781,   1575869280,  -1343221123,  -1766469944,    284991647,
            -717305514,    892653651,  -1368347075,   -615701972,   -730369849,   1360396003,  -1869287623,   1778269052,
            -586061545,   -699517114,     61530249,  -1860611767,   -519660852,   1841085925,   1555610093,   -399979337,
            -790345742,    422355947,   2007965433,   2044952550,  -1712164595,   -102915702,   -693865324,  -1894042487,
           -1285020072,   -215883074,     95833252,   1625818040,  -1055951680,    513067085,   1825246558,   -553461652,
           -1923361799,  -1869480206,    567232636,  -1751727150,  -1832301399,   -108136455,  -1312244126,     14006795,
             850221366,   -382389732,  -1741556188,  -1317464467,   1948314870,    753994471,   1028235947,    342494132,
           -1862256693,    723808794,   -234257642,   1609928369,   -802733456,   1315831915,   1436072885,   1224767136,
            2144557791,  -1839965886,    224821018,  -1461697757,  -1080386760,   1638573498,  -1188173812,   -325181523,
           -1750676219,  -1780415850,    698793362,   -908352052,    299746482,   -161660934,   1938166833,    800297005,
              56640033,  -1214932666,  -1248124842,   1822796868,   1777615881,   -718517774,   1908159957,   1733053281,
            1851844331,   1283519375,  -1771494956,   2060179999,   1666129209,   1919453531,   -498145770,    697567008,
            1855487148,  -1587163491,    565216434,  -1477877933,   -925662919,   -806492585,  -1201439047,  -1424534232,
            1788616523,     69414717,    655893636,  -1175978556,     24787512,   -861550001,    439525754,   -190433174,
            -383811606,   -508589783,   1441608687,    608181366,   1539467064,    925903122,    697209654,   1878283393,
           -1967567432,  -1659677763,   -249658183,    847096354,    397741956,   -125334541,  -1286840731,   1016461908,
            -997968592,   1795331475,   1856856501,  -1716726445,   -582181331,   -887091847,    426964855,   -609219941,
           -1456232632,   -483467616,   1069260754,    972242064,  -1406786247,   1954194029,     52627891,   1212755081,
            2117436668,    281073392,    741537353,   -483063506,   1850906286,   -244876135,   -270818140,   1817568823
        };

        for (int i = 0; i < refInt.length; ++i) {
            Assertions.assertEquals(refInt[i], mt.nextInt());
        }

    }

}

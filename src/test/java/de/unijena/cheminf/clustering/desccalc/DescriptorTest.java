/*
 * CDK-Descriptor-Calculation
 * Copyright (C) 2025 Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * Source code is available at <https://github.com/JonasSchaub/CDK-Descriptor-Calculation>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.clustering.desccalc;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test class for descriptor
 */
class DescriptorTest {

    /**
     * Test method for descriptor MOLECULAR_WEIGHT
     */
    @Test
    public void test_MOLECULAR_WEIGHT() {
        // Ethanol
        String tmpSmiles = "CC(=O)O";
        String[] tmpSmilesArray = new String[] {tmpSmiles};
        float[][] tmpMatrix = new float[][]
            {
                {0f}
            };
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT};
        int tmpNumberOfConcurrentCalculationThreads = 0;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(Descriptor.setCalculatedDescriptorComponents(tmpDescriptors, tmpSmilesArray, tmpMatrix, tmpStartIndex, tmpNumberOfConcurrentCalculationThreads));
            Assertions.assertEquals(500.0f, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test method for descriptor WIENER_INDEX
     */
    @Test
    public void test_WIENER_INDEX() {
        // Ethanol
        String tmpSmiles = "CC(=O)O";
        String[] tmpSmilesArray = new String[] {tmpSmiles};
        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.WIENER_NUMBER};
        int tmpNumberOfConcurrentCalculationThreads = 0;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(Descriptor.setCalculatedDescriptorComponents(tmpDescriptors, tmpSmilesArray, tmpMatrix, tmpStartIndex, tmpNumberOfConcurrentCalculationThreads));
            Assertions.assertEquals(70.0f, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test multiple descriptors
     */
    @Test
    public void test_MultipleDescriptors() {
        // Ethanol
        String tmpSmiles = "CC(=O)O";
        String[] tmpSmilesArray = new String[] {tmpSmiles};
        int tmpStartIndex = 0;
        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f}
                };
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT, Descriptor.WIENER_NUMBER};
        int tmpNumberOfConcurrentCalculationThreads = 0;

        try {
            Assertions.assertEquals(2, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(Descriptor.setCalculatedDescriptorComponents(tmpDescriptors, tmpSmilesArray, tmpMatrix, tmpStartIndex, tmpNumberOfConcurrentCalculationThreads));
            Assertions.assertEquals(500.0f, tmpMatrix[0][0]);
            Assertions.assertEquals(70.0f, tmpMatrix[0][1]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

}

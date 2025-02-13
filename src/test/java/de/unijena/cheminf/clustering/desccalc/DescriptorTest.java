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
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

/**
 * Test class for Descriptor class.
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 */
class DescriptorTest {

    /**
     * Test method for descriptor MOLECULAR_WEIGHT
     */
    @Test
    public void test_MOLECULAR_WEIGHT() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        float[][] tmpMatrix = new float[][]
            {
                {0f}
            };
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT};
        int tmpNumberOfConcurrentCalculationThreads = 0;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(
                Descriptor.setCalculatedDescriptorComponents(
                    tmpDescriptors,
                    tmpMoleculesArray,
                    tmpMatrix,
                    tmpStartIndex,
                    tmpNumberOfConcurrentCalculationThreads
                )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test method for descriptor WIENER_INDEX
     */
    @Test
    public void test_WIENER_INDEX() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        float[][] tmpMatrix = new float[][]
                {
                    {0f, 0f}
                };
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.WIENER_NUMBER};
        int tmpNumberOfConcurrentCalculationThreads = 0;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(2, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(
                Descriptor.setCalculatedDescriptorComponents(
                    tmpDescriptors,
                    tmpMoleculesArray,
                    tmpMatrix,
                    tmpStartIndex,
                    tmpNumberOfConcurrentCalculationThreads
                )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test multiple descriptors
     */
    @Test
    public void test_MultipleDescriptors() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT, Descriptor.WIENER_NUMBER};
        int tmpNumberOfConcurrentCalculationThreads = 0;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(
                Descriptor.setCalculatedDescriptorComponents(
                    tmpDescriptors,
                    tmpMoleculesArray,
                    tmpMatrix,
                    tmpStartIndex,
                    tmpNumberOfConcurrentCalculationThreads
                )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[0][2]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test multiple molecules with multiple descriptors: Sequential and parallel
     */
    @Test
    public void test_MultipleMoleculesAndDescriptors_SequentialParallel() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        int tmpNumberOfMolecules = 10;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        float[][] tmpMatrix = new float[tmpNumberOfMolecules][];
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            tmpMoleculesArray[i] = tmpMolecule;
            tmpMatrix[i] = new float[] {0f, 0f, 0f};
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT, Descriptor.WIENER_NUMBER};
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

        try {
            // Sequential code
            int tmpNumberOfConcurrentCalculationThreads = 0;
            Assertions.assertTrue(
                Descriptor.setCalculatedDescriptorComponents(
                    tmpDescriptors,
                    tmpMoleculesArray,
                    tmpMatrix,
                    tmpStartIndex,
                    tmpNumberOfConcurrentCalculationThreads
                )
            );
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[i][0]));
                Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[i][1]));
                Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[i][2]));
            }
        } catch (Exception anException) {
            Assertions.fail();
        }

        // Reset descriptor matrix
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            tmpMoleculesArray[i] = tmpMolecule;
            tmpMatrix[i] = new float[] {0f, 0f, 0f};
        }

        try {
            // Parallel code
            int tmpNumberOfConcurrentCalculationThreads = 8;
            Assertions.assertTrue(
                Descriptor.setCalculatedDescriptorComponents(
                    tmpDescriptors,
                    tmpMoleculesArray,
                    tmpMatrix,
                    tmpStartIndex,
                    tmpNumberOfConcurrentCalculationThreads
                )
            );
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[i][0]));
                Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[i][1]));
                Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[i][2]));
            }
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

}

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
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Locale;

/**
 * Test class for Descriptor class.
 * Note: For adding tests of a new descriptor goto "Add new descriptor tests here!"
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 */
class DescriptorTest {

    //<editor-fold desc="Single descriptor tests">
    /**
     * Tests method for descriptor MOLECULAR_WEIGHT
     */
    @Test
    public void test_MOLECULAR_WEIGHT() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor WIENER_NUMBER
     */
    @Test
    public void test_WIENER_NUMBER() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.WIENER_NUMBER};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(2, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart

            tmpMatrix = new float[][]
                    {
                            {0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart

            tmpMatrix = new float[][]
                    {
                            {0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor ATOM_COUNT
     */
    @Test
    public void test_ATOM_COUNT() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ATOM_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };

            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor H_BOND_ACCEPTOR_COUNT_1
     */
    @Test
    public void test_H_BOND_ACCEPTOR_COUNT_1() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_ACCEPTOR_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 2 hydrogen bond acceptors

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 2 hydrogen bond acceptors

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 2 hydrogen bond acceptors
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor H_BOND_ACCEPTOR_COUNT_2
     */
    @Test
    public void test_H_BOND_ACCEPTOR_COUNT_2() throws Exception {
        // Nitro-indazole derivative
        String tmpSmiles = "O=N(=O)c1cccc2cn[nH]c12";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_ACCEPTOR_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Nitro-indazole derivative has 1 hydrogen bond acceptor

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Nitro-indazole derivative has 1 hydrogen bond acceptor

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Nitro-indazole derivative has 1 hydrogen bond acceptor
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor H_BOND_DONOR_COUNT_1
     */
    @Test
    public void test_H_BOND_DONOR_COUNT_1() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_DONOR_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor H_BOND_DONOR_COUNT_2
     */
    @Test
    public void test_H_BOND_DONOR_COUNT_2() throws Exception {
        // Phenol
        String tmpSmiles = "Oc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_DONOR_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Phenol has 1 hydrogen bond donor

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Phenol has 1 hydrogen bond donor

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Phenol has 1 hydrogen bond donor
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor TPSA
     */
    @Test
    public void test_TPSA_1() throws Exception {
        String tmpSmiles = "C=NC(CC#N)N(C)C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.TPSA};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor TPSA
     */
    @Test
    public void test_TPSA_2() throws Exception {
        String tmpSmiles = "CCCN(=O)=O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.TPSA};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor LARGEST_CHAIN
     */
    @Test
    public void test_LARGEST_CHAIN_1() throws Exception {
        String tmpSmiles = "C=CC=Cc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LARGEST_CHAIN};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor LARGEST_CHAIN
     */
    @Test
    public void test_LARGEST_CHAIN_2() throws Exception {
        String tmpSmiles = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LARGEST_CHAIN};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor LONGEST_ALIPHATIC_CHAIN
     */
    @Test
    public void test_LONGEST_ALIPHATIC_CHAIN_1() throws Exception {
        String tmpSmiles = "CCCCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LONGEST_ALIPHATIC_CHAIN};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor LONGEST_ALIPHATIC_CHAIN with a more complex structure
     */
    @Test
    public void test_LONGEST_ALIPHATIC_CHAIN_2() throws Exception {
        String tmpSmiles = "CC(C)(C)c2ccc(OCCCC(=O)Nc1nccs1)cc2";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LONGEST_ALIPHATIC_CHAIN};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor MANNHOLD_LOGP
     */
    @Test
    public void test_MANNHOLD_LOGP() throws Exception {
        String tmpSmiles = "C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MANNHOLD_LOGP};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BCUT
     */
    @Test
    public void test_BCUT() throws Exception {
        String tmpSmiles = "CC(=O)N";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BCUT};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.00001; // tolerance range

        try {
            Assertions.assertEquals(6, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(11.881587, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(16.005958, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(-0.381844, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.325509, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(3.374638, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(5.033583, tmpMatrix[0][5], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(11.881587, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(16.005958, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(-0.381844, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.325509, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(3.374638, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(5.033583, tmpMatrix[0][5], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(11.881587, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(16.005958, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(-0.381844, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.325509, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(3.374638, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(5.033583, tmpMatrix[0][5], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BOND_COUNT_ALL
     */
    @Test
    public void test_BOND_COUNT_ALL_1() throws Exception {
        String tmpSmiles = "CCO";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_ALL};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BOND_COUNT_ALL
     */
    @Test
    public void test_BOND_COUNT_ALL_2() throws Exception {
        String tmpSmiles = "C=C=C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_ALL};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BOND_COUNT_ALL
     */
    @Test
    public void test_BOND_COUNT_ALL_3() throws Exception {
        String tmpSmiles = "CC#N";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_ALL};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BOND_COUNT_SPECIFIED
     */
    @Test
    public void test_BOND_COUNT_SPECIFIED_1() throws Exception {
        String tmpSmiles = "C=C=C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_SPECIFIED};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BOND_COUNT_SPECIFIED
     */
    @Test
    public void test_BOND_COUNT_SPECIFIED_2() throws Exception {
        String tmpSmiles = "CCO";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_SPECIFIED};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BOND_COUNT_SPECIFIED
     */
    @Test
    public void test_BOND_COUNT_SPECIFIED_3() throws Exception {
        String tmpSmiles = "C#N";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_SPECIFIED};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor B_POL
     */
    @Test
    public void test_B_POL() throws Exception {
        String tmpSmiles = "O=C(O)CC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.B_POL};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor RULE_OF_FIVE
     */
    @Test
    public void test_RULE_OF_FIVE() throws Exception {
        String tmpSmiles = "CCCC(OCC)OCC(c1cccc2ccccc12)C4CCC(CCCO)C(CC3CNCNC3)C4";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.RULE_OF_FIVE};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor AROMATIC_ATOMS_COUNT
     */
    @Test
    public void test_AROMATIC_ATOMS_COUNT() throws Exception {
        String tmpSmiles = "c1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AROMATIC_ATOMS_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor AROMATIC_BONDS_COUNT with benzene
     */
    @Test
    public void test_AROMATIC_BONDS_COUNT() throws Exception {
        String tmpSmiles = "c1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AROMATIC_BONDS_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor ROTATABLE_BONDS_COUNT with an amide bond
     */
    @Test
    public void test_ROTATABLE_BONDS_COUNT_amide() throws Exception {
        String tmpSmiles = "CCNC(=O)CC(C)C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ROTATABLE_BONDS_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor FMF
     */
    @Test
    public void test_FMF() throws Exception {
        String tmpSmiles = "Clc1cc(cc(Cl)c1N)C(O)CNC(C)(C)C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FMF};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor FRACTIONAL_CSP3
     */
    @Test
    public void test_FRACTIONAL_CSP3() throws Exception {
        String tmpSmiles = "CC1=CC=CC(C)=N1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FRACTIONAL_CSP3};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor HYBRIDIZATION_RATIO
     */
    @Test
    public void test_HYBRIDIZATION_RATIO() throws Exception {
        String tmpSmiles = "CCC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.HYBRIDIZATION_RATIO};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor KAPPA_SHAPE_INDICES
     */
    @Test
    public void test_KAPPA_SHAPE_INDICES() throws Exception {
        String tmpSmiles = "O=C(O)CC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.KAPPA_SHAPE_INDICES};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(5.0, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(2.25, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(4.0, tmpMatrix[0][2], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(5.0, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(2.25, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(4.0, tmpMatrix[0][2], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(5.0, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(2.25, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(4.0, tmpMatrix[0][2], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor PETITJEAN_NUMBER
     */
    @Test
    public void test_PETITJEAN_NUMBER() throws Exception {
        String tmpSmiles = "O=C(O)CC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.PETITJEAN_NUMBER};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor SPIRO_ATOM_COUNT
     */
    @Test
    public void test_SPIRO_ATOM_COUNT() throws Exception {
        String tmpSmiles = "C1CCC2(CC1)CC=C1C=CC=CC1=C2";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.SPIRO_ATOM_COUNT};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor V_ADJ_MAT
     */
    @Test
    public void test_V_ADJ_MAT() throws Exception {
        String tmpSmiles = "C1CCC2CCCCC2C1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.V_ADJ_MAT};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.001; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor WEIGHTED_PATH
     */
    @Test
    public void test_WEIGHTED_PATH() throws Exception {
        String tmpSmiles = "CCCC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.WEIGHTED_PATH};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.00001; // tolerance range

        try {
            Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(6.87132, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.71783, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][4], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(6.87132, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.71783, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][4], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(6.87132, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.71783, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0, tmpMatrix[0][4], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor ZAGREB_INDEX
     */
    @Test
    public void test_ZAGREB_INDEX() throws Exception {
        String tmpSmiles = "O=C(O)CC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ZAGREB_INDEX};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CARBON_TYPES
     */
    @Test
    public void test_CARBON_TYPES() throws Exception {
        String tmpSmiles = "CCCC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CARBON_TYPES};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(9, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][3]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][4]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][5]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][6]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][7]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][8]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][3]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][4]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][5]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][6]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][7]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][8]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][2]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][3]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][4]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][5]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][6]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][7]));
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][8]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor A_LOG_P
     */
    @Test
    public void test_A_LOG_P() throws Exception {
        String tmpSmiles = "CCCCl";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
        // Aromaticity detection and marking
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.A_LOG_P};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.000", tmpSymbols);

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.719", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2.955", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("20.584", tmpFormat.format(tmpMatrix[0][2]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.719", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2.955", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("20.584", tmpFormat.format(tmpMatrix[0][2]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("1.719", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2.955", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("20.584", tmpFormat.format(tmpMatrix[0][2]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor X_LOG_P
     */
    @Test
    public void test_X_LOG_P() throws Exception {
        String tmpSmiles = "O=C(O)C(N)CCCN";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
        // Aromaticity detection and marking
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.X_LOG_P};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.1; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor JP_LOG_P
     */
    @Test
    public void test_JP_LOG_P() throws Exception {
        String tmpSmiles = "CCC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.JP_LOG_P};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.1; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor A_POL
     */
    @Test
    public void test_A_POL() throws Exception {
        String tmpSmiles = "O=C(O)CC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.A_POL};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    // Add new descriptor tests here!

    //</editor-fold>

    //<editor-fold desc="Combined descriptor tests">
    /**
     * Tests combined descriptors
     */
    @Test
    public void test_CombinedDescriptors() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        // Add new descriptor tests here!
        Descriptor[] tmpDescriptors =
                new Descriptor[]
                        {
                                Descriptor.MOLECULAR_WEIGHT,
                                Descriptor.WIENER_NUMBER,
                                Descriptor.ATOM_COUNT,
                                Descriptor.H_BOND_ACCEPTOR_COUNT,
                                Descriptor.H_BOND_DONOR_COUNT

                        };
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        try {
            Assertions.assertEquals(6, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[1][];
            tmpMatrix[0] = new float[Descriptor.getNumberOfComponents(tmpDescriptors)];
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[0][2]));
            Assertions.assertEquals("8.00", tmpFormat.format(tmpMatrix[0][3]));
            Assertions.assertEquals("2.00", tmpFormat.format(tmpMatrix[0][4]));
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][5]));

            tmpMatrix = new float[1][];
            tmpMatrix[0] = new float[Descriptor.getNumberOfComponents(tmpDescriptors)];
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[0][2]));
            Assertions.assertEquals("8.00", tmpFormat.format(tmpMatrix[0][3]));
            Assertions.assertEquals("2.00", tmpFormat.format(tmpMatrix[0][4]));
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][5]));

            tmpMatrix = new float[1][];
            tmpMatrix[0] = new float[Descriptor.getNumberOfComponents(tmpDescriptors)];
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[0][2]));
            Assertions.assertEquals("8.00", tmpFormat.format(tmpMatrix[0][3]));
            Assertions.assertEquals("2.00", tmpFormat.format(tmpMatrix[0][4]));
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][5]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }
    //</editor-fold>

    //<editor-fold desc="Tests with all implemented descriptors">
    /**
     * Tests parallelization
     */
    @Test
    public void test_Parallelization() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);

        int tmpNumberOfMolecules = 10;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        Arrays.fill(tmpMoleculesArray, tmpMolecule);
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getAllDescriptors();
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);

        try {

            float[][] tmpMatrixSequential = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixSequential[i] = new float[tmpNumberOfComponents];
            }
            boolean tmpIsParallelCalculation = false;
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            float[][] tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                for (int j = 0; j < tmpNumberOfComponents; j++) {
                    Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
                }
            }

            tmpMatrixSequential = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixSequential[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = false;
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                for (int j = 0; j < tmpNumberOfComponents; j++) {
                    Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
                }
            }

            tmpMatrixSequential = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixSequential[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = false;
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                for (int j = 0; j < tmpNumberOfComponents; j++) {
                    Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
                }
            }

        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests integrity
     */
    @Test
    public void test_Integrity() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        // Aromaticity detection and marking
        Cycles.markRingAtomsAndBonds((tmpMolecule));
        Aromaticity.apply(Aromaticity.Model.Daylight, tmpMolecule);
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        Arrays.fill(tmpMoleculesArray, tmpMolecule);
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getAllDescriptors();
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);
        boolean tmpIsParallelCalculation = false;

        try {

            float[][] tmpMatrix1 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix1[i] = new float[tmpNumberOfComponents];
            }
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix1,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            float[][] tmpMatrix2 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix2[i] = new float[tmpNumberOfComponents];
            }
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix2,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            float[][] tmpMatrix3 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix3[i] = new float[tmpNumberOfComponents];
            }
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix3,
                            tmpStartIndex,
                            tmpIsParallelCalculation
                    )
            );

            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                for (int j = 0; j < tmpNumberOfComponents; j++) {
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix2[i][j]);
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix3[i][j]);
                }
            }

        } catch (Exception anException) {
            Assertions.fail();
        }
    }
    //</editor-fold>

}

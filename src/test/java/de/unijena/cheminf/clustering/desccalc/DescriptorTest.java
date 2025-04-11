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
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT};
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
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.WIENER_NUMBER};
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
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.ATOM_COUNT};
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
     * Acetic Acid
     */
    @Test
    public void test_H_BOND_ACCEPTOR_COUNT_1() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.H_BOND_ACCEPTOR_COUNT};
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
     * Nitro-indazole derivative
     */
    @Test
    public void test_H_BOND_ACCEPTOR_COUNT_2() throws Exception {
        // Nitro-indazole derivative
        String tmpSmiles = "O=N(=O)c1cccc2cn[nH]c12";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.H_BOND_ACCEPTOR_COUNT};
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
     * Acetic acid
     */
    @Test
    public void test_H_BOND_DONOR_COUNT_1() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.H_BOND_DONOR_COUNT};
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
     * Phenol
     */
    @Test
    public void test_H_BOND_DONOR_COUNT_2() throws Exception {
        // Phenol
        String tmpSmiles = "Oc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.H_BOND_DONOR_COUNT};
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
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        // Add new descriptor tests here!
        Descriptor[] tmpDescriptors =
            new Descriptor[]
                {
                    Descriptor.MOLECULER_WEIGHT,
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

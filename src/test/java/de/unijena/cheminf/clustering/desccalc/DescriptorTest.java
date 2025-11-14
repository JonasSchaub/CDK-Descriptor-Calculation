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
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.BitSet;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;

/**
 * Test class for Descriptor class.
 * Note: For adding tests of a new descriptor goto "Add new descriptor tests here!"
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
class DescriptorTest {
    //<editor-fold desc="Single descriptor tests">
    /**
     * Test method for descriptor MOLECULAR_WEIGHT.
     */
    @Test
    public void test_MOLECULAR_WEIGHT() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O"; //Acetic Acid CID: 176
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor WIENER_NUMBER.
     */
    @Test
    public void test_WIENER_NUMBER() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O"; //Acetic Acid CID: 176
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart

            tmpMatrix = new float[][]
                    {
                            {0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart

            tmpMatrix = new float[][]
                    {
                            {0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart

            tmpMatrix = new float[][]
                    {
                            {0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("9", tmpFormat.format(tmpMatrix[0][0])); // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
            Assertions.assertEquals("0", tmpFormat.format(tmpMatrix[0][1])); // there are no atoms that are 3 bonds apart
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor ATOM_COUNT.
     */
    @Test
    public void test_ATOM_COUNT() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O"; //Acetic Acid CID: 176
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

            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("8", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 8 atoms
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for the organic subset of individual atom counts (C, H, N, O, S, P, F, Br, Cl, I) in a complex molecule.
     */
    @Test
    public void test_ATOM_COUNT_ORGANIC_SUBSET() throws Exception {
        String tmpSmiles = "CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)OP(=O)([O-])OC(C)CNC(=O)CCC4(C(C5C6(C(C(C(=N6)C(=C7C(C(C(=N7)C=C8C(C(C(=N8)C(=C4[N-]5)C)CCC(=O)N)(C)C)CCC(=O)N)(C)CC(=O)N)C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[Co+3]"; // Cobalamin CID: 74413906
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{
                Descriptor.ATOM_COUNT_C,
                Descriptor.ATOM_COUNT_H,
                Descriptor.ATOM_COUNT_N,
                Descriptor.ATOM_COUNT_O,
                Descriptor.ATOM_COUNT_S,
                Descriptor.ATOM_COUNT_P,
                Descriptor.ATOM_COUNT_F,
                Descriptor.ATOM_COUNT_BR,
                Descriptor.ATOM_COUNT_CL,
                Descriptor.ATOM_COUNT_I
        };
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(10, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            // Expected values for the organic subset atoms in this complex molecule
            Assertions.assertEquals(62, tmpMatrix[0][0]); // C count
            Assertions.assertEquals(88, tmpMatrix[0][1]); // H count
            Assertions.assertEquals(13, tmpMatrix[0][2]); // N count
            Assertions.assertEquals(14, tmpMatrix[0][3]);  // O count
            Assertions.assertEquals(0, tmpMatrix[0][4]);  // S count
            Assertions.assertEquals(1, tmpMatrix[0][5]);  // P count
            Assertions.assertEquals(0, tmpMatrix[0][6]);  // F count
            Assertions.assertEquals(0, tmpMatrix[0][7]);  // Br count
            Assertions.assertEquals(0, tmpMatrix[0][8]);  // Cl count
            Assertions.assertEquals(0, tmpMatrix[0][9]);  // I count

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(62, tmpMatrix[0][0]); // C count
            Assertions.assertEquals(88, tmpMatrix[0][1]); // H count
            Assertions.assertEquals(13, tmpMatrix[0][2]); // N count
            Assertions.assertEquals(14, tmpMatrix[0][3]);  // O count
            Assertions.assertEquals(0, tmpMatrix[0][4]);  // S count
            Assertions.assertEquals(1, tmpMatrix[0][5]);  // P count
            Assertions.assertEquals(0, tmpMatrix[0][6]);  // F count
            Assertions.assertEquals(0, tmpMatrix[0][7]);  // Br count
            Assertions.assertEquals(0, tmpMatrix[0][8]);  // Cl count
            Assertions.assertEquals(0, tmpMatrix[0][9]);  // I count

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f,}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(62, tmpMatrix[0][0]); // C count
            Assertions.assertEquals(88, tmpMatrix[0][1]); // H count
            Assertions.assertEquals(13, tmpMatrix[0][2]); // N count
            Assertions.assertEquals(14, tmpMatrix[0][3]);  // O count
            Assertions.assertEquals(0, tmpMatrix[0][4]);  // S count
            Assertions.assertEquals(1, tmpMatrix[0][5]);  // P count
            Assertions.assertEquals(0, tmpMatrix[0][6]);  // F count
            Assertions.assertEquals(0, tmpMatrix[0][7]);  // Br count
            Assertions.assertEquals(0, tmpMatrix[0][8]);  // Cl count
            Assertions.assertEquals(0, tmpMatrix[0][9]);  // I count

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f,}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(62, tmpMatrix[0][0]); // C count
            Assertions.assertEquals(88, tmpMatrix[0][1]); // H count
            Assertions.assertEquals(13, tmpMatrix[0][2]); // N count
            Assertions.assertEquals(14, tmpMatrix[0][3]);  // O count
            Assertions.assertEquals(0, tmpMatrix[0][4]);  // S count
            Assertions.assertEquals(1, tmpMatrix[0][5]);  // P count
            Assertions.assertEquals(0, tmpMatrix[0][6]);  // F count
            Assertions.assertEquals(0, tmpMatrix[0][7]);  // Br count
            Assertions.assertEquals(0, tmpMatrix[0][8]);  // Cl count
            Assertions.assertEquals(0, tmpMatrix[0][9]);  // I count
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor H_BOND_ACCEPTOR_COUNT.
     */
    @Test
    public void test_H_BOND_ACCEPTOR_COUNT() throws Exception {
        // Acetic acid
        String tmpSmiles1 = "CC(=O)O"; //Acetic Acid CID: 176
        String tmpSmiles2 = "O=N(=O)c1cccc2cn[nH]c12"; // 7-Nitroindole CID: 1893
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ElectronDonation[] models = {
                // Mdl and PiBonds lead to a calculation of 2 rather than 1 H-bond acceptor which is correct
                // because of the different handling of aromaticity in the second molecule but this behavior
                // is not tested here
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                //Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles
                //Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model for the first molecule
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_ACCEPTOR_COUNT};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor H_BOND_DONOR_COUNT.
     */
    @Test
    public void test_H_BOND_DONOR_COUNT() throws Exception {
        String tmpSmiles1 = "CC(=O)O"; //Acetic Acid CID: 176
        String tmpSmiles2 = "Oc1ccccc1"; // Phenol CID: 996
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model for the first molecule
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_DONOR_COUNT};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0])); // Phenol has 1 hydrogen bond donor

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0])); // Phenol has 1 hydrogen bond donor

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0])); // Phenol has 1 hydrogen bond donor

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Acetic acid has 1 hydrogen bond donor
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[1][0])); // Phenol has 1 hydrogen bond donor
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor TPSA.
     */
    @Test
    public void test_TPSA() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Preparation of the first molecule
        String tmpSmiles1 = "C=NC(CC#N)N(C)C"; // Not known in the PubChem database, but a valid SMILES
        // Preparation of the second Molecule
        String tmpSmiles2 = "CCCN(=O)=O"; // 1-Nitropropane CID: 7903

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model for the first molecule
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.TPSA};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("39.39", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("45.82", tmpFormat.format(tmpMatrix[1][0]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor LARGEST_CHAIN.
     */
    @Test
    public void test_LARGEST_CHAIN() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Preparation of the first molecule
        String tmpSmiles1 = "C=CC=Cc1ccccc1"; // 1-Phenylbutadiene CID: 137048
        // Preparation of the second molecule
        String tmpSmiles2 = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2"; // Not known in the PubChem database, but a valid SMILES
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model for the first molecule
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LARGEST_CHAIN};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[1][0]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor LONGEST_ALIPHATIC_CHAIN.
     */
    @Test
    public void test_LONGEST_ALIPHATIC_CHAIN() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Preparation of the first molecule
        String tmpSmiles1 = "CCCCc1ccccc1"; // Butylbenzene CID: 7705
        // Preparation of the second molecule
        String tmpSmiles2 = "CC(C)(C)c2ccc(OCCCC(=O)Nc1nccs1)cc2"; // 4-(4-tert-butylphenoxy)-N-(1,3-thiazol-2-yl)butanamide CID: 1565007
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model for the first molecule
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LONGEST_ALIPHATIC_CHAIN};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[1][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}, {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[1][0]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor MANNHOLD_LOGP.
     */
    @Test
    public void test_MANNHOLD_LOGP() throws Exception {
        String tmpSmiles = "C"; // Methane CID: 297
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.57", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor BCUT.
     */
    @Test
    public void test_BCUT() throws Exception {
        String tmpSmiles = "CC(=O)N"; // Acetamide CID: 178
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
     * Tests method for descriptor BOND_COUNT_ALL.
     */
    @Test
    public void test_BOND_COUNT_ALL() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Preparation of the first molecule
        String tmpSmiles1 = "CCO"; // Ethanol CID: 702
        IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
        // Preparation of the second molecule
        String tmpSmiles2 = "C=C=C"; // Allene CID: 10037
        IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
        // Preparation of the third molecule
        String tmpSmiles3 = "CC#N"; // Acetonitrile CID: 6342
        IAtomContainer tmpMolecule3 = tmpSmilesParser.parseSmiles(tmpSmiles3);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2, tmpMolecule3};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_ALL};
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}, {0f}, {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[1][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[2][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}, {0f}, {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[1][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[2][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}, {0f}, {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[1][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[2][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}, {0f}, {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[1][0]));
            Assertions.assertEquals("2", tmpFormat.format(tmpMatrix[2][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptors BOND_COUNT_SINGLE, BOND_COUNT_DOUBLE, and BOND_COUNT_TRIPLE.
     */
    @Test
    public void test_BOND_COUNT_SPECIFIED() throws Exception {
        String tmpSmiles = "C=CC#N"; // Acrylonitrile CID: 7855
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{
                Descriptor.BOND_COUNT_SINGLE,
                Descriptor.BOND_COUNT_DOUBLE,
                Descriptor.BOND_COUNT_TRIPLE
        };
        boolean tmpIsParallelCalculation = false;
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        try {
            Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0])); // Single bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][1])); // Double bonds
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2])); // Triple bonds
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor B_POL.
     */
    @Test
    public void test_B_POL() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(7.517242, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor RULE_OF_FIVE.
     */
    @Test
    public void test_RULE_OF_FIVE() throws Exception {
        String tmpSmiles = "CCCC(OCC)OCC(c1cccc2ccccc12)C4CCC(CCCO)C(CC3CNCNC3)C4"; // Not known in the PubChem database, but a valid SMILES
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.RULE_OF_FIVE};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("3", tmpFormat.format(tmpMatrix[0][0]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor AROMATIC_ATOMS_COUNT.
     */
    @Test
    public void test_AROMATIC_ATOMS_COUNT() throws Exception {
        String tmpSmiles = "c1ccccc1"; // Benzene CID: 241
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Array of all available electron donation models
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AROMATIC_ATOMS_COUNT};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));
            } catch (Exception anException) {
                Assertions.fail("Failed with model " + model.getClass().getSimpleName() + ": " + anException.getMessage());
            }
        }
    }

    /**
     * Tests method for descriptor AROMATIC_BONDS_COUNT with benzene.
     */
    @Test
    public void test_AROMATIC_BONDS_COUNT() throws Exception {
        String tmpSmiles = "c1ccccc1"; // Benzene CID: 241
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AROMATIC_BONDS_COUNT};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("6", tmpFormat.format(tmpMatrix[0][0]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor ROTATABLE_BONDS_COUNT.
     */
    @Test
    public void test_ROTATABLE_BONDS_COUNT() throws Exception {
        String tmpSmiles = "CCNC(=O)CC(C)C"; // N-ethyl-3-methylbutanamide CID: 528605
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("4", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor FMF.
     */
    @Test
    public void test_FMF() throws Exception {
        String tmpSmiles = "Clc1cc(cc(Cl)c1N)C(O)CNC(C)(C)C"; // Clenbuterol CID: 2783
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.353, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor FRACTIONAL_CSP3.
     */
    @Test
    public void test_FRACTIONAL_CSP3() throws Exception {
        String tmpSmiles = "CC1=CC=CC(C)=N1"; // 2,6-Dimethylpyridine CID: 7937
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.29, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor HYBRIDIZATION_RATIO.
     */
    @Test
    public void test_HYBRIDIZATION_RATIO() throws Exception {
        String tmpSmiles = "CCC"; // Propane CID: 6334
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor KAPPA_SHAPE_INDICES.
     */
    @Test
    public void test_KAPPA_SHAPE_INDICES() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(5.0, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(2.25, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(4.0, tmpMatrix[0][2], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(5.0, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(2.25, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(4.0, tmpMatrix[0][2], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(5.0, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(2.25, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(4.0, tmpMatrix[0][2], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
     * Tests method for descriptor PETITJEAN_NUMBER.
     */
    @Test
    public void test_PETITJEAN_NUMBER() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.33333334, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor SPIRO_ATOM_COUNT.
     */
    @Test
    public void test_SPIRO_ATOM_COUNT() throws Exception {
        String tmpSmiles = "C1CCC2(CC1)CC=C1C=CC=CC1=C2"; // Not known in the PubChem database, but a valid SMILES
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor V_ADJ_MAT.
     */
    @Test
    public void test_V_ADJ_MAT() throws Exception {
        String tmpSmiles = "C1CCC2CCCCC2C1"; // Decalin CID: 7044
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(4.459, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor WEIGHTED_PATH.
     */
    @Test
    public void test_WEIGHTED_PATH() throws Exception {
        String tmpSmiles = "CCCC"; // Butane CID: 7843
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
     * Tests method for descriptor ZAGREB_INDEX.
     */
    @Test
    public void test_ZAGREB_INDEX() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(16, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CARBON_TYPES.
     */
    @Test
    public void test_CARBON_TYPES() throws Exception {
        String tmpSmiles = "CCCC"; // Butane CID: 7843
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
     * Tests method for descriptor A_LOG_P.
     */
    @Test
    public void test_A_LOG_P() throws Exception {
        String tmpSmiles = "CCCCl"; // 1-Chloropropane CID: 10899
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.719", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2.955", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("20.584", tmpFormat.format(tmpMatrix[0][2]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.719", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2.955", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("20.584", tmpFormat.format(tmpMatrix[0][2]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals("1.719", tmpFormat.format(tmpMatrix[0][0]));
            Assertions.assertEquals("2.955", tmpFormat.format(tmpMatrix[0][1]));
            Assertions.assertEquals("20.584", tmpFormat.format(tmpMatrix[0][2]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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
     * Tests method for descriptor X_LOG_P.
     */
    @Test
    public void test_X_LOG_P() throws Exception {
        String tmpSmiles = "O=C(O)C(N)CCCN"; // 2,5-Diaminopentanoic Acid CID: 389
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.X_LOG_P};
            boolean tmpIsParallelCalculation = false;
            double epsilon = 0.1; // tolerance range

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(-3.30, tmpMatrix[0][0], epsilon);
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor JP_LOG_P.
     */
    @Test
    public void test_JP_LOG_P() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.3, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor A_POL.
     */
    @Test
    public void test_A_POL() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
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
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(10.88, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor AUTOCORRELATION_CHARGE.
     */
    @Test
    public void test_AUTOCORRELATION_CHARGE() throws Exception {
        String tmpSmiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_CHARGE};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor AUTOCORRELATION_MASS.
     */
    @Test
    public void test_AUTOCORRELATION_MASS() throws Exception {
        String tmpSmiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_MASS};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotEquals(0f, tmpMatrix[0][0]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][1]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][2]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][3]);
            Assertions.assertNotEquals(0f, tmpMatrix[0][4]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor AUTOCORRELATION_POLARIZABILITY.
     * No validated result because the descriptor itself is not validated in the CDK.
     * Result can be printed to see if descriptor calculates values.
     */
    @Test
    public void test_AUTOCORRELATION_POLARIZABILITY() throws Exception {
        String tmpSmiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_POLARIZABILITY};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            //System.out.println("Synchronized calculation results:");
            //System.out.println(Arrays.toString(tmpMatrix[0]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            //System.out.println("New calculation results:");
            //System.out.println(Arrays.toString(tmpMatrix[0]));

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            //System.out.println("Descriptor parallelization results:");
            //System.out.println(Arrays.toString(tmpMatrix[0]));
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor FRAGMENT_COMPLEXITY.
     */
    @Test
    public void test_FRAGMENT_COMPLEXITY() throws Exception {
        String tmpSmiles = "c1ccc(CCc2ccccc2)cc1"; // 1,2-Diphenylethane CID: 7647
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FRAGMENT_COMPLEXITY};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(659.00, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(659.00, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(659.00, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(659.00, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CHI_CHAIN.
     */
    @Test
    public void test_CHI_CHAIN() throws Exception {
        String tmpSmiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_CHAIN};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(10, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.2887, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][8], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.2887, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][8], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.2887, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][8], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.2887, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][8], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CHI_CLUSTER.
     */
    @Test
    public void test_CHI_CLUSTER() throws Exception {
        String tmpSmiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_CLUSTER};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(8, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.2887, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.1667, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CHI_PATH_CLUSTER.
     */
    @Test
    public void test_CHI_PATH_CLUSTER() throws Exception {
        String tmpSmiles = "C1=C(Cl)C=CC=C1(Cl)"; // 1,3-Dichlorobenzene CID: 10943
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_PATH_CLUSTER};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(6, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.7416, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.0934, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.0202, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.4072, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.5585, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.4376, tmpMatrix[0][5], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.7416, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.0934, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.0202, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.4072, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.5585, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.4376, tmpMatrix[0][5], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.7416, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.0934, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.0202, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.4072, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.5585, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.4376, tmpMatrix[0][5], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.7416, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.0934, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.0202, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.4072, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.5585, tmpMatrix[0][4], epsilon);
            Assertions.assertEquals(0.4376, tmpMatrix[0][5], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CHI_PATH.
     */
    @Test
    public void test_CHI_PATH() throws Exception {
        String tmpSmiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_PATH};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(16, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2.9916, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.8938, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.6825, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.5773, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(2.6927, tmpMatrix[0][8], epsilon);
            Assertions.assertEquals(1.5099, tmpMatrix[0][9], epsilon);
            Assertions.assertEquals(1.1439, tmpMatrix[0][10], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2.9916, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.8938, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.6825, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.5773, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(2.6927, tmpMatrix[0][8], epsilon);
            Assertions.assertEquals(1.5099, tmpMatrix[0][9], epsilon);
            Assertions.assertEquals(1.1439, tmpMatrix[0][10], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2.9916, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.8938, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.6825, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.5773, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(2.6927, tmpMatrix[0][8], epsilon);
            Assertions.assertEquals(1.5099, tmpMatrix[0][9], epsilon);
            Assertions.assertEquals(1.1439, tmpMatrix[0][10], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2.9916, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(1.8938, tmpMatrix[0][1], epsilon);
            Assertions.assertEquals(1.6825, tmpMatrix[0][2], epsilon);
            Assertions.assertEquals(0.5773, tmpMatrix[0][3], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][5], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][6], epsilon);
            Assertions.assertEquals(0.0000, tmpMatrix[0][7], epsilon);
            Assertions.assertEquals(2.6927, tmpMatrix[0][8], epsilon);
            Assertions.assertEquals(1.5099, tmpMatrix[0][9], epsilon);
            Assertions.assertEquals(1.1439, tmpMatrix[0][10], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor FRACTIONAL_PSA.
     * Expected results were calculated by TPSADescriptor / MolecularWeightDescriptor.
     */
    @Test
    public void test_FRACTIONAL_PSA() throws Exception {
        String tmpSmiles = "O=C(O)c1ccncc1"; // Isonicotinic Acid CID: 5922
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FRACTIONAL_PSA};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.001; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.4077, tmpMatrix[0][0], epsilon);


            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.4077, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.4077, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.4077, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor LARGEST_PI_SYSTEM.
     */
    @Test
    public void test_LARGEST_PI_SYSTEM() throws Exception {
        String tmpSmiles = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2"; // Not known in the PubChem database, but a valid SMILES
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LARGEST_PI_SYSTEM};
        boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(8, tmpMatrix[0][0]);

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(8, tmpMatrix[0][0]);

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(8, tmpMatrix[0][0]);

                tmpMatrix = new float[][]
                        {
                                {0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals(8, tmpMatrix[0][0]);

            } catch (Exception anException) {
                Assertions.fail();
            }
    }

    /**
     * Tests method for descriptor SMALL_RING.
     */
    @Test
    public void test_SMALL_RING() throws Exception {
        String tmpSmiles = "O=C1c2ccccc2C(=O)c2cc3cc4ccccc4cc3cc21"; // 5,14-Pentacenedione CID: 10686237
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);

            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
            int tmpStartIndex = 0;
            Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.SMALL_RING};
            boolean tmpIsParallelCalculation = false;

            try {
                Assertions.assertEquals(11, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[][]
                        {
                                {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                        };
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][1]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][3]));

                tmpMatrix = new float[][]
                        {
                                {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][1]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][3]));

                tmpMatrix = new float[][]
                        {
                                {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][1]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][3]));

                tmpMatrix = new float[][]
                        {
                                {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                        };
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][0]));
                Assertions.assertEquals("5", tmpFormat.format(tmpMatrix[0][1]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][2]));
                Assertions.assertEquals("1", tmpFormat.format(tmpMatrix[0][3]));
            } catch (Exception anException) {
                Assertions.fail();
            }
        }
    }

    /**
     * Tests method for descriptor BASIC_GROUP_COUNT.
     */
    @Test
    public void test_BASIC_GROUP_COUNT() throws Exception {
        String tmpSmiles = "NC"; // Methylamine CID: 6329
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BASIC_GROUP_COUNT};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor ACIDIC_GROUP_COUNT.
     */
    @Test
    public void test_ACIDIC_GROUP_COUNT() throws Exception {
        String tmpSmiles = "CC(=O)O"; // Acetic acid CID: 176
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ACIDIC_GROUP_COUNT};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor AMINO_ACID_COUNT.
     */
    @Test
    public void test_AMINO_ACID_COUNT() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Preparation of the first molecule
        String tmpSmiles1 = "N[C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])([C@]([H])(O)C)C(=O)O"; // L-threonyl-L-threonine CID: 11321969
        IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
        // Preparation of the second molecule
        String tmpSmiles2 = "C(C(=O)NCC(=O)O)N"; // Glycylglycine CID: 11163
        IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AMINO_ACID_COUNT};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(20, Descriptor.getNumberOfComponents(tmpDescriptors));


            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][8]);
            Assertions.assertEquals(2, tmpMatrix[0][16]);
            Assertions.assertEquals(2, tmpMatrix[1][8]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][8]);
            Assertions.assertEquals(2, tmpMatrix[0][16]);
            Assertions.assertEquals(2, tmpMatrix[1][8]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][8]);
            Assertions.assertEquals(2, tmpMatrix[0][16]);
            Assertions.assertEquals(2, tmpMatrix[1][8]);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][8]);
            Assertions.assertEquals(2, tmpMatrix[0][16]);
            Assertions.assertEquals(2, tmpMatrix[1][8]);

        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor KIER_HALL_SMARTS.
     */
    @Test
    public void test_KIER_HALL_SMARTS() throws Exception {
        String tmpSmiles = "c1c(CN)cc(CCNC)cc1C(CO)CC(=O)CCOCCCO"; // Not known in the PubChem database, but a valid SMILES
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.KIER_HALL_SMARTS};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(79, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            new float[79]
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][33]);
            Assertions.assertEquals(1, tmpMatrix[0][34]);
            Assertions.assertEquals(1, tmpMatrix[0][35]);
            Assertions.assertEquals(1, tmpMatrix[0][20]);
            Assertions.assertEquals(1, tmpMatrix[0][23]);

            tmpMatrix = new float[][]
                    {
                            new float[79]
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][33]);
            Assertions.assertEquals(1, tmpMatrix[0][34]);
            Assertions.assertEquals(1, tmpMatrix[0][35]);
            Assertions.assertEquals(1, tmpMatrix[0][20]);
            Assertions.assertEquals(1, tmpMatrix[0][23]);

            tmpMatrix = new float[][]
                    {
                            new float[79]
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][33]);
            Assertions.assertEquals(1, tmpMatrix[0][34]);
            Assertions.assertEquals(1, tmpMatrix[0][35]);
            Assertions.assertEquals(1, tmpMatrix[0][20]);
            Assertions.assertEquals(1, tmpMatrix[0][23]);

            tmpMatrix = new float[][]
                    {
                            new float[79]
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(2, tmpMatrix[0][33]);
            Assertions.assertEquals(1, tmpMatrix[0][34]);
            Assertions.assertEquals(1, tmpMatrix[0][35]);
            Assertions.assertEquals(1, tmpMatrix[0][20]);
            Assertions.assertEquals(1, tmpMatrix[0][23]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor ECCENTRIC_CONNECTIVITY_INDEX.
     */
    @Test
    public void test_ECCENTRIC_CONNECTIVITY_INDEX() throws Exception {
        String tmpSmiles = "C[C@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4"; // (1R,4S,5S,8S,9R,12S,13R)-1,5,9-trimethyl-11,14,15,16-tetraoxatetracyclo[10.3.1.04,13.08,13]hexadecan-10-one CID: 98047509
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ECCENTRIC_CONNECTIVITY_INDEX};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(254, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(254, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(254, tmpMatrix[0][0]);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(254, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor MDE.
     */
    @Test
    public void test_MDE() throws Exception {
        String tmpSmiles = "COOC(C)(CO)OOC"; // Not known in the PubChem database, but a valid SMILES
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MDE};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        try {
            Assertions.assertEquals(19, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.0000, tmpMatrix[0][10], epsilon);
            Assertions.assertEquals(1.1547, tmpMatrix[0][11], epsilon);
            Assertions.assertEquals(2.9416, tmpMatrix[0][12], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.0000, tmpMatrix[0][10], epsilon);
            Assertions.assertEquals(1.1547, tmpMatrix[0][11], epsilon);
            Assertions.assertEquals(2.9416, tmpMatrix[0][12], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.0000, tmpMatrix[0][10], epsilon);
            Assertions.assertEquals(1.1547, tmpMatrix[0][11], epsilon);
            Assertions.assertEquals(2.9416, tmpMatrix[0][12], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(0.0000, tmpMatrix[0][10], epsilon);
            Assertions.assertEquals(1.1547, tmpMatrix[0][11], epsilon);
            Assertions.assertEquals(2.9416, tmpMatrix[0][12], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor VABC.
     */
    @Test
    public void test_VABC() throws Exception {
        String tmpSmiles = "COc2ccc1[nH]c(nc1c2)S(=O)Cc3ncc(C)c(OC)c3C"; // Prilosec CID: 4594
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.VABC};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(292.23, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(292.23, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(292.23, tmpMatrix[0][0], epsilon);

            tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(292.23, tmpMatrix[0][0], epsilon);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests PubChem fingerprint functionality.
     */
    @Test
    public void test_PUBCHEM_FINGERPRINT() throws Exception {
        String tmpSmiles = "C1=CC=C(C=C1)C[N+]2=C(C=C(C=C2C=CC3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5"; // 1-Benzyl-2,4-diphenyl-6-(2-phenylethenyl)pyridin-1-ium CID: 3828524
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.PUBCHEM_FINGERPRINTER};
        boolean tmpIsParallelCalculation = false;
        // Verify that expected bits are set (convert reference fingerprint to expected values)
        IFingerprinter printer = new PubchemFingerprinter(SilentChemObjectBuilder.getInstance());
        BitSet ref = PubchemFingerprinter
                .decode("AAADceB+AAAAAAAAAAAAAAAAAAAAAAAAAAA8YMGCAAAAAAAB1AAAHAAAAAAADAjBHgQwgJMMEACgAyRiRACCgCAhAiAI2CA4ZJgIIOLAkZGEIAhggADIyAcQgMAOgAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==");

        try {
            // Test descriptor component count
            Assertions.assertEquals(881, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[1][881];
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());

            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            for (int i = 0; i < 881; i++) {
                float expected = ref.get(i) ? 1.0f : 0.0f;
                Assertions.assertEquals(expected, tmpMatrix[0][i]);
            }

            tmpMatrix = new float[1][881];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());

            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            for (int i = 0; i < 881; i++) {
                float expected = ref.get(i) ? 1.0f : 0.0f;
                Assertions.assertEquals(expected, tmpMatrix[0][i]);
            }

            tmpMatrix = new float[1][881];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());

            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            for (int i = 0; i < 881; i++) {
                float expected = ref.get(i) ? 1.0f : 0.0f;
                Assertions.assertEquals(expected, tmpMatrix[0][i]);
            }

        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor MACCS_FINGERPRINTER.
     */
    @Test
    public void test_MACCS_FINGERPRINTER() throws Exception {
        String tmpSmiles = "c1ccccc1CCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MACCS_FINGERPRINTER};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(166, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[1][166];
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][124]);
            Assertions.assertEquals(0, tmpMatrix[0][165]);

            tmpMatrix = new float[1][166];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][124]);
            Assertions.assertEquals(0, tmpMatrix[0][165]);

            tmpMatrix = new float[1][166];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertEquals(1, tmpMatrix[0][124]);
            Assertions.assertEquals(0, tmpMatrix[0][165]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CIRCULAR_FINGERPRINTER_ECFP.
     */
    @Test
    public void test_CIRCULAR_FINGERPRINTER_ECFP() throws Exception {
        String tmpSmiles = "c1ccccc1CCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0,Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(1024 * 4, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[1][1024 * 4];
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotNull(tmpMatrix);

            tmpMatrix = new float[1][1024 * 4];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotNull(tmpMatrix);

            tmpMatrix = new float[1][1024 * 4];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotNull(tmpMatrix);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Tests method for descriptor CIRCULAR_FINGERPRINTER_FCFP.
     */
    @Test
    public void test_CIRCULAR_FINGERPRINTER_FCFP() throws Exception {
        String tmpSmiles = "c1ccccc1CCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6};
        boolean tmpIsParallelCalculation = false;

        try {
            Assertions.assertEquals(4 * 1024, Descriptor.getNumberOfComponents(tmpDescriptors));

            float[][] tmpMatrix = new float[1][4 * 1024];
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotNull(tmpMatrix);

            tmpMatrix = new float[1][4 * 1024];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotNull(tmpMatrix);

            tmpMatrix = new float[1][4 * 1024];
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );
            Assertions.assertNotNull(tmpMatrix);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    // Add new descriptor tests here!

    //</editor-fold>

    //<editor-fold desc="Combined descriptor tests">
    /**
     * Tests combined descriptors.
     */
    @Test
    public void test_CombinedDescriptorsSequential() throws Exception {
        // Acetic acid
        String tmpSmiles = "CC(=O)O"; // Acetic acid CID: 176
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };
        DecimalFormatSymbols tmpSymbols = new DecimalFormatSymbols(Locale.US);
        DecimalFormat tmpFormat = new DecimalFormat("0.00", tmpSymbols);

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

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

            try {
                Assertions.assertEquals(6, Descriptor.getNumberOfComponents(tmpDescriptors));

                float[][] tmpMatrix = new float[1][];
                tmpMatrix[0] = new float[Descriptor.getNumberOfComponents(tmpDescriptors)];
                List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
                        )
                );
                Assertions.assertEquals("60.05", tmpFormat.format(tmpMatrix[0][0])); //molecular weight
                Assertions.assertEquals("9.00", tmpFormat.format(tmpMatrix[0][1])); //wiener path number
                Assertions.assertEquals("0.00", tmpFormat.format(tmpMatrix[0][2])); //wiener polarity number
                Assertions.assertEquals("8.00", tmpFormat.format(tmpMatrix[0][3])); //atom count
                Assertions.assertEquals("2.00", tmpFormat.format(tmpMatrix[0][4])); //H-bond acceptor count
                Assertions.assertEquals("1.00", tmpFormat.format(tmpMatrix[0][5])); //H-bond donor count

                tmpMatrix = new float[1][];
                tmpMatrix[0] = new float[Descriptor.getNumberOfComponents(tmpDescriptors)];
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
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
                aNanPositions = Collections.synchronizedList(new LinkedList<>());
                Assertions.assertTrue(
                        Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                                tmpDescriptors,
                                tmpMoleculesArray,
                                tmpMatrix,
                                tmpStartIndex,
                                tmpIsParallelCalculation,
                                aNanPositions
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
    }
    //</editor-fold>

    //<editor-fold desc="Tests with all implemented descriptors">
    /**
     * Tests parallelization.
     */
    @Test
    public void test_Parallelization() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];

        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            tmpMoleculesArray[i] = tmpMolecule;
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getAllDescriptors();
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);

        try {

            float[][] tmpMatrixSequential = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixSequential[i] = new float[tmpNumberOfComponents];
            }
            boolean tmpIsParallelCalculation = false;
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            List<int[]> aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
     * Tests parallelization without fingerprints.
     */
    @Test
    public void test_Parallelization_without_Fingerprints() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        String[] tmpMoleculeStringsArray = new String[tmpNumberOfMolecules];

        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            tmpMoleculesArray[i] = tmpMolecule;
            tmpMoleculeStringsArray[i] = tmpSmiles;
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);

        try {

            float[][] tmpMatrixSequential = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixSequential[i] = new float[tmpNumberOfComponents];
            }
            boolean tmpIsParallelCalculation = false;
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            List<int[]> aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            100,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            100,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculeSmilesStringsByBatchParallelization(
                            tmpDescriptors,
                            tmpMoleculeStringsArray,
                            tmpMatrixSequential,
                            tmpStartIndex,
                            100,
                            null,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            tmpMatrixParallel = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrixParallel[i] = new float[tmpNumberOfComponents];
            }
            tmpIsParallelCalculation = true;
            aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculeSmilesStringsByBatchParallelization(
                            tmpDescriptors,
                            tmpMoleculeStringsArray,
                            tmpMatrixParallel,
                            tmpStartIndex,
                            100,
                            null,
                            tmpIsParallelCalculation,
                            aNanPositionsParallel
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
     * Tests integrity.
     */
    @Test
    public void test_Integrity() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];

        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            tmpMoleculesArray[i] = tmpMolecule;
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getAllDescriptors();
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);
        boolean tmpIsParallelCalculation = false;

        try {

            float[][] tmpMatrix1 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix1[i] = new float[tmpNumberOfComponents];
            }
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix1,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix2 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix2[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix2,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix3 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix3[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix3,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
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

    /**
     * Tests integrity without Fingerprints because Descriptor.setDescriptorsForMoleculesByMoleculeParallelization does not support the calculation of Fingerprints.
     */
    @Test
    public void test_Integrity_without_Fingerprints() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        String[] tmpMoleculeStringsArray = new String[tmpNumberOfMolecules];

        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            tmpMoleculesArray[i] = tmpMolecule;
            tmpMoleculeStringsArray[i] = tmpSmiles;
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);
        boolean tmpIsParallelCalculation = false;

        try {

            float[][] tmpMatrix1 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix1[i] = new float[tmpNumberOfComponents];
            }
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix1,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix2 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix2[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix2,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix3 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix3[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByDescriptorParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix3,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix4 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix4[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix4,
                            tmpStartIndex,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix5 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix5[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                            tmpDescriptors,
                            tmpMoleculesArray,
                            tmpMatrix5,
                            tmpStartIndex,
                            100,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            float[][] tmpMatrix6 = new float[tmpNumberOfMolecules][];
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                tmpMatrix6[i] = new float[tmpNumberOfComponents];
            }
            aNanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculeSmilesStringsByBatchParallelization(
                            tmpDescriptors,
                            tmpMoleculeStringsArray,
                            tmpMatrix6,
                            tmpStartIndex,
                            100,
                            null,
                            tmpIsParallelCalculation,
                            aNanPositions
                    )
            );

            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                for (int j = 0; j < tmpNumberOfComponents; j++) {
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix2[i][j]);
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix3[i][j]);
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix4[i][j]);
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix5[i][j]);
                    Assertions.assertEquals(tmpMatrix1[i][j], tmpMatrix6[i][j]);
                }
            }

        } catch (Exception anException) {
            Assertions.fail();
        }
    }
    //</editor-fold>

    //<editor-fold desc="Tests preparation Methods">
    /**
     * Tests createMoleculeWithExplicitHydrogens.
     */
    @Test
    public void testCreateMoleculeWithExplicitHydrogens() throws Exception {
        // Create a simple molecule (methane) with implicit hydrogen atoms
        String methanSmiles = "C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer methaneImplicit = smilesParser.parseSmiles(methanSmiles);

        // Ensure the input is correct (1 atom, no bonds to H)
        Assertions.assertEquals(1, methaneImplicit.getAtomCount(),
                "Original methane should only have 1 atom");
        Assertions.assertEquals(0, methaneImplicit.getBondCount(),
                "Original methane should have no bonds");

        // Convert to molecule with explicit hydrogen atoms
        IAtomContainer methaneExplicit = Descriptor.createMoleculeWithExplicitHydrogens(methaneImplicit);

        // Check that the new molecule has the expected number of atoms (C + 4H = 5)
        Assertions.assertEquals(5, methaneExplicit.getAtomCount(),
                "Methane with explicit H should have 5 atoms (C + 4H)");

        // Check that the new molecule has the expected number of bonds (4 C-H bonds)
        Assertions.assertEquals(4, methaneExplicit.getBondCount(),
                "Methane with explicit H should have 4 bonds");
        // Check count of hydrogens
        int hydrogenCount = 0;
        for (IAtom atom : methaneExplicit.atoms()) {
            if ("H".equals(atom.getSymbol())) {
                hydrogenCount++;
            }
        }
        Assertions.assertEquals(4, hydrogenCount,
                "There should be 4 explicit hydrogen atoms");

        // Test with a more complex molecule
        String ethanolSmiles = "CCO";
        IAtomContainer ethanolImplicit = smilesParser.parseSmiles(ethanolSmiles);

        // Check if input is correct
        Assertions.assertEquals(3, ethanolImplicit.getAtomCount(),
                "Ethanol should have 3 atoms (implicit H)");

        // Convert to molecule with explicit hydrogens
        IAtomContainer ethanolExplicit = Descriptor.createMoleculeWithExplicitHydrogens(ethanolImplicit);

        // Check atom count (C + C + O + 6H = 9)
        Assertions.assertEquals(9, ethanolExplicit.getAtomCount(),
                "Ethanol with explicit H should have 9 atoms");
    }

    /**
     * Tests the copyMolecule method.
     */
    @Test
    public void testCopyMolecule() throws Exception {
        // Create a test molecule (benzene)
        String benzeneSmiles = "c1ccccc1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer benzene = smilesParser.parseSmiles(benzeneSmiles);

        // Verify initial state
        Assertions.assertEquals(6, benzene.getAtomCount(), "Original benzene should have 6 atoms");
        Assertions.assertEquals(6, benzene.getBondCount(), "Original benzene should have 6 bonds");

        // Create a copy of the molecule
        IAtomContainer benzeneClone = Descriptor.copyMolecule(benzene);

        // Verify the copy has the same structure
        Assertions.assertEquals(benzene.getAtomCount(), benzeneClone.getAtomCount(),
                "Clone should have the same number of atoms as the original");
        Assertions.assertEquals(benzene.getBondCount(), benzeneClone.getBondCount(),
                "Clone should have the same number of bonds as the original");

        // Ensure they are separate objects (deep copy)
        Assertions.assertNotSame(benzene, benzeneClone, "Clone should be a different object instance");

        // Modify the original to confirm the clone is independent
        IAtom newAtom = benzene.getBuilder().newInstance(IAtom.class, "O");
        benzene.addAtom(newAtom);

        // Verify the clone remains unchanged
        Assertions.assertEquals(7, benzene.getAtomCount(), "Modified original should have 7 atoms");
        Assertions.assertEquals(6, benzeneClone.getAtomCount(), "Clone should still have 6 atoms");
    }

    /**
     * Tests setAromaticity.
     */
    @Test
    public void testSetAromaticity() throws Exception {
        // Test molecule: Benzene
        String benzeneSmiles = "c1ccccc1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer benzene = smilesParser.parseSmiles(benzeneSmiles);

        // Aromaticity models
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Fresh molecule for each test
            IAtomContainer testMolecule = benzene.clone();

            // Critical step: Perceive atom types and configure atoms
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(testMolecule);

            // Reset aromaticity flags
            Aromaticity.clear(testMolecule);

            // Set aromaticity using the current model
            // Note: This test still uses the old Descriptor.setAromaticity which takes ElectronDonation
            // It might need updating if Descriptor.setAromaticity now expects Aromaticity.Model
            Descriptor.setAromaticity(testMolecule, model);

            // Check if aromatic elements are present
            boolean hasAromaticAtoms = false;
            boolean hasAromaticBonds = false;

            for (IAtom atom : testMolecule.atoms()) {
                if (atom.isAromatic()) {
                    hasAromaticAtoms = true;
                    break;
                }
            }

            for (IBond bond : testMolecule.bonds()) {
                if (bond.isAromatic()) {
                    hasAromaticBonds = true;
                    break;
                }
            }

            String modelName = model.getClass().getSimpleName();
            Assertions.assertTrue(hasAromaticAtoms,
                    "Model " + modelName + " should identify aromatic atoms");
            Assertions.assertTrue(hasAromaticBonds,
                    "Model " + modelName + " should identify aromatic bonds");
        }
    }

    /**
     * Tests getAllFingerprints.
     */
    @Test
    public void test_getAllFingerprints() throws Exception {
        Descriptor[] allFingerprints = Descriptor.getAllFingerprints();

        // Check if at least one fingerprint is returned
        Assertions.assertTrue(allFingerprints.length > 0);

        // Check if all returned descriptors are indeed fingerprints
        for (Descriptor descriptor : allFingerprints) {
            Assertions.assertTrue(Descriptor.isFingerprint(descriptor));
        }

        // Check for presence of specific known fingerprints
        List<Descriptor> fingerprintList = java.util.Arrays.asList(allFingerprints);
        Assertions.assertTrue(fingerprintList.contains(Descriptor.PUBCHEM_FINGERPRINTER));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.MACCS_FINGERPRINTER));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6));
        Assertions.assertTrue(fingerprintList.contains(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6));
    }

    //</editor-fold>
}

/*
 * CDK-Descriptor-Calculation
 * Copyright (C) 2026 Manuel Schauer, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

//TODO test descriptors for how they handle empty molecules and empty SMILES strings.
/**
 * Test class for Descriptor class. This class first tests all the descriptors included in the {@link Descriptor} class
 * individually, i.e. whether they produce the expected results for some example molecules (in most cases taken from the
 * respective CDK descriptor test class) in different parallelization settings.
 * Note: For adding tests of a new descriptor goto "Add new descriptor tests here!"
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
class DescriptorTest {
    /**
     * Test method for descriptor MOLECULAR_WEIGHT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MOLECULAR_WEIGHT() throws Exception {
        //TODO for CDK integration: remove tmp-prefixes
        String tmpSmiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 60.05f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor WIENER_NUMBER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_WIENER_NUMBER() throws Exception {
        String tmpSmiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.WIENER_NUMBER};
        boolean tmpIsParallelCalculation = false;
        float tmpExpectedWienerPath = 9.0f; // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
        float tmpExpectedWienerPolarity = 0.0f; // there are no atoms that are 3 bonds apart

        Assertions.assertEquals(2, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpectedWienerPath, tmpMatrix[0][0]);
        Assertions.assertEquals(tmpExpectedWienerPolarity, tmpMatrix[0][1]);

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
        Assertions.assertEquals(tmpExpectedWienerPath, tmpMatrix[0][0]);
        Assertions.assertEquals(tmpExpectedWienerPolarity, tmpMatrix[0][1]);
    }

    //TODO: the CDK definitely needs a heavy atom count descriptor or a way to configure this descriptor to only count heavy atoms. Let's discuss how to best realise this.
    /**
     * Test method for descriptor ATOM_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ATOM_COUNT() throws Exception {
        String tmpSmiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ATOM_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 8.0f; // Acetic acid has 8 atoms (implicit Hs included)

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };

        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for the organic subset of individual atom counts (C, H, N, O, S, P, F, Br, Cl, I) in a complex molecule.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ATOM_COUNT_ORGANIC_SUBSET() throws Exception {
        // Cobalamin CID: 74413906
        String tmpSmiles = "CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)OP(=O)([O-])OC(C)CNC(=O)CCC4(C(C5C6(C(C(C(=N6)C(=C7C(C(C(=N7)C=C8C(C(C(=N8)C(=C4[N-]5)C)CCC(=O)N)(C)C)CCC(=O)N)(C)CC(=O)N)C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[Co+3]";
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

        Assertions.assertEquals(10, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[] tmpExpected = new float[]{
                62, // C count
                88, // H count
                13, // N count
                14, // O count
                0, // S count
                1, // P count
                0, // F count
                0, // Br count
                0, // Cl count
                0 // I count
        };

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertArrayEquals(tmpExpected, tmpMatrix[0]);

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
        Assertions.assertArrayEquals(tmpExpected, tmpMatrix[0]);
    }

    /**
     * Test method for descriptor H_BOND_ACCEPTOR_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_H_BOND_ACCEPTOR_COUNT() throws Exception {
        String tmpSmiles1 = "CC(=O)O"; //Acetic Acid CID: 176
        String tmpSmiles2 = "O=N(=O)c1cccc2cn[nH]c12"; // 7-Nitroindole CID: 1893
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_ACCEPTOR_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpectedMol1 = 2f;
        float tmpExpectedMol2 = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model for the first molecule
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            //TODO: setAromaticity clears all existing flags, right? In that case, we do not need a new instance for every iteration.
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);

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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);
        }
    }

    /**
     * Test method for descriptor H_BOND_DONOR_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_H_BOND_DONOR_COUNT() throws Exception {
        String tmpSmiles1 = "CC(=O)O"; //Acetic Acid CID: 176
        String tmpSmiles2 = "Oc1ccccc1"; // Phenol CID: 996
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.H_BOND_DONOR_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpectedMol1 = 1f; // Acetic acid has 1 hydrogen bond donor
        float tmpExpectedMol2 = 1f; // Phenol has 1 hydrogen bond donor

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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
            // Parse fresh molecule for each model for the first molecule TODO: see above
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);

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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);
        }
    }

    /**
     * Test method for descriptor TPSA.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_TPSA() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // 3-(dimethylamino)-3-(methyleneamino)propanenitrile (not in PubChem)
        String tmpSmiles1 = "C=NC(CC#N)N(C)C";
        // 1-Nitropropane CID: 7903
        String tmpSmiles2 = "CCCN(=O)=O";

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.TPSA};
        boolean tmpIsParallelCalculation = false;

        double epsilon = 0.01; //tolerance range
        float tmpExpectedMol1 = 39.39f;
        float tmpExpectedMol2 = 45.82f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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
            // Parse fresh molecule for each model for the first molecule TODO: see above
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0], epsilon);

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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0], epsilon);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0], epsilon);
        }
    }

    /**
     * Test method for descriptor LARGEST_CHAIN.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LARGEST_CHAIN() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // 1-Phenylbutadiene CID: 137048
        String tmpSmiles1 = "C=CC=Cc1ccccc1";
        // 4-(4-(penta-2,4-dien-1-yl)benzyl)-3-vinylpyridine (not in PubChem)
        String tmpSmiles2 = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2";

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LARGEST_CHAIN};
        boolean tmpIsParallelCalculation = false;
        float tmpExpectedMol1 = 4f;
        float tmpExpectedMol2 = 5f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        //TODO: does this descriptor really need aromaticity info?
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
            // Parse fresh molecule for each model for the first molecule TODO see bove
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);

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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);
        }
    }

    /**
     * Test method for descriptor LONGEST_ALIPHATIC_CHAIN.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LONGEST_ALIPHATIC_CHAIN() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Butylbenzene CID: 7705
        String tmpSmiles1 = "CCCCc1ccccc1";
        // 4-(4-tert-butylphenoxy)-N-(1,3-thiazol-2-yl)butanamide CID: 1565007
        String tmpSmiles2 = "CC(C)(C)c2ccc(OCCCC(=O)Nc1nccs1)cc2";

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LONGEST_ALIPHATIC_CHAIN};
        boolean tmpIsParallelCalculation = false;
        float tmpExpectedMol1 = 4f;
        float tmpExpectedMol2 = 4f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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
            // Parse fresh molecule for each model for the first molecule TODO see above
            IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(tmpMolecule1, model);
            // Parse fresh molecule for each model for the second molecule
            IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(tmpMolecule2, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);

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
            Assertions.assertEquals(tmpExpectedMol1, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpectedMol2, tmpMatrix[1][0]);
        }
    }

    /**
     * Test method for descriptor MANNHOLD_LOGP.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MANNHOLD_LOGP() throws Exception {
        String tmpSmiles = "C"; // Methane CID: 297
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MANNHOLD_LOGP};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; //tolerance range
        float tmpExpected = 1.57f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor BCUT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BCUT() throws Exception {
        String tmpSmiles = "CC(=O)N"; // Acetamide CID: 178
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BCUT};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.00001; // tolerance range

        float[] tmpExpected = new float[] {11.881587f, 16.005958f, -0.381844f, 0.325509f, 3.374638f, 5.033583f};

        Assertions.assertEquals(6, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor BOND_COUNT_ALL.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BOND_COUNT_ALL() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Ethanol CID: 702
        String tmpSmiles1 = "CCO";
        IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
        // Allene CID: 10037
        String tmpSmiles2 = "C=C=C";
        IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
        // Acetonitrile CID: 6342
        String tmpSmiles3 = "CC#N";
        IAtomContainer tmpMolecule3 = tmpSmilesParser.parseSmiles(tmpSmiles3);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2, tmpMolecule3};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BOND_COUNT_ALL};
        boolean tmpIsParallelCalculation = false;
        //same total bond count for all three molecules
        float tmpExpected = 2f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}, {0f}, {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
        Assertions.assertEquals(tmpExpected, tmpMatrix[1][0]);
        Assertions.assertEquals(tmpExpected, tmpMatrix[2][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
        Assertions.assertEquals(tmpExpected, tmpMatrix[1][0]);
        Assertions.assertEquals(tmpExpected, tmpMatrix[2][0]);
    }

    /**
     * Test method for descriptors BOND_COUNT_SINGLE, BOND_COUNT_DOUBLE, and BOND_COUNT_TRIPLE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BOND_COUNT_SPECIFIED() throws Exception {
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

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

        float tmpExpectedSingleBonds = 1f;
        float tmpExpectedDoubleBonds = 1f;
        float tmpExpectedTripleBonds = 1f;

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpectedSingleBonds, tmpMatrix[0][0]);
        Assertions.assertEquals(tmpExpectedDoubleBonds, tmpMatrix[0][1]);
        Assertions.assertEquals(tmpExpectedTripleBonds, tmpMatrix[0][2]);

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
        Assertions.assertEquals(tmpExpectedSingleBonds, tmpMatrix[0][0]);
        Assertions.assertEquals(tmpExpectedDoubleBonds, tmpMatrix[0][1]);
        Assertions.assertEquals(tmpExpectedTripleBonds, tmpMatrix[0][2]);
    }

    /**
     * Test method for descriptor B_POL.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_B_POL() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.B_POL};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 7.517242f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor RULE_OF_FIVE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_RULE_OF_FIVE() throws Exception {
        // 3-(4-(2-(1-ethoxybutoxy)-1-(naphthalen-1-yl)ethyl)-2-((hexahydropyrimidin-5-yl)methyl)cyclohexyl)propan-1-ol (not in PubChem)
        String tmpSmiles = "CCCC(OCC)OCC(c1cccc2ccccc12)C4CCC(CCCO)C(CC3CNCNC3)C4";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.RULE_OF_FIVE};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 3f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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
            // Parse fresh molecule for each model TODO: again, necessary?
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
        }
    }

    /**
     * Test method for descriptor AROMATIC_ATOMS_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AROMATIC_ATOMS_COUNT() throws Exception {
        String tmpSmiles = "c1ccccc1"; // Benzene CID: 241
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AROMATIC_ATOMS_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 6f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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

        for (ElectronDonation model : models) {
            // Parse fresh molecule for each model TODO see above
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
        }
    }

    /**
     * Test method for descriptor AROMATIC_BONDS_COUNT with benzene.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AROMATIC_BONDS_COUNT() throws Exception {
        String tmpSmiles = "c1ccccc1"; // Benzene CID: 241
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AROMATIC_BONDS_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 6f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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
            // Parse fresh molecule for each model TODO see above
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
        }
    }

    /**
     * Test method for descriptor ROTATABLE_BONDS_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ROTATABLE_BONDS_COUNT() throws Exception {
        // N-ethyl-3-methylbutanamide CID: 528605
        String tmpSmiles = "CCNC(=O)CC(C)C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ROTATABLE_BONDS_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 4f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for descriptor FMF.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FMF() throws Exception {
        // Clenbuterol CID: 2783
        String tmpSmiles = "Clc1cc(cc(Cl)c1N)C(O)CNC(C)(C)C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FMF};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 0.353f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor FRACTIONAL_CSP3.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FRACTIONAL_CSP3() throws Exception {
        // 2,6-Dimethylpyridine CID: 7937
        String tmpSmiles = "CC1=CC=CC(C)=N1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FRACTIONAL_CSP3};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 0.29f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor HYBRIDIZATION_RATIO.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_HYBRIDIZATION_RATIO() throws Exception {
        String tmpSmiles = "CCC"; // Propane CID: 6334
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.HYBRIDIZATION_RATIO};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 1.00f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for descriptor KAPPA_SHAPE_INDICES.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_KAPPA_SHAPE_INDICES() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.KAPPA_SHAPE_INDICES};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] tmpExpected = new float[] {5.0f, 2.25f, 4.0f};

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor PETITJEAN_NUMBER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_PETITJEAN_NUMBER() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.PETITJEAN_NUMBER};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 0.33333334f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor SPIRO_ATOM_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_SPIRO_ATOM_COUNT() throws Exception {
        // 3'H-spiro[cyclohexane-1,2'-naphthalene] (not in PubChem)
        String tmpSmiles = "C1CCC2(CC1)CC=C1C=CC=CC1=C2";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.SPIRO_ATOM_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for descriptor V_ADJ_MAT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_V_ADJ_MAT() throws Exception {
        String tmpSmiles = "C1CCC2CCCCC2C1"; // Decalin CID: 7044
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.V_ADJ_MAT};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.001; // tolerance range TODO: I'm only realising now that the epsilon is different in most cases; what is the rationale behind choosing the value?
        float tmpExpected = 4.459f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor WEIGHTED_PATH.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_WEIGHTED_PATH() throws Exception {
        String tmpSmiles = "CCCC"; // Butane CID: 7843
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.WEIGHTED_PATH};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.00001; // tolerance range TODO: a very small tolerance range, is there a reason for it? Or are the others just really large?

        float[] tmpExpected = new float[] {6.87132f, 1.71783f, 0.0f, 0.0f, 0.0f};

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor ZAGREB_INDEX.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ZAGREB_INDEX() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ZAGREB_INDEX};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range
        float tmpExpected = 16f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor CARBON_TYPES.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CARBON_TYPES() throws Exception {
        String tmpSmiles = "CCCC"; // Butane CID: 7843
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CARBON_TYPES};
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(9, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[] tmpExpected = new float[] {0f, 0f, 0f, 0f, 0f, 2f, 2f, 0f, 0f};

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertArrayEquals(tmpExpected, tmpMatrix[0]);

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
        Assertions.assertArrayEquals(tmpExpected, tmpMatrix[0]);
    }

    /**
     * Test method for descriptor A_LOG_P.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_A_LOG_P() throws Exception {
        String tmpSmiles = "CCCCl"; // 1-Chloropropane CID: 10899
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[] {tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.A_LOG_P};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.001; //tolerance range

        float[] tmpExpected = new float[] {1.719f, 2.955f, 20.584f};

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor X_LOG_P.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_X_LOG_P() throws Exception {
        String tmpSmiles = "O=C(O)C(N)CCCN"; // 2,5-Diaminopentanoic Acid CID: 389
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.X_LOG_P};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.1; // tolerance range TODO: here, for example, is the tolerance so big because of the different models used?
        float tmpExpected = -3.30f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

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
            // Parse fresh molecule for each model TODO see above
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
            Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
        }
    }

    /**
     * Test method for descriptor JP_LOG_P.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_JP_LOG_P() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.JP_LOG_P};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.1; // tolerance range
        float tmpExpected = 0.3f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor A_POL.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_A_POL() throws Exception {
        String tmpSmiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.A_POL};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 10.88f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

        tmpMatrix = new float[][]
                {
                        {0f} //TODO: BTW, does the matrix necessarily have to be initialised with zeros to use the library?
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor AUTOCORRELATION_CHARGE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_CHARGE() throws Exception {
        String tmpSmiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_CHARGE};
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        //TODO: I know, this is taken like that from the CDK tests; but could we maybe test for the actual values instead of != 0?
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
    }

    /**
     * Test method for descriptor AUTOCORRELATION_MASS.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_MASS() throws Exception {
        String tmpSmiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_MASS};
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        //TODO see above
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
    }

    /**
     * Test method for descriptor AUTOCORRELATION_POLARIZABILITY.
     * No validated result because the descriptor itself is not validated in the CDK.
     * Result can be printed to see if descriptor calculates values.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_POLARIZABILITY() throws Exception {
        String tmpSmiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_POLARIZABILITY};
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        //TODO: we need to find a better solution for this, there should be no print-outs; maybe there are example values in the original publication?
        System.out.println("New calculation results:");
        System.out.println(Arrays.toString(tmpMatrix[0]));

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
        System.out.println("Descriptor parallelization results:");
        System.out.println(Arrays.toString(tmpMatrix[0]));
    }

    /**
     * Test method for descriptor FRAGMENT_COMPLEXITY.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FRAGMENT_COMPLEXITY() throws Exception {
        String tmpSmiles = "c1ccc(CCc2ccccc2)cc1"; // 1,2-Diphenylethane CID: 7647
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FRAGMENT_COMPLEXITY};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 659.00f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor CHI_CHAIN.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_CHAIN() throws Exception {
        String tmpSmiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_CHAIN};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] tmpExpected = new float[] {0.2887f, 0.2887f, 0.0000f, 0.0000f, 0.0000f, 0.1667f, 0.1667f, 0.0000f, 0.0000f};

        Assertions.assertEquals(10, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor CHI_CLUSTER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_CLUSTER() throws Exception {
        String tmpSmiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_CLUSTER};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] tmpExpected = new float[]{0.2887f, 0.0000f, 0.0000f, 0.0000f, 0.1667f, 0.0000f, 0.0000f, 0.0000f};

        Assertions.assertEquals(8, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor CHI_PATH_CLUSTER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_PATH_CLUSTER() throws Exception {
        String tmpSmiles = "C1=C(Cl)C=CC=C1(Cl)"; // 1,3-Dichlorobenzene CID: 10943
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_PATH_CLUSTER};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] tmpExpected = new float[]{0.7416f, 1.0934f, 1.0202f, 0.4072f, 0.5585f, 0.4376f};

        Assertions.assertEquals(6, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor CHI_PATH.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_PATH() throws Exception {
        String tmpSmiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.CHI_PATH};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] tmpExpected = new float[] {2.9916f, 1.8938f, 1.6825f, 0.5773f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 2.6927f, 1.5099f, 1.1439f};

        Assertions.assertEquals(16, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }

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
        for (int i = 0; i < tmpExpected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(tmpExpected[i], tmpMatrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor FRACTIONAL_PSA.
     * Expected results were calculated by TPSADescriptor / MolecularWeightDescriptor. TODO: nice!
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FRACTIONAL_PSA() throws Exception {
        String tmpSmiles = "O=C(O)c1ccncc1"; // Isonicotinic Acid CID: 5922
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.FRACTIONAL_PSA};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.001; // tolerance range
        float tmpExpected = 0.4077f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor LARGEST_PI_SYSTEM.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LARGEST_PI_SYSTEM() throws Exception {
        //4-(4-(penta-2,4-dien-1-yl)benzyl)-3-vinylpyridine (not in PubChem)
        String tmpSmiles = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.LARGEST_PI_SYSTEM};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 8f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

    }

    /**
     * Test method for descriptor SMALL_RING.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_SMALL_RING() throws Exception {
        // 5,14-Pentacenedione CID: 10686237
        String tmpSmiles = "O=C1c2ccccc2C(=O)c2cc3cc4ccccc4cc3cc21";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.SMALL_RING};
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(11, Descriptor.getNumberOfComponents(tmpDescriptors));

        //numbers indicate the tested descriptor value, not all values are tested here, only the first 4 positions
        float tmpExpected0 = 5f;
        float tmpExpected1 = 5f;
        float tmpExpected2 = 1f;
        float tmpExpected3 = 1f;

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
            // Parse fresh molecule for each model TODO see above
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);

            // Apply aromaticity with the current model
            Descriptor.setAromaticity(tmpMolecule, model);

            IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};

            float[][] tmpMatrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
            Assertions.assertEquals(tmpExpected0, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpected1, tmpMatrix[0][1]);
            Assertions.assertEquals(tmpExpected2, tmpMatrix[0][2]);
            Assertions.assertEquals(tmpExpected3, tmpMatrix[0][3]);

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
            Assertions.assertEquals(tmpExpected0, tmpMatrix[0][0]);
            Assertions.assertEquals(tmpExpected1, tmpMatrix[0][1]);
            Assertions.assertEquals(tmpExpected2, tmpMatrix[0][2]);
            Assertions.assertEquals(tmpExpected3, tmpMatrix[0][3]);
        }
    }

    /**
     * Test method for descriptor BASIC_GROUP_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BASIC_GROUP_COUNT() throws Exception {
        String tmpSmiles = "NC"; // Methylamine CID: 6329
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.BASIC_GROUP_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for descriptor ACIDIC_GROUP_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ACIDIC_GROUP_COUNT() throws Exception {
        String tmpSmiles = "CC(=O)O"; // Acetic acid CID: 176
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ACIDIC_GROUP_COUNT};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for descriptor AMINO_ACID_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AMINO_ACID_COUNT() throws Exception {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // L-threonyl-L-threonine CID: 11321969
        String tmpSmiles1 = "N[C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])([C@]([H])(O)C)C(=O)O";
        IAtomContainer tmpMolecule1 = tmpSmilesParser.parseSmiles(tmpSmiles1);
        // Glycylglycine CID: 11163
        String tmpSmiles2 = "C(C(=O)NCC(=O)O)N";
        IAtomContainer tmpMolecule2 = tmpSmilesParser.parseSmiles(tmpSmiles2);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule1, tmpMolecule2};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.AMINO_ACID_COUNT};
        boolean tmpIsParallelCalculation = false;
        //expected threonine count does not apply to the 2nd mol
        float tmpExpectedGlycineAndThreonineCount = 2f;

        Assertions.assertEquals(20, Descriptor.getNumberOfComponents(tmpDescriptors));


        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpectedGlycineAndThreonineCount, tmpMatrix[0][8]); //TODO: is this the inherent behaviour of this descriptor, that all amino acid moieties are also detected as Glycine?
        Assertions.assertEquals(tmpExpectedGlycineAndThreonineCount, tmpMatrix[0][16]);
        Assertions.assertEquals(tmpExpectedGlycineAndThreonineCount, tmpMatrix[1][8]);

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
        Assertions.assertEquals(tmpExpectedGlycineAndThreonineCount, tmpMatrix[0][8]);
        Assertions.assertEquals(tmpExpectedGlycineAndThreonineCount, tmpMatrix[0][16]);
        Assertions.assertEquals(tmpExpectedGlycineAndThreonineCount, tmpMatrix[1][8]);

    }

    /**
     * Test method for descriptor KIER_HALL_SMARTS.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_KIER_HALL_SMARTS() throws Exception {
        //5-(3-(aminomethyl)-5-(2-(methylamino)ethyl)phenyl)-6-hydroxy-1-(3-hydroxypropoxy)hexan-3-one (not in PubChem)
        String tmpSmiles = "c1c(CN)cc(CCNC)cc1C(CO)CC(=O)CCOCCCO";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.KIER_HALL_SMARTS};
        boolean tmpIsParallelCalculation = false;
        // numbers indicate the vector positions of the tested descriptor values; not all return values are tested here
        float tmpExpected33 = 2f;
        float tmpExpected34 = 1f;
        float tmpExpected35 = 1f;
        float tmpExpected20 = 1f;
        float tmpExpected23 = 1f;

        Assertions.assertEquals(79, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        new float[79]
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected33, tmpMatrix[0][33]);
        Assertions.assertEquals(tmpExpected34, tmpMatrix[0][34]);
        Assertions.assertEquals(tmpExpected35, tmpMatrix[0][35]);
        Assertions.assertEquals(tmpExpected20, tmpMatrix[0][20]);
        Assertions.assertEquals(tmpExpected23, tmpMatrix[0][23]);

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
        Assertions.assertEquals(tmpExpected33, tmpMatrix[0][33]);
        Assertions.assertEquals(tmpExpected34, tmpMatrix[0][34]);
        Assertions.assertEquals(tmpExpected35, tmpMatrix[0][35]);
        Assertions.assertEquals(tmpExpected20, tmpMatrix[0][20]);
        Assertions.assertEquals(tmpExpected23, tmpMatrix[0][23]);
    }

    /**
     * Test method for descriptor ECCENTRIC_CONNECTIVITY_INDEX.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ECCENTRIC_CONNECTIVITY_INDEX() throws Exception {
        // (1R,4S,5S,8S,9R,12S,13R)-1,5,9-trimethyl-11,14,15,16-tetraoxatetracyclo[10.3.1.04,13.08,13]hexadecan-10-one CID: 98047509
        String tmpSmiles = "C[C@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.ECCENTRIC_CONNECTIVITY_INDEX};
        boolean tmpIsParallelCalculation = false;
        float tmpExpected = 254f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0]);
    }

    /**
     * Test method for descriptor MDE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MDE() throws Exception {
        //2,2-bis(methylperoxy)propan-1-ol (not in PubChem)
        String tmpSmiles = "COOC(C)(CO)OOC";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MDE};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        //numbers indicate the position of the tested descriptor result; not all result values are tested here
        float tmpExpected10 = 0.0000f;
        float tmpExpected11 = 1.1547f;
        float tmpExpected12 = 2.9416f;

        Assertions.assertEquals(19, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected10, tmpMatrix[0][10], epsilon);
        Assertions.assertEquals(tmpExpected11, tmpMatrix[0][11], epsilon);
        Assertions.assertEquals(tmpExpected12, tmpMatrix[0][12], epsilon);

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
        Assertions.assertEquals(tmpExpected10, tmpMatrix[0][10], epsilon);
        Assertions.assertEquals(tmpExpected11, tmpMatrix[0][11], epsilon);
        Assertions.assertEquals(tmpExpected12, tmpMatrix[0][12], epsilon);
    }

    /**
     * Test method for descriptor VABC.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_VABC() throws Exception {
        // Prilosec CID: 4594
        String tmpSmiles = "COc2ccc1[nH]c(nc1c2)S(=O)Cc3ncc(C)c(OC)c3C";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.VABC};
        boolean tmpIsParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float tmpExpected = 292.23f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);

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
        Assertions.assertEquals(tmpExpected, tmpMatrix[0][0], epsilon);
    }

    /**
     * Tests PubChem fingerprint descriptor.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_PUBCHEM_FINGERPRINT() throws Exception {
        // 1-Benzyl-2,4-diphenyl-6-(2-phenylethenyl)pyridin-1-ium CID: 3828524
        String tmpSmiles = "C1=CC=C(C=C1)C[N+]2=C(C=C(C=C2C=CC3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.PUBCHEM_FINGERPRINTER};
        boolean tmpIsParallelCalculation = false;
        // Verify that expected bits are set (convert reference fingerprint to expected values)
        BitSet ref = PubchemFingerprinter.decode(
                "AAADceB+AAAAAAAAAAAAAAAAAAAAAAAAAAA8YMGCAAAAAAAB1AAAHAAAAAAADAjBHgQwgJMMEACgAyRiRACCgCAhAiAI2CA4ZJgIIOLAkZGEIAhggADIyAcQgMAOgAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=="
        );

        // Test descriptor component count
        Assertions.assertEquals(881, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[1][881];
        Arrays.fill(tmpMatrix[0], 0f);
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());

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

        for (int i = 0; i < 881; i++) {
            float expected = ref.get(i) ? 1.0f : 0.0f;
            Assertions.assertEquals(expected, tmpMatrix[0][i]);
        }

        tmpMatrix = new float[1][881];
        Arrays.fill(tmpMatrix[0], 0f);
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

    }

    /**
     * Test method for descriptor MACCS_FINGERPRINTER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MACCS_FINGERPRINTER() throws Exception {
        //1,2-Diphenylethane CID: 7647
        String tmpSmiles = "c1ccccc1CCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MACCS_FINGERPRINTER};
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(166, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[1][166];
        Arrays.fill(tmpMatrix[0], 0f);
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        //TODO: is there no reference fingerprint?
        Assertions.assertEquals(1, tmpMatrix[0][124]);
        Assertions.assertEquals(0, tmpMatrix[0][165]);

        tmpMatrix = new float[1][166];
        Arrays.fill(tmpMatrix[0], 0f);
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
    }

    /**
     * Test method for all CIRCULAR_FINGERPRINTER_ECFP descriptors of different diameters.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CIRCULAR_FINGERPRINTER_ECFP() throws Exception {
        //1,2-Diphenylethane CID: 7647
        String tmpSmiles = "c1ccccc1CCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6
        };
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(1024 * 4, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[1][1024 * 4];
        Arrays.fill(tmpMatrix[0], 0f);
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        //TODO see above
        Assertions.assertEquals(0, aNanPositions.size());

        tmpMatrix = new float[1][1024 * 4];
        Arrays.fill(tmpMatrix[0], 0f);
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
        Assertions.assertEquals(0, aNanPositions.size());
    }

    /**
     * Test method for all CIRCULAR_FINGERPRINTER_FCFP descriptors of different diameters.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CIRCULAR_FINGERPRINTER_FCFP() throws Exception {
        //1,2-Diphenylethane CID: 7647
        String tmpSmiles = "c1ccccc1CCc1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
        Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[]{tmpMolecule};
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6
        };
        boolean tmpIsParallelCalculation = false;

        Assertions.assertEquals(4 * 1024, Descriptor.getNumberOfComponents(tmpDescriptors));

        float[][] tmpMatrix = new float[1][4 * 1024];
        Arrays.fill(tmpMatrix[0], 0f);
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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
        //TODO see above
        Assertions.assertEquals(0, aNanPositions.size());

        tmpMatrix = new float[1][4 * 1024];
        Arrays.fill(tmpMatrix[0], 0f);
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
        Assertions.assertEquals(0, aNanPositions.size());
    }

    // Add new descriptor tests here!



    /**
     * Tests parallelization. TODO: this test can split up into 5. And please add a bit more doc.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    //TODO for CDK integration: tag @Tag("SlowTest") needs to be added here
    void test_Parallelization() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        String[] tmpMoleculeStringsArray = new String[tmpNumberOfMolecules];

        //fill the arrays with 1000 (independent!) instances of propionic acid
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            tmpMoleculesArray[i] = tmpMolecule;
            tmpMoleculeStringsArray[i] = tmpSmiles;
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getAllDescriptors();
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);

        //first, calculate the results sequentially
        float[][] tmpMatrixSequential = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        //TODO: I again have the question here whether the matrix needs to be initialized with zeros or not; if not, we have to remove all my fill() calls above again and you should not initialize the arrays in the test methods above
        boolean tmpIsParallelCalculation = false;
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
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

        //second, test the parallel calculation that uses a new descriptor instance for every calculation
        float[][] tmpMatrixParallel = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = true;
        List<int[]> aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
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

        //compare the results
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            for (int j = 0; j < tmpNumberOfComponents; j++) {
                Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where single molecules are distributed onto threads
        tmpMatrixSequential = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
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

        tmpMatrixParallel = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
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

        //compare the results
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            for (int j = 0; j < tmpNumberOfComponents; j++) {
                Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where batches of molecules are distributed onto threads
        tmpMatrixSequential = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = false;
        int tmpBatchSize = 100;
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        tmpDescriptors,
                        tmpMoleculesArray,
                        tmpMatrixSequential,
                        tmpStartIndex,
                        tmpBatchSize,
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );

        tmpMatrixParallel = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = true;
        aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        tmpDescriptors,
                        tmpMoleculesArray,
                        tmpMatrixParallel,
                        tmpStartIndex,
                        tmpBatchSize,
                        tmpIsParallelCalculation,
                        aNanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            for (int j = 0; j < tmpNumberOfComponents; j++) {
                Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where batches of SMILES strings are distributed onto threads
        tmpMatrixSequential = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = false;
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        tmpDescriptors,
                        tmpMoleculeStringsArray,
                        tmpMatrixSequential,
                        tmpStartIndex,
                        tmpBatchSize,
                        null, //use default aromaticity model
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );

        tmpMatrixParallel = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = true;
        aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        tmpDescriptors,
                        tmpMoleculeStringsArray,
                        tmpMatrixParallel,
                        tmpStartIndex,
                        tmpBatchSize,
                        null, //use default aromaticity model
                        tmpIsParallelCalculation,
                        aNanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            for (int j = 0; j < tmpNumberOfComponents; j++) {
                Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where single SMILES strings are distributed onto threads
        tmpMatrixSequential = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = false;
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
                        tmpDescriptors,
                        tmpMoleculeStringsArray,
                        tmpMatrixSequential,
                        tmpStartIndex,
                        null, //use default aromaticity model
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );

        tmpMatrixParallel = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        tmpIsParallelCalculation = true;
        aNanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
                        tmpDescriptors,
                        tmpMoleculeStringsArray,
                        tmpMatrixParallel,
                        tmpStartIndex,
                        null, //use default aromaticity model
                        tmpIsParallelCalculation,
                        aNanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            for (int j = 0; j < tmpNumberOfComponents; j++) {
                Assertions.assertEquals(tmpMatrixSequential[i][j], tmpMatrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests integrity. TODO: please add a bit more doc.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    //TODO for CDK integration: tag @Tag("SlowTest") needs to be added here
    void test_Integrity() throws Exception {
        String tmpSmiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int tmpNumberOfMolecules = 1000;
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpNumberOfMolecules];
        String[] tmpMoleculeStringsArray = new String[tmpNumberOfMolecules];

        //fill the arrays with 1000 (independent!) instances of propionic acid
        for (int i = 0; i < tmpNumberOfMolecules; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpSmiles);
            Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            tmpMoleculesArray[i] = tmpMolecule;
            tmpMoleculeStringsArray[i] = tmpSmiles;
        }
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = Descriptor.getAllDescriptors();
        int tmpNumberOfComponents = Descriptor.getNumberOfComponents(tmpDescriptors);
        boolean tmpIsParallelCalculation = false;
        int tmpBatchSize = 100;

        float[][] tmpMatrix1 = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        List<int[]> aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        tmpDescriptors,
                        tmpMoleculesArray,
                        tmpMatrix1,
                        tmpStartIndex,
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );

        float[][] tmpMatrix2 = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        tmpDescriptors,
                        tmpMoleculesArray,
                        tmpMatrix2,
                        tmpStartIndex,
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );

        float[][] tmpMatrix3 = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        tmpDescriptors,
                        tmpMoleculesArray,
                        tmpMatrix3,
                        tmpStartIndex,
                        tmpBatchSize,
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );

        float[][] tmpMatrix4 = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        tmpDescriptors,
                        tmpMoleculeStringsArray,
                        tmpMatrix4,
                        tmpStartIndex,
                        tmpBatchSize,
                        null, //use default aromaticity model
                        tmpIsParallelCalculation,
                        aNanPositions
                )
        );
        float[][] tmpMatrix5 = new float[tmpNumberOfMolecules][tmpNumberOfComponents];
        aNanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
                        tmpDescriptors,
                        tmpMoleculeStringsArray,
                        tmpMatrix5,
                        tmpStartIndex,
                        null, //use default aromaticity model
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
            }
        }
    }

    /**
     * Test method for createMoleculeWithExplicitHydrogens.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testCreateMoleculeWithExplicitHydrogens() throws Exception {
        // Create a simple molecule (methane) with implicit hydrogen atoms
        String methaneSmiles = "C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer methaneImplicit = smilesParser.parseSmiles(methaneSmiles);

        // Ensure the input is correct (1 atom, no bonds to H)
        Assertions.assertEquals(1, methaneImplicit.getAtomCount(),
                "Original methane should only have 1 heavy/explicit atom at this point");
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
                "Ethanol should have 3 heavy/explicit atoms here (all Hs are implicit)");

        // Convert to molecule with explicit hydrogens
        IAtomContainer ethanolExplicit = Descriptor.createMoleculeWithExplicitHydrogens(ethanolImplicit);

        // Check atom count (C + C + O + 6H = 9)
        Assertions.assertEquals(9, ethanolExplicit.getAtomCount(),
                "Ethanol with explicit H should have 9 atoms");
    }

    /**
     * Test method for copyMolecule method. TODO: look at my comment on the copy method and extend the test accordingly.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testCopyMolecule() throws Exception {
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
     * Test method for getDescriptorAndComponentIndex and getDescriptorAndComponentInfo methods.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_getDescriptorAndComponentIndexAndInfo() throws Exception {
        Descriptor[] tmpDescriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.BCUT, Descriptor.WIENER_NUMBER};
        // MOLECULAR_WEIGHT (1 component), BCUT (6 components), WIENER_NUMBER (2 components)

        int[] resultIndex = Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, 0);
        Assertions.assertArrayEquals(new int[]{0, 0}, resultIndex);
        String[] resultInfo = Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, 0);
        Assertions.assertArrayEquals(new String[]{Descriptor.MOLECULAR_WEIGHT.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, 1);
        Assertions.assertArrayEquals(new int[]{1, 0}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, 1);
        Assertions.assertArrayEquals(new String[]{Descriptor.BCUT.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, 6);
        Assertions.assertArrayEquals(new int[]{1, 5}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, 6);
        Assertions.assertArrayEquals(new String[]{Descriptor.BCUT.getName(), "5"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, 7);
        Assertions.assertArrayEquals(new int[]{2, 0}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, 7);
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, 8);
        Assertions.assertArrayEquals(new int[]{2, 1}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, 8);
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "1"}, resultInfo);

        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, -1);
        });
        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, -1);
        });

        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentIndex(tmpDescriptors, 9);
        });
        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentInfo(tmpDescriptors, 9);
        });
    }
}

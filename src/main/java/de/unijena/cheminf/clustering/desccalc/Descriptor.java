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

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SpiroAtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.EnumMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * Descriptor related calculations based on the CDK for the enrichment of data vectors.
 * Note. There are 3 "public static boolean setDescriptorsForMolecules...()" methods with different forms of
 * (parallelized) calculation since individual CDK descriptor calculation classes are unfortunately NOT thread-safe.
 * Note: For adding a new descriptor goto "Add new descriptor information here!"
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 */
public enum Descriptor {

    //<editor-fold desc="Descriptor enumeration and initialization">
    /**
     * Molecular weight, adds up the natural masses (weighted average of all known
     * isotopes of the particular element based on their natural abundances) of every atom
     * in the given molecule, it is NOT the exact mass.
     */
    MOLECULAR_WEIGHT,
    /**
     * Wiener number, returns Wiener path number and Wiener polarity number.
     * Path number: sum of the distances between any two atoms in the molecule.
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     */
    WIENER_NUMBER,
    /**
     * Atom count, counts the number of all atoms in the given molecule.
     */
    ATOM_COUNT,
    /**
     * HBondAcceptorCount, counts hydrogen bond acceptors based on a simplified PHACIR scheme.
     * It includes:
     * Oxygen atoms with formal charge ≤ 0, excluding:
     * 		 Aromatic ether oxygens
     *       Oxygens adjacent to nitrogen
     * Nitrogen atoms with formal charge ≤ 0, excluding:
     *       Nitrogens adjacent to oxygen
     */
    H_BOND_ACCEPTOR_COUNT,
    /**
     * HBondDonorCount, counts hydrogen bond donors based on a simplified PHACIR classification.
     * It includes:
     * 		OH groups where the oxygen has a formal charge ≥ 0
     * 		NH groups where the nitrogen has a formal charge ≥ 0
     */
    H_BOND_DONOR_COUNT,
    /**
     * TPSADescriptor, calculates the topological polar surface area (TPSA) of a molecule.
     * TPSA is the sum of the surface areas of polar atoms (typically oxygen and nitrogen)
     * and their attached hydrogens, based on a topological approximation (2D structure only).
     */
    TPSA,
    /**
     * LargestChain descriptor, calculates the number of atoms in the longest chain in the molecule.
     * This is a simple topological descriptor that provides a measure of molecular linearity.
     */
    LARGEST_CHAIN,
    /**
     * LongestAliphaticChain descriptor, calculates the number of atoms in the longest aliphatic chain.
     * This descriptor provides information about the maximum linear extent of non-aromatic
     * portions of the molecular structure, which relates to molecular shape properties.
     */
    LONGEST_ALIPHATIC_CHAIN,
    /**
     * MannholdLogPDescriptor, calculates the octanol-water partition coefficient (logP) using the Mannhold method.
     * LogP describes the hydrophilicity or lipophilicity of a compound and is crucial for
     * predicting solubility, permeability, and bioavailability.
     */
    MANNHOLD_LOGP,
    /**
     * BCUT descriptor, calculates Burden matrix modified eigenvalues with different weighting schemes.
     * 1. BCUTw-1l, BCUTw-2l ... - nhigh lowest atom weighted BCUTS
     * 2. BCUTw-1h, BCUTw-2h ... - nlow highest atom weighted BCUTS
     * 3. BCUTc-1l, BCUTc-2l ... - nhigh lowest partial charge weighted BCUTS
     * 4. BCUTc-1h, BCUTc-2h ... - nlow highest partial charge weighted BCUTS
     * 5. BCUTp-1l, BCUTp-2l ... - nhigh lowest polarizability weighted BCUTS
     * 6. BCUTp-1h, BCUTp-2h ... - nlow highest polarizability weighted BCUTS
     * Note: No array for one parameter is returned, just the highest and lowest numbers
     * nhigh = 1 as default parameter
     * nlow = 1 as default parameter
     */
    BCUT,
    /**
     * BondCount, counts the number of bonds in a molecule with a specific bond order.
     * Default: counts all bonds (total bond count), no bonds to hydrogen atoms are counted.
     */
    BOND_COUNT_ALL,
    /**
     * BondCount Specified, counts the number of bonds in a molecule with specified bond orders.
     * Returns an array with counts for single, double and triple bonds.
     * For aromatic bonds counts use AROMATIC_BONDS_COUNT
     * No bonds to hydrogen atoms are counted.
     */
    BOND_COUNT_SPECIFIED,
    /**
     * Bond polarizability descriptor.
     * The BPolDescriptor calculates the bond polarizability of a molecule.
     * Bond polarizability is a simple sum of polarizability contributions from all bonds, based on bond types and involved atoms.
     * It provides a rough estimate of how easily the electron cloud in a molecule can be distorted,
     * which relates to intermolecular interactions, polarizability and refractive behavior.
     */
    B_POL,
    /**
     * RuleOfFive descriptor, calculates the number of failures of Lipinski's Rule of Five.
     * The descriptor returns the number of violations (0-4).
     */
    RULE_OF_FIVE,
    /**
     * AromaticAtomsCount, counts the number of aromatic atoms in a molecule.
     * Note: Requires that aromatic atoms in the molecule have already been detected and marked.
     */
    AROMATIC_ATOMS_COUNT,
    /**
     * AromaticBondsCount, counts the number of aromatic bonds in a molecule.
     * Note: Requires that aromatic bonds in the molecule have already been detected and marked.
     */
    AROMATIC_BONDS_COUNT,
    /**
     * RotatableBondsCount, counts the number of rotatable bonds in a molecule.
     * A rotatable bond is defined as any single non-ring bond, where atoms on both sides
     * have at least two heavy-atom neighbors. Amide C-N bonds are not counted as rotatable.
     * Excluding terminal bonds.
     * TODO: Include terminal bonds? Exclude Amide C-N bonds?
     */
    ROTATABLE_BONDS_COUNT,
    /**
     * FMF (Framework Match Fraction) descriptor, calculates the ratio of heavy atoms in
     * the framework to the total number of heavy atoms in the molecule.
     * This provides an indication of the proportion of the molecule that is part of the
     * scaffold or core structure, versus the proportion that is in side chains.
     */
    FMF,
    /**
     * FractionalCSP3 descriptor, characterizes the non-flatness of a molecule by calculating
     * the fraction of sp3 hybridized carbon atoms over the total carbon count.
     * This provides information about the three-dimensionality and complexity
     * of a molecule, which relates to drug-likeness properties.
     */
    FRACTIONAL_CSP3,
    /**
     * HybridizationRatio descriptor, calculates the ratio of sp3 carbons to sp2 carbons.
     * This provides valuable information about the three-dimensionality and flatness
     * of a molecule, which can be useful for predicting drug-like properties and
     * comparing structural characteristics.
     */
    HYBRIDIZATION_RATIO,
    /**
     * KappaShapeIndices descriptor, calculates Kier and Hall kappa molecular shape indices.
     * These indices compare the molecular graph with minimal and maximal molecular graphs:
     * Kier1 - First kappa shape index
     * Kier2 - Second kappa shape index
     * Kier3 - Third kappa shape index
     * Note: Hydrogens are ignored in the calculation.
     */
    KAPPA_SHAPE_INDICES,
    /**
     * PetitjeanNumber descriptor, calculates an index characterizing molecular graph topology.
     * This topological descriptor is based on the calculation of the graph eccentricity
     * and provides information about the molecular shape and branching pattern.
     */
    PETITJEAN_NUMBER,
    /**
     * Spiro atom count descriptor
     * Returns the number of spiro atoms in a molecule.
     */
    SPIRO_ATOM_COUNT,
    /**
     * VAdjMa descriptor, calculates the Vertex adjacency information (magnitude).
     * This is calculated as 1 + log2 m, where m is the number of heavy-heavy bonds.
     * If m is zero, then zero is returned.
     * This descriptor characterizes molecular complexity in terms of edge connectivity.
     */
    V_ADJ_MAT,
    /**
     * WeightedPath descriptor, evaluates the weighted path descriptors for a molecule.
     * Returns five values:
     * WTPT1 - molecular ID
     * WTPT2 - molecular ID / number of atoms
     * WTPT3 - sum of path lengths starting from heteroatoms
     * WTPT4 - sum of path lengths starting from oxygens
     * WTPT5 - sum of path lengths starting from nitrogens
     */
    WEIGHTED_PATH,
    /**
     * ZagrebIndex descriptor, calculates the Zagreb index of a molecule.
     * The Zagreb index is the sum of the squares of atom degrees over all heavy atoms,
     * which provides information about the molecular complexity and topological structure.
     */
    ZAGREB_INDEX,
    /**
     * CarbonTypes descriptor, calculates the frequency of occurrence of 9 different types of carbon atoms:
     * C1SP1 - triply bound carbon bound to one other carbon
     * C2SP1 - triply bound carbon bound to two other carbons
     * C1SP2 - doubly bound carbon bound to one other carbon
     * C2SP2 - doubly bound carbon bound to two other carbons
     * C3SP2 - doubly bound carbon bound to three other carbons
     * C1SP3 - singly bound carbon bound to one other carbon
     * C2SP3 - singly bound carbon bound to two other carbons
     * C3SP3 - singly bound carbon bound to three other carbons
     * C4SP3 - singly bound carbon bound to four other carbons
     */
    CARBON_TYPES,
    /**
     * ALogP descriptor, calculates Ghose-Crippen LogP values, molar refractivity values
     * and ALogP squared values.
     * 1. ALogP (logP value) is the Ghose-Crippen octanol-water partition coefficient.
     * 2. ALogP² is the squared ALogP value.
     * 3. Molar Refractivity (MR) measures the volume occupied by an atom or group of atoms.
     */
    A_LOG_P,
    /**
     * XLogP descriptor
     * Prediction of logP based on the atom-type method called XLogP.
     * Requires all hydrogens to be explicit.
     */
    X_LOG_P,
    /**
     * Calculates the JP_LOG_P descriptor (octanol-water partition coefficient based on JPlogP method)
     * Original publication: Junghwan Lee et al. "Estimation of partition coefficients...".
     */
    JP_LOG_P,
    /**
     * APol descriptor, calculates the sum of the atomic polarizabilities (including implicit hydrogens).
     */
    A_POL;

    // Add new descriptor information here!

    /**
     * EnumMap that maps a descriptor to its number of calculated components
     */
    private static final EnumMap<Descriptor, Integer> descriptorToComponentNumberMap = new EnumMap<>(Descriptor.class);
    //<editor-fold desc="Private static final LOGGER">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(Descriptor.class.getName());
    //</editor-fold>
    /**
     * EnumMap that maps a descriptor to an instance of its CDK descriptor class
     */
    private static final EnumMap<Descriptor, IMolecularDescriptor> descriptorToCdkObjectMap = new EnumMap<>(Descriptor.class);
    static {
        try {
            // MOLECULER_WEIGHT has 1 component
            descriptorToComponentNumberMap.put(MOLECULAR_WEIGHT, 1);
            descriptorToCdkObjectMap.put(MOLECULAR_WEIGHT, new WeightDescriptor());

            // WIENER_NUMBER has 2 components, Wiener path number and Wiener polarity number
            descriptorToComponentNumberMap.put(WIENER_NUMBER, 2);
            descriptorToCdkObjectMap.put(WIENER_NUMBER, new WienerNumbersDescriptor());

            // ATOM_COUNT has 1 component
            descriptorToComponentNumberMap.put(ATOM_COUNT, 1);
            descriptorToCdkObjectMap.put(ATOM_COUNT, new AtomCountDescriptor());

            // H_BOND_ACCEPTOR_COUNT has 1 component
            descriptorToComponentNumberMap.put(H_BOND_ACCEPTOR_COUNT, 1);
            descriptorToCdkObjectMap.put(H_BOND_ACCEPTOR_COUNT, new HBondAcceptorCountDescriptor());

            // H_BOND_DONOR_COUNT has 1 component
            descriptorToComponentNumberMap.put(H_BOND_DONOR_COUNT, 1);
            descriptorToCdkObjectMap.put(H_BOND_DONOR_COUNT, new HBondDonorCountDescriptor());

            // TPSA has 1 component
            descriptorToComponentNumberMap.put(TPSA, 1);
            descriptorToCdkObjectMap.put(TPSA, new TPSADescriptor());

            // LARGEST_CHAIN has 1 component
            descriptorToComponentNumberMap.put(LARGEST_CHAIN, 1);
            descriptorToCdkObjectMap.put(LARGEST_CHAIN, new LargestChainDescriptor());

            // LONGEST_ALIPHATIC_CHAIN has 1 component
            descriptorToComponentNumberMap.put(LONGEST_ALIPHATIC_CHAIN, 1);
            descriptorToCdkObjectMap.put(LONGEST_ALIPHATIC_CHAIN, new LongestAliphaticChainDescriptor());

            // MANNHOLD_LOGP has 1 component
            descriptorToComponentNumberMap.put(MANNHOLD_LOGP, 1);
            descriptorToCdkObjectMap.put(MANNHOLD_LOGP, new MannholdLogPDescriptor());

            // BCUT has 6 components
            descriptorToComponentNumberMap.put(BCUT, 6);
            descriptorToCdkObjectMap.put(BCUT, new BCUTDescriptor());

            // BOND_COUNT has 1 component
            descriptorToComponentNumberMap.put(BOND_COUNT_ALL, 1);
            descriptorToCdkObjectMap.put(BOND_COUNT_ALL, new BondCountDescriptor());

            // BOND_COUNT_SPECIFIED has 3 components (single, double, triple bonds)
            descriptorToComponentNumberMap.put(BOND_COUNT_SPECIFIED, 3);
            descriptorToCdkObjectMap.put(BOND_COUNT_SPECIFIED, new BondCountDescriptor());

            // B_POL has 1 component
            descriptorToComponentNumberMap.put(B_POL, 1);
            descriptorToCdkObjectMap.put(B_POL, new BPolDescriptor());

            // RULE_OF_FIVE has 1 component
            descriptorToComponentNumberMap.put(RULE_OF_FIVE, 1);
            descriptorToCdkObjectMap.put(RULE_OF_FIVE, new RuleOfFiveDescriptor());

            // AROMATIC_ATOMS_COUNT has 1 component
            descriptorToComponentNumberMap.put(AROMATIC_ATOMS_COUNT, 1);
            descriptorToCdkObjectMap.put(AROMATIC_ATOMS_COUNT, new AromaticAtomsCountDescriptor());

            // AROMATIC_BONDS_COUNT has 1 component
            descriptorToComponentNumberMap.put(AROMATIC_BONDS_COUNT, 1);
            descriptorToCdkObjectMap.put(AROMATIC_BONDS_COUNT, new AromaticBondsCountDescriptor());

            // ROTATABLE_BONDS_COUNT has 1 component
            descriptorToComponentNumberMap.put(ROTATABLE_BONDS_COUNT, 1);
            descriptorToCdkObjectMap.put(ROTATABLE_BONDS_COUNT, new RotatableBondsCountDescriptor());

            // FMF has 1 component
            descriptorToComponentNumberMap.put(FMF, 1);
            descriptorToCdkObjectMap.put(FMF, new FMFDescriptor());

            // FRACTIONAL_CSP3 has 1 component
            descriptorToComponentNumberMap.put(FRACTIONAL_CSP3, 1);
            descriptorToCdkObjectMap.put(FRACTIONAL_CSP3, new FractionalCSP3Descriptor());

            // HYBRIDIZATION_RATIO has 1 component
            descriptorToComponentNumberMap.put(HYBRIDIZATION_RATIO, 1);
            descriptorToCdkObjectMap.put(HYBRIDIZATION_RATIO, new HybridizationRatioDescriptor());

            // KAPPA_SHAPE_INDICES has 3 components
            descriptorToComponentNumberMap.put(KAPPA_SHAPE_INDICES, 3);
            descriptorToCdkObjectMap.put(KAPPA_SHAPE_INDICES, new KappaShapeIndicesDescriptor());

            // PETITJEAN_NUMBER has 1 component
            descriptorToComponentNumberMap.put(PETITJEAN_NUMBER, 1);
            descriptorToCdkObjectMap.put(PETITJEAN_NUMBER, new PetitjeanNumberDescriptor());

            // SPIRO_ATOM_COUNT has 1 component
            descriptorToComponentNumberMap.put(SPIRO_ATOM_COUNT, 1);
            descriptorToCdkObjectMap.put(SPIRO_ATOM_COUNT, new SpiroAtomCountDescriptor());

            // V_ADJ_MAT has 1 component
            descriptorToComponentNumberMap.put(V_ADJ_MAT, 1);
            descriptorToCdkObjectMap.put(V_ADJ_MAT, new VAdjMaDescriptor());

            // WEIGHTED_PATH has 5 components
            descriptorToComponentNumberMap.put(WEIGHTED_PATH, 5);
            descriptorToCdkObjectMap.put(WEIGHTED_PATH, new WeightedPathDescriptor());

            // ZAGREB_INDEX has 1 component
            descriptorToComponentNumberMap.put(ZAGREB_INDEX, 1);
            descriptorToCdkObjectMap.put(ZAGREB_INDEX, new ZagrebIndexDescriptor());

            // CARBON_TYPES has 9 components
            descriptorToComponentNumberMap.put(CARBON_TYPES, 9);
            descriptorToCdkObjectMap.put(CARBON_TYPES, new CarbonTypesDescriptor());

            // A_LOG_P has 3 components
            descriptorToComponentNumberMap.put(A_LOG_P, 3);
            descriptorToCdkObjectMap.put(A_LOG_P, new ALOGPDescriptor());

            // X_LOG_P has 1 component
            descriptorToComponentNumberMap.put(X_LOG_P, 1);
            descriptorToCdkObjectMap.put(X_LOG_P, new XLogPDescriptor());

            // JP_LOG_P has 1 component
            descriptorToComponentNumberMap.put(JP_LOG_P, 1);
            descriptorToCdkObjectMap.put(JP_LOG_P, new JPlogPDescriptor());

            // A_POL has 1 component
            descriptorToComponentNumberMap.put(A_POL, 1);
            descriptorToCdkObjectMap.put(A_POL, new APolDescriptor());

            // Add new descriptor information here!

        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Failed to initialize descriptors", anException
            );
        }

    }
    //</editor-fold>

    //<editor-fold desc="Public static methods">
    /**
     * Returns all available descriptors
     *
     * @return All available descriptors
     */
    public static Descriptor[] getAllDescriptors(){
        return values();
    }

    /**
     * Returns specified available descriptors.
     * Note: If both flags are false, then null is returned. If both flags are true, then all descriptors are returned.
     *
     * @param isQuicklyCalculableDescriptorInclusion True: Quickly calculable descriptors are returned, false: Otherwise.
     * @param isSlowlyCalculableDescriptorInclusion True: Slowly calculable descriptors are returned, false: Otherwise.
     * @return Specified descriptors
     */
    public static Descriptor[] getSpecifiedDescriptors(
        boolean isQuicklyCalculableDescriptorInclusion,
        boolean isSlowlyCalculableDescriptorInclusion
    ){
        if (isQuicklyCalculableDescriptorInclusion && isSlowlyCalculableDescriptorInclusion) {
            return values();
        } else if (!isQuicklyCalculableDescriptorInclusion && !isSlowlyCalculableDescriptorInclusion) {
            return null;
        } else if (isQuicklyCalculableDescriptorInclusion) {
            // TODO: Implement after analytical data are evaluated.
            return null;
        } else {
            // MUST be isSlowlyCalculableDescriptorInclusion only
            // TODO: Implement after analytical data are evaluated.
            return null;
        }
    }

    /**
     * Returns sum of number of calculated components of an array of defined descriptors
     *
     * @param aDescriptors Array of descriptors (IS NOT CHANGED)
     * @return Sum of number of calculated components of array of descriptors
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if number of components could not be evaluated: This should never happen.
     */
    public static int getNumberOfComponents(
        Descriptor[] aDescriptors
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        if (aDescriptors == null || aDescriptors.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.getNumberOfComponents: aDescriptors is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.getNumberOfComponents: aDescriptor is null or has length 0.");
        }
        for (Descriptor tmpDescriptor : aDescriptors) {
            if (tmpDescriptor == null) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.getNumberOfComponents: A descriptor in aDescriptors is null."
                );
                throw new IllegalArgumentException("Descriptor.getNumberOfComponents: A single descriptor in aDescriptors is null.");
            }
        }
        //</editor-fold>

        try {
            int tmpTotalNumberOfComponents = 0;
            for (Descriptor tmpDescriptor : aDescriptors) {
                tmpTotalNumberOfComponents += descriptorToComponentNumberMap.get(tmpDescriptor);
            }
            return tmpTotalNumberOfComponents;
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.getNumberOfComponents: Number of components could not be evaluated: This should never happen."
            );
            throw new Exception("Descriptor.getNumberOfComponents: An exception occurred: This should never happen.");
        }
    }

    /**
     * Sets calculated descriptor components in vectors of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecules.
     * Note: Uses synchronized descriptor calculation methods which slows down the calculation.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector. (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *               anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculesByMoleculeParallelizationSynchronized(
        Descriptor[] aDescriptors,
        IAtomContainer[] anAtomContainerArray,
        float[][] aMatrix,
        int aStartIndex,
        boolean anIsParallelCalculation
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        if (aDescriptors == null || aDescriptors.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aDescriptors is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aDescriptor is null or has length 0.");
        }
        for (Descriptor tmpDescriptor : aDescriptors) {
            if (tmpDescriptor == null) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A descriptor in aDescriptors is null."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A descriptor in aDescriptors is null.");
            }
        }
        if (anAtomContainerArray == null || anAtomContainerArray.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: anAtomContainerArray is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: anAtomContainerArray is null or has length 0.");
        }
        for (IAtomContainer tmpMolecule : anAtomContainerArray) {
            if (tmpMolecule == null || tmpMolecule.isEmpty()) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A molecule in anAtomContainerArray is null or empty."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A molecule in anAtomContainerArray is null or empty.");
            }
        }
        if (aMatrix == null || aMatrix.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix is null or has length 0.");
        }
        if (aMatrix.length != anAtomContainerArray.length) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix and anAtomContainerArray must have the same length."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix and anAtomContainerArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A vector is null or has length 0."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A vector is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
            }
            try {
                if (aStartIndex + Descriptor.getNumberOfComponents(aDescriptors) > tmpVector.length) {
                    Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                    );
                    throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
                }
            } catch (Exception anException) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
            }
        }
        //</editor-fold>

        try {
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += descriptorToComponentNumberMap.get(aDescriptors[i]);
            }
            if (anIsParallelCalculation) {
                try {
                    boolean[] tmpIsDescriptorCalculations = new boolean[anAtomContainerArray.length];
                    // Advise by Oracle: Parallel streams should use the common Fork-join pool
                    IntStream.range(0, anAtomContainerArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                tmpIsDescriptorCalculations[i] =
                                    Descriptor.setDescriptorsForSingleMoleculeSynchronized(
                                        aDescriptors,
                                        anAtomContainerArray[i],
                                        aMatrix[i],
                                        tmpStartIndices
                                    );
                            } catch (Exception anException) {
                                tmpIsDescriptorCalculations[i] = false;
                            }
                        }
                    );
                    for (int i = 0; i < anAtomContainerArray.length; i++) {
                        if (!tmpIsDescriptorCalculations[i]) {
                            return false;
                        }
                    }
                    return true;
                } catch (Exception anException) {
                    return false;
                }
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < anAtomContainerArray.length; i++) {
                    if (!Descriptor.setDescriptorsForSingleMoleculeSynchronized(aDescriptors, anAtomContainerArray[i], aMatrix[i], tmpStartIndices)) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationSynchronized: An exception occurred: This should never happen."
            );
            throw anException;
        }
    }

    /**
     * Sets calculated descriptor components in vectors of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecules.
     * Note: Uses a new descriptor instance for EVERY descriptor calculation which slows down the calculation.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector. (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *               anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculesByMoleculeParallelizationNew(
            Descriptor[] aDescriptors,
            IAtomContainer[] anAtomContainerArray,
            float[][] aMatrix,
            int aStartIndex,
            boolean anIsParallelCalculation
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        if (aDescriptors == null || aDescriptors.length == 0) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aDescriptors is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aDescriptor is null or has length 0.");
        }
        for (Descriptor tmpDescriptor : aDescriptors) {
            if (tmpDescriptor == null) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A descriptor in aDescriptors is null."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A descriptor in aDescriptors is null.");
            }
        }
        if (anAtomContainerArray == null || anAtomContainerArray.length == 0) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: anAtomContainerArray is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: anAtomContainerArray is null or has length 0.");
        }
        for (IAtomContainer tmpMolecule : anAtomContainerArray) {
            if (tmpMolecule == null || tmpMolecule.isEmpty()) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A molecule in anAtomContainerArray is null or empty."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A molecule in anAtomContainerArray is null or empty.");
            }
        }
        if (aMatrix == null || aMatrix.length == 0) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix is null or has length 0.");
        }
        if (aMatrix.length != anAtomContainerArray.length) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix and anAtomContainerArray must have the same length."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix and anAtomContainerArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A vector is null or has length 0."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: A vector is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
            }
            try {
                if (aStartIndex + Descriptor.getNumberOfComponents(aDescriptors) > tmpVector.length) {
                    Descriptor.LOGGER.log(
                            Level.SEVERE,
                            "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                    );
                    throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
                }
            } catch (Exception anException) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
            }
        }
        //</editor-fold>

        try {
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += descriptorToComponentNumberMap.get(aDescriptors[i]);
            }
            if (anIsParallelCalculation) {
                try {
                    boolean[] tmpIsDescriptorCalculations = new boolean[anAtomContainerArray.length];
                    // Advise by Oracle: Parallel streams should use the common Fork-join pool
                    IntStream.range(0, anAtomContainerArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                tmpIsDescriptorCalculations[i] =
                                    Descriptor.setDescriptorsForSingleMoleculeNew(
                                        aDescriptors,
                                        anAtomContainerArray[i],
                                        aMatrix[i],
                                        tmpStartIndices
                                    );
                            } catch (Exception anException) {
                                tmpIsDescriptorCalculations[i] = false;
                            }
                        }
                    );
                    for (int i = 0; i < anAtomContainerArray.length; i++) {
                        if (!tmpIsDescriptorCalculations[i]) {
                            return false;
                        }
                    }
                    return true;
                } catch (Exception anException) {
                    return false;
                }
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < anAtomContainerArray.length; i++) {
                    if (!Descriptor.setDescriptorsForSingleMoleculeNew(aDescriptors, anAtomContainerArray[i], aMatrix[i], tmpStartIndices)) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew: An exception occurred: This should never happen."
            );
            throw anException;
        }
    }

    /**
     * Sets calculated descriptor components in vectors of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of descriptor calculations.
     * Note: Parallelization of descriptor calculation is fastest since new descriptor instances for every calculation
     * or synchronized descriptor calculation are avoided.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector. (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *               anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculesByDescriptorParallelization(
            Descriptor[] aDescriptors,
            IAtomContainer[] anAtomContainerArray,
            float[][] aMatrix,
            int aStartIndex,
            boolean anIsParallelCalculation
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        if (aDescriptors == null || aDescriptors.length == 0) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aDescriptors is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aDescriptor is null or has length 0.");
        }
        for (Descriptor tmpDescriptor : aDescriptors) {
            if (tmpDescriptor == null) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: A descriptor in aDescriptors is null."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: A descriptor in aDescriptors is null.");
            }
        }
        if (anAtomContainerArray == null || anAtomContainerArray.length == 0) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: anAtomContainerArray is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: anAtomContainerArray is null or has length 0.");
        }
        for (IAtomContainer tmpMolecule : anAtomContainerArray) {
            if (tmpMolecule == null || tmpMolecule.isEmpty()) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: A molecule in anAtomContainerArray is null or empty."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: A molecule in anAtomContainerArray is null or empty.");
            }
        }
        if (aMatrix == null || aMatrix.length == 0) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aMatrix is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aMatrix is null or has length 0.");
        }
        if (aMatrix.length != anAtomContainerArray.length) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aMatrix and anAtomContainerArray must have the same length."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aMatrix and anAtomContainerArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: A vector is null or has length 0."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: A vector is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aStartIndex is illegal.");
            }
            try {
                if (aStartIndex + Descriptor.getNumberOfComponents(aDescriptors) > tmpVector.length) {
                    Descriptor.LOGGER.log(
                            Level.SEVERE,
                            "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aStartIndex is illegal."
                    );
                    throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aStartIndex is illegal.");
                }
            } catch (Exception anException) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponentsByDescriptorParallelization: aStartIndex is illegal.");
            }
        }
        //</editor-fold>

        try {
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += descriptorToComponentNumberMap.get(aDescriptors[i]);
            }
            if (anIsParallelCalculation) {
                try {
                    boolean[] tmpIsMoleculeCalculations = new boolean[aDescriptors.length];
                    // Advise by Oracle: Parallel streams should use the common Fork-join pool
                    IntStream.range(0, aDescriptors.length).parallel().forEach(
                        i ->
                        {
                            try {
                                tmpIsMoleculeCalculations[i] =
                                    Descriptor.setSingleDescriptorForMolecules(
                                        aDescriptors[i],
                                        anAtomContainerArray,
                                        aMatrix,
                                        tmpStartIndices[i]
                                    );
                            } catch (Exception anException) {
                                tmpIsMoleculeCalculations[i] = false;
                            }
                        }
                    );
                    for (int i = 0; i < aDescriptors.length; i++) {
                        if (!tmpIsMoleculeCalculations[i]) {
                            return false;
                        }
                    }
                    return true;
                } catch (Exception anException) {
                    return false;
                }
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < aDescriptors.length; i++) {
                    if (!Descriptor.setSingleDescriptorForMolecules(aDescriptors[i], anAtomContainerArray, aMatrix, tmpStartIndices[i])) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setDescriptorsForMoleculesByDescriptorParallelization: An exception occurred: This should never happen."
            );
            throw anException;
        }
    }
    //</editor-fold>

    //<editor-fold desc="Private static methods">
    /**
     * Sets calculated descriptor components in aVector (that corresponds to anAtomContainer) at aStartIndices
     * Note: Uses synchronized descriptor calculation methods.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setDescriptorsForSingleMoleculeSynchronized (
        Descriptor[] aDescriptors,
        IAtomContainer anAtomContainer,
        float[] aVector,
        int[] aStartIndices
    ) throws Exception {
        try {
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < aDescriptors.length; i++) {
                if (!Descriptor.setDescriptorSynchronized(aDescriptors[i], anAtomContainer, aVector, aStartIndices[i])) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: An exception occurred: This should never happen."
            );
            throw new Exception("Descriptor.setCalculatedDescriptorComponents: An exception occurred: This should never happen.");
        }
    }

    /**
     * Sets calculated descriptor components in aVector (that corresponds to anAtomContainer) at aStartIndices
     * Note: Uses a new descriptor instance for EVERY descriptor calculation.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setDescriptorsForSingleMoleculeNew(
        Descriptor[] aDescriptors,
        IAtomContainer anAtomContainer,
        float[] aVector,
        int[] aStartIndices
    ) throws Exception {
        try {
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < aDescriptors.length; i++) {
                if (!Descriptor.setDescriptorNew(aDescriptors[i], anAtomContainer, aVector, aStartIndices[i])) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setDescriptorsForSingleMoleculeNew: An exception occurred: This should never happen."
            );
            throw new Exception("Descriptor.setDescriptorsForSingleMoleculeNew: An exception occurred: This should never happen.");
        }
    }

    /**
     * Sets calculated descriptor components in aMatrix (that corresponds to anAtomContainerArray) beginning with
     * aStartIndex.
     * Note: Fast implementation without any new descriptor instances or locks/synchronization.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector. (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *               anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setSingleDescriptorForMolecules(
        Descriptor aDescriptor,
        IAtomContainer[] anAtomContainerArray,
        float[][] aMatrix,
        int aStartIndex
    ) throws Exception {
        try {
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < anAtomContainerArray.length; i++) {
                if (!Descriptor.setDescriptor(aDescriptor, anAtomContainerArray[i], aMatrix[i], aStartIndex)) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setSingleDescriptorForMolecules: An exception occurred: This should never happen."
            );
            throw new Exception("Descriptor.setSingleDescriptorForMolecules: An exception occurred: This should never happen.");
        }
    }

    /**
     * Sets component values of aDescriptor for anAtomContainer in aVector beginning with aStartIndex.
     * Note: This method uses synchronized descriptor calculation methods and is thread-safe in concurrent computing.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setDescriptorSynchronized(
        Descriptor aDescriptor,
        IAtomContainer anAtomContainer,
        float[] aVector,
        int aStartIndex
    ) {
        try {
            switch (aDescriptor) {
                case MOLECULAR_WEIGHT:
                    setMolecularWeight(anAtomContainer, aVector, aStartIndex);
                    break;
                case WIENER_NUMBER:
                    setWienerNumber(anAtomContainer, aVector, aStartIndex);
                    break;
                case ATOM_COUNT:
                    setAtomCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case H_BOND_ACCEPTOR_COUNT:
                    setHBondAcceptorCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case H_BOND_DONOR_COUNT:
                    setHBondDonorCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case TPSA:
                    setTPSA(anAtomContainer, aVector, aStartIndex);
                    break;
                case LARGEST_CHAIN:
                    setLargestChain(anAtomContainer, aVector, aStartIndex);
                    break;
                case LONGEST_ALIPHATIC_CHAIN:
                    setLongestAliphaticChain(anAtomContainer, aVector, aStartIndex);
                    break;
                case MANNHOLD_LOGP:
                    setMannholdLogP(anAtomContainer, aVector, aStartIndex);
                    break;
                case BCUT:
                    setBCUT(anAtomContainer, aVector, aStartIndex);
                    break;
                case BOND_COUNT_ALL:
                    setBondCountAll(anAtomContainer, aVector, aStartIndex);
                    break;
                case BOND_COUNT_SPECIFIED:
                    setBondCountSpecified(anAtomContainer, aVector, aStartIndex);
                    break;
                case B_POL:
                    setBPol(anAtomContainer, aVector, aStartIndex);
                    break;
                case RULE_OF_FIVE:
                    setRuleOfFive(anAtomContainer, aVector, aStartIndex);
                    break;
                case AROMATIC_ATOMS_COUNT:
                    setAromaticAtomsCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case AROMATIC_BONDS_COUNT:
                    setAromaticBondsCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case ROTATABLE_BONDS_COUNT:
                    setRotatableBondsCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case FMF:
                    setFMF(anAtomContainer, aVector, aStartIndex);
                    break;
                case FRACTIONAL_CSP3:
                    setFractionalCSP3(anAtomContainer, aVector, aStartIndex);
                    break;
                case HYBRIDIZATION_RATIO:
                    setHybridizationRatio(anAtomContainer, aVector, aStartIndex);
                    break;
                case KAPPA_SHAPE_INDICES:
                    setKappaShapeIndices(anAtomContainer, aVector, aStartIndex);
                    break;
                case PETITJEAN_NUMBER:
                    setPetitjeanNumber(anAtomContainer, aVector, aStartIndex);
                    break;
                case SPIRO_ATOM_COUNT:
                    setSpiroAtomCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case V_ADJ_MAT:
                    setVAdjMat(anAtomContainer, aVector, aStartIndex);
                    break;
                case WEIGHTED_PATH:
                    setWeightedPath(anAtomContainer, aVector, aStartIndex);
                    break;
                case ZAGREB_INDEX:
                    setZagrebIndex(anAtomContainer, aVector, aStartIndex);
                    break;
                case CARBON_TYPES:
                    setCarbonTypesDescriptor(anAtomContainer, aVector, aStartIndex);
                    break;
                case A_LOG_P:
                    setALogP(anAtomContainer, aVector, aStartIndex);
                    break;
                case X_LOG_P:
                    setXLogP(anAtomContainer, aVector, aStartIndex);
                    break;
                case JP_LOG_P:
                    setJPLogP(anAtomContainer, aVector, aStartIndex);
                    break;
                case A_POL:
                    setAPol(anAtomContainer, aVector, aStartIndex);
                    break;
                // Add new descriptor information here!
                default:
                    throw new UnsupportedOperationException("This descriptor does not have a routine yet!");
            }
            return true;
        } catch (Exception anException) {
            for (int i = aStartIndex; i < descriptorToComponentNumberMap.get(aDescriptor); i++) {
                aVector[i] = Float.NaN;
            }
            return false;
        }
    }

    /**
     * Sets component values of aDescriptor for anAtomContainer in aVector beginning with aStartIndex.
     * Note: This method instantiates a CDK descriptor calculator for EVERY calculation and is thread-safe in
     * concurrent computing.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setDescriptorNew(
            Descriptor aDescriptor,
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            switch (aDescriptor) {
                case MOLECULAR_WEIGHT:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new WeightDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case WIENER_NUMBER:
                    DoubleArrayResult tmpResult = (DoubleArrayResult) (new WienerNumbersDescriptor()).calculate(anAtomContainer).getValue();
                    aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
                    aVector[aStartIndex + 1] = (float) tmpResult.get(1); //Wiener polarity number
                    break;
                case ATOM_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new AtomCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case H_BOND_ACCEPTOR_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new HBondAcceptorCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case H_BOND_DONOR_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new HBondDonorCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case TPSA:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new TPSADescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case LARGEST_CHAIN:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new LargestChainDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case LONGEST_ALIPHATIC_CHAIN:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new LongestAliphaticChainDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case MANNHOLD_LOGP:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new MannholdLogPDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case BCUT:
                    DoubleArrayResult bcutResult = (DoubleArrayResult) (new BCUTDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        aVector[aStartIndex + i] = (float) bcutResult.get(i);
                    }
                    break;
                case BOND_COUNT_ALL:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new BondCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_SPECIFIED:
                    BondCountDescriptor singleBondDesc = new BondCountDescriptor();
                    singleBondDesc.setParameters(new Object[]{"s"});
                    aVector[aStartIndex] = (float) ((IntegerResult) singleBondDesc.calculate(anAtomContainer).getValue()).intValue();

                    BondCountDescriptor doubleBondDesc = new BondCountDescriptor();
                    doubleBondDesc.setParameters(new Object[]{"d"});
                    aVector[aStartIndex + 1] = (float) ((IntegerResult) doubleBondDesc.calculate(anAtomContainer).getValue()).intValue();

                    BondCountDescriptor tripleBondDesc = new BondCountDescriptor();
                    tripleBondDesc.setParameters(new Object[]{"t"});
                    aVector[aStartIndex + 2] = (float) ((IntegerResult) tripleBondDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case B_POL:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new BPolDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case RULE_OF_FIVE:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new RuleOfFiveDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case AROMATIC_ATOMS_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new AromaticAtomsCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case AROMATIC_BONDS_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new AromaticBondsCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ROTATABLE_BONDS_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new RotatableBondsCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case FMF:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new FMFDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case FRACTIONAL_CSP3:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new FractionalCSP3Descriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case HYBRIDIZATION_RATIO:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new HybridizationRatioDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case KAPPA_SHAPE_INDICES:
                    DoubleArrayResult kappaResult = (DoubleArrayResult) (new KappaShapeIndicesDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 3; i++) {
                        aVector[aStartIndex + i] = (float) kappaResult.get(i);
                    }
                    break;
                case PETITJEAN_NUMBER:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new PetitjeanNumberDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case SPIRO_ATOM_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new SpiroAtomCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case V_ADJ_MAT:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new VAdjMaDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case WEIGHTED_PATH:
                    DoubleArrayResult tmpWeightedPathResultNew = (DoubleArrayResult) (new WeightedPathDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpWeightedPathResultNew.get(i);
                    }
                    break;
                case ZAGREB_INDEX:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new ZagrebIndexDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case CARBON_TYPES:
                    IntegerArrayResult carbonTypesResultNew = (IntegerArrayResult) (new CarbonTypesDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 9; i++) {
                        aVector[aStartIndex + i] = carbonTypesResultNew.get(i);
                    }
                    break;
                case A_LOG_P:
                    IAtomContainer ALogPMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    DoubleArrayResult aLogPResult = (DoubleArrayResult) (new ALOGPDescriptor()).calculate(ALogPMoleculeWithExplicitH).getValue();
                    aVector[aStartIndex] = (float) aLogPResult.get(0);// ALogP
                    aVector[aStartIndex + 1] = (float) aLogPResult.get(1);  // ALogP squared
                    aVector[aStartIndex + 2] = (float) aLogPResult.get(2);  // Molar Refractivity
                    break;
                case X_LOG_P:
                    IAtomContainer xLogPMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    aVector[aStartIndex] = (float) ((DoubleResult) (new XLogPDescriptor()).calculate(xLogPMoleculeWithExplicitH).getValue()).doubleValue();
                    break;
                case JP_LOG_P:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new JPlogPDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case A_POL:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new APolDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;

                // Add new descriptor information here!
                default:
                    throw new UnsupportedOperationException("This descriptor does not have a routine yet!");
            }
            return true;
        } catch (Exception anException) {
            for (int i = aStartIndex; i < descriptorToComponentNumberMap.get(aDescriptor); i++) {
                aVector[i] = Float.NaN;
            }
            return false;
        }
    }

    /**
     * Sets component values of aDescriptor for anAtomContainer in aVector beginning with aStartIndex.
     * Note: This method is NOT thread-safe in concurrent computing and is NOT made for parallelized access of the
     * SAME non-thread-safe descriptor, only parallelized access of DIFFERENT non-thread-safe descriptors is advised.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setDescriptor(
            Descriptor aDescriptor,
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            switch (aDescriptor) {
                case MOLECULAR_WEIGHT:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MOLECULAR_WEIGHT).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case WIENER_NUMBER:
                    DoubleArrayResult tmpResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(WIENER_NUMBER).calculate(anAtomContainer).getValue();
                    aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
                    aVector[aStartIndex + 1] = (float) tmpResult.get(1); //Wiener polarity number
                    break;
                case ATOM_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(ATOM_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case H_BOND_ACCEPTOR_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(H_BOND_ACCEPTOR_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case H_BOND_DONOR_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(H_BOND_DONOR_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case TPSA:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(TPSA).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case LARGEST_CHAIN:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(LARGEST_CHAIN).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case LONGEST_ALIPHATIC_CHAIN:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(LONGEST_ALIPHATIC_CHAIN).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case MANNHOLD_LOGP:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MANNHOLD_LOGP).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case BCUT:
                    DoubleArrayResult bcutResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(BCUT).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        aVector[aStartIndex + i] = (float) bcutResult.get(i);
                    }
                    break;
                case BOND_COUNT_ALL:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(BOND_COUNT_ALL).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_SPECIFIED:
                    // Single bonds (s)
                    BondCountDescriptor singleBondDesc = (BondCountDescriptor)descriptorToCdkObjectMap.get(BOND_COUNT_SPECIFIED).getClass().getDeclaredConstructor().newInstance();
                    singleBondDesc.setParameters(new Object[]{"s"});
                    aVector[aStartIndex] = (float) ((IntegerResult) singleBondDesc.calculate(anAtomContainer).getValue()).intValue();

                    // Double bonds (d)
                    BondCountDescriptor doubleBondDesc = (BondCountDescriptor)descriptorToCdkObjectMap.get(BOND_COUNT_SPECIFIED).getClass().getDeclaredConstructor().newInstance();
                    doubleBondDesc.setParameters(new Object[]{"d"});
                    aVector[aStartIndex + 1] = (float) ((IntegerResult) doubleBondDesc.calculate(anAtomContainer).getValue()).intValue();

                    // Triple bonds (t)
                    BondCountDescriptor tripleBondDesc = (BondCountDescriptor)descriptorToCdkObjectMap.get(BOND_COUNT_SPECIFIED).getClass().getDeclaredConstructor().newInstance();
                    tripleBondDesc.setParameters(new Object[]{"t"});
                    aVector[aStartIndex + 2] = (float) ((IntegerResult) tripleBondDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case B_POL:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(B_POL).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case RULE_OF_FIVE:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(RULE_OF_FIVE).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case AROMATIC_ATOMS_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(AROMATIC_ATOMS_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case AROMATIC_BONDS_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(AROMATIC_BONDS_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ROTATABLE_BONDS_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(ROTATABLE_BONDS_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case FMF:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FMF).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case FRACTIONAL_CSP3:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FRACTIONAL_CSP3).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case HYBRIDIZATION_RATIO:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(HYBRIDIZATION_RATIO).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case KAPPA_SHAPE_INDICES:
                    DoubleArrayResult kappaResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(KAPPA_SHAPE_INDICES).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 3; i++) {
                        aVector[aStartIndex + i] = (float) kappaResult.get(i);
                    }
                    break;
                case PETITJEAN_NUMBER:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(PETITJEAN_NUMBER).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case SPIRO_ATOM_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(SPIRO_ATOM_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case V_ADJ_MAT:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(V_ADJ_MAT).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case WEIGHTED_PATH:
                    DoubleArrayResult tmpWeightedPathResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(WEIGHTED_PATH).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpWeightedPathResult.get(i);
                    }
                    break;
                case ZAGREB_INDEX:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(ZAGREB_INDEX).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case CARBON_TYPES:
                    IntegerArrayResult carbonTypesResult = (IntegerArrayResult) descriptorToCdkObjectMap.get(CARBON_TYPES).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 9; i++) {
                        aVector[aStartIndex + i] = carbonTypesResult.get(i);
                    }
                    break;
                case A_LOG_P:
                    IAtomContainer aLogPMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    DoubleArrayResult alogpResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(A_LOG_P).calculate(aLogPMoleculeWithExplicitH).getValue();
                    aVector[aStartIndex] = (float) alogpResult.get(0);      // ALogP
                    aVector[aStartIndex + 1] = (float) alogpResult.get(1);  // ALogP squared
                    aVector[aStartIndex + 2] = (float) alogpResult.get(2);  // Molar Refractivity
                    break;
                case X_LOG_P:
                    IAtomContainer xLogPMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(X_LOG_P).calculate(xLogPMoleculeWithExplicitH).getValue()).doubleValue();
                    break;
                case JP_LOG_P:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(JP_LOG_P).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case A_POL:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(A_POL).calculate(anAtomContainer).getValue()).doubleValue();
                    break;

                // Add new descriptor information here!
                default:
                    throw new UnsupportedOperationException("This descriptor does not have a routine yet!");
            }
            return true;
        } catch (Exception anException) {
            for (int i = aStartIndex; i < descriptorToComponentNumberMap.get(aDescriptor); i++) {
                aVector[i] = Float.NaN;
            }
            return false;
        }
    }
    //</editor-fold>
    //<editor-fold desc="Private static synchronized descriptor calculation methods">
    /**
     * Sets molecular weight
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setMolecularWeight(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MOLECULAR_WEIGHT).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets Wiener number(s)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setWienerNumber(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult tmpResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(WIENER_NUMBER).calculate(anAtomContainer).getValue();
        aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
        aVector[aStartIndex + 1] = (float) tmpResult.get(1); //Wiener polarity number
    }

    /**
     * Sets atom count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAtomCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(ATOM_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets hydrogen bond acceptor count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setHBondAcceptorCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(H_BOND_ACCEPTOR_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets hydrogen bond donor count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setHBondDonorCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(H_BOND_DONOR_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the topological polar surface area (TPSA)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setTPSA(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(TPSA).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets largest chain size
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setLargestChain(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(LARGEST_CHAIN).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets longest aliphatic chain size
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setLongestAliphaticChain(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(LONGEST_ALIPHATIC_CHAIN).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the Mannhold LogP value (octanol-water partition coefficient)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods.
     * Note: Method must be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of the molecule to be filled with calculated descriptor components (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector where the calculated descriptor components will be stored
     */
    private static synchronized void setMannholdLogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MANNHOLD_LOGP).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets BCUT descriptor values
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setBCUT(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(BCUT).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 6; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets bond count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setBondCountAll(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(BOND_COUNT_ALL).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets specific bond counts (single, double, triple)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setBondCountSpecified(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            BondCountDescriptor singleBondDesc = BondCountDescriptor.class.getDeclaredConstructor().newInstance();
            singleBondDesc.setParameters(new Object[]{"s"});
            aVector[aStartIndex] = (float) ((IntegerResult) singleBondDesc.calculate(anAtomContainer).getValue()).intValue();

            // Double bonds (d)
            BondCountDescriptor doubleBondDesc = BondCountDescriptor.class.getDeclaredConstructor().newInstance();
            doubleBondDesc.setParameters(new Object[]{"d"});
            aVector[aStartIndex + 1] = (float) ((IntegerResult) doubleBondDesc.calculate(anAtomContainer).getValue()).intValue();

            // Triple bonds (t)
            BondCountDescriptor tripleBondDesc = BondCountDescriptor.class.getDeclaredConstructor().newInstance();
            tripleBondDesc.setParameters(new Object[]{"t"});
            aVector[aStartIndex + 2] = (float) ((IntegerResult) tripleBondDesc.calculate(anAtomContainer).getValue()).intValue();
        } catch (Exception e) {
            for (int i = 0; i < 3; i++) {
                aVector[aStartIndex + i] = Float.NaN;
            }
        }
    }

    /**
     * Sets bond polarizability value
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setBPol(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(B_POL).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets Lipinski's Rule of Five violations count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     * TODO: Test works check if we need to addExplicitHydrogens to the molecule
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setRuleOfFive(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(RULE_OF_FIVE).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the number of aromatic atoms
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAromaticAtomsCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(AROMATIC_ATOMS_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the number of aromatic bonds
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAromaticBondsCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(AROMATIC_BONDS_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets rotatable bonds count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setRotatableBondsCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(ROTATABLE_BONDS_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the FMF (Framework Match Fraction) value
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setFMF(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FMF).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets the FractionalCSP3 value (fraction of sp3 hybridized carbon atoms)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setFractionalCSP3(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FRACTIONAL_CSP3).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets hybridization ratio (sp3 carbons to sp2 carbons)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setHybridizationRatio(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(HYBRIDIZATION_RATIO).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets Kappa shape indices
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setKappaShapeIndices(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(KAPPA_SHAPE_INDICES).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 3; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets Petitjean number value
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setPetitjeanNumber(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(PETITJEAN_NUMBER).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets spiro atom count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setSpiroAtomCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(SPIRO_ATOM_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the vertex adjacency information (magnitude)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setVAdjMat(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(V_ADJ_MAT).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets weighted path descriptor values
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setWeightedPath(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(WEIGHTED_PATH).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets Zagreb index value
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setZagrebIndex(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(ZAGREB_INDEX).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets carbon types descriptor values
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setCarbonTypesDescriptor(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult result = (IntegerArrayResult) descriptorToCdkObjectMap.get(CARBON_TYPES).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 9; i++) {
            aVector[aStartIndex + i] = result.get(i);
        }
    }

    /**
     * Sets ALogP values (Ghose-Crippen LogP,ALogP squared and molar refractivity)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setALogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            // Create a copy of the molecule with explicit hydrogen atoms
            IAtomContainer moleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);

            // Calculate ALogP using the molecule copy containing explicit hydrogen atoms
            DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(A_LOG_P).calculate(moleculeWithExplicitH).getValue();

            aVector[aStartIndex] = (float) result.get(0);      // ALogP
            aVector[aStartIndex + 1] = (float) result.get(1);  // ALogP squared
            aVector[aStartIndex + 2] = (float) result.get(2);  // Molar Refractivity
        } catch (Exception anException) {
            aVector[aStartIndex] = Float.NaN;
            aVector[aStartIndex + 1] = Float.NaN;
            aVector[aStartIndex + 2] = Float.NaN;
            Descriptor.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setALogP: An exception occurred: " + anException.getMessage(),
                    anException
            );
        }
    }

    /**
     * Sets XLogP value (prediction of logP based on the atom-type method)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     * Note: XLogP requires explicit hydrogens for correct calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setXLogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            // Create a copy of the molecule with explicit hydrogen atoms
            IAtomContainer moleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);

            // Calculate XLogP using the molecule copy containing explicit hydrogen atoms
            aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(X_LOG_P).calculate(moleculeWithExplicitH).getValue()).doubleValue();
        } catch (Exception anException) {
            aVector[aStartIndex] = Float.NaN;
            Descriptor.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setXLogP: An exception occurred: " + anException.getMessage(),
                    anException
            );
        }
    }
    /**
     * Sets JP LogP value (octanol-water partition coefficient based on JPlogP method)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setJPLogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(JP_LOG_P).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets APol value (sum of the atomic polarizabilities)
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAPol(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(A_POL).calculate(anAtomContainer).getValue()).doubleValue();
    }

    // Add new descriptor information here!
    //</editor-fold>

    //<editor-fold desc="Public static molecule preparation methods">
    /**
     * Creates a new molecule with explicit hydrogens from a molecule that has implicit hydrogens.
     * Note: This method requires that the implicit hydrogen counts on atoms are properly set
     *  before calling this method. If implicit hydrogen counts are not set or are all zero,
     *  no explicit hydrogens will be added to the resulting molecule. Checks are NOT performed here.
     *  All necessary checks have already been made in public methods above.
     *
     *
     * @param aMolecule Molecule with implicit hydrogens (IS NOT CHANGED)
     * @return New molecule with explicit hydrogens (Not allowed to be null or empty)
     * @throws Exception Thrown if the molecular structure is invalid or hydrogens cannot be added
     */
    public static IAtomContainer createMoleculeWithExplicitHydrogens(
            IAtomContainer aMolecule
    ) throws Exception {

        try {
            // Create a new empty atom container with the same properties
            IAtomContainer tmpMoleculeWithExplicitH = aMolecule.getBuilder().newInstance(IAtomContainer.class);

            // Copy properties
            for (Object key : aMolecule.getProperties().keySet()) {
                tmpMoleculeWithExplicitH.setProperty(key, aMolecule.getProperty(key));
            }

            // Copy atoms
            for (IAtom tmpAtom : aMolecule.atoms()) {
                IAtom tmpNewAtom = tmpAtom.getBuilder().newInstance(IAtom.class);

                // Copy atom properties
                tmpNewAtom.setSymbol(tmpAtom.getSymbol());
                tmpNewAtom.setAtomicNumber(tmpAtom.getAtomicNumber());
                tmpNewAtom.setMassNumber(tmpAtom.getMassNumber());
                tmpNewAtom.setFormalCharge(tmpAtom.getFormalCharge());
                tmpNewAtom.setImplicitHydrogenCount(tmpAtom.getImplicitHydrogenCount());
                tmpNewAtom.setCharge(tmpAtom.getCharge());

                // Copy atom flags
                if (tmpAtom.isAromatic()) {
                    tmpNewAtom.setIsAromatic(true);
                }
                if (tmpAtom.isInRing()) {
                    tmpNewAtom.setIsInRing(true);
                }

                // Copy atom properties
                for (Object key : tmpAtom.getProperties().keySet()) {
                    tmpNewAtom.setProperty(key, tmpAtom.getProperty(key));
                }

                // Add atom to new container
                tmpMoleculeWithExplicitH.addAtom(tmpNewAtom);
            }

            // Copy bonds
            for (IBond tmpBond : aMolecule.bonds()) {
                IBond tmpNewBond = tmpBond.getBuilder().newInstance(IBond.class);

                // Get atoms for this bond in the new molecule
                IAtom tmpAtom1 = tmpMoleculeWithExplicitH.getAtom(aMolecule.indexOf(tmpBond.getBegin()));
                IAtom tmpAtom2 = tmpMoleculeWithExplicitH.getAtom(aMolecule.indexOf(tmpBond.getEnd()));

                // Set bond properties
                tmpNewBond.setOrder(tmpBond.getOrder());
                tmpNewBond.setAtoms(new IAtom[]{tmpAtom1, tmpAtom2});

                // Copy bond flags
                if (tmpBond.isAromatic()) {
                    tmpNewBond.setIsAromatic(true);
                }
                if (tmpBond.isInRing()) {
                    tmpNewBond.setIsInRing(true);
                }

                // Copy bond properties
                for (Object key : tmpBond.getProperties().keySet()) {
                    tmpNewBond.setProperty(key, tmpBond.getProperty(key));
                }

                // Add bond to new container
                tmpMoleculeWithExplicitH.addBond(tmpNewBond);
            }

            // Now that we have a copy of the original molecule, add explicit hydrogens
            try {
                AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMoleculeWithExplicitH);
            } catch (Exception anException) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.createMoleculeWithExplicitHydrogens: Error converting implicit to explicit hydrogens: " + anException.getMessage(),
                        anException
                );
                throw new Exception("Failed to convert implicit to explicit hydrogens: " + anException.getMessage());
            }

            return tmpMoleculeWithExplicitH;

        } catch (Exception anException) {
            if (!(anException instanceof IllegalArgumentException)) {
                Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.createMoleculeWithExplicitHydrogens: Error creating molecule with explicit hydrogens: " + anException.getMessage(),
                        anException
                );
            }
            throw anException;
        }
    }
    //</editor-fold>

}

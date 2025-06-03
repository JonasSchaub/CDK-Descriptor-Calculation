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

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability;
import org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor;
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
import org.openscience.cdk.silent.SilentChemObjectBuilder;
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
     * Path number: sum of the distances between any two atoms in the molecule.<br>
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.<br>
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     */
    WIENER_NUMBER,
    /**
     * Atom count, counts the number of all atoms in the given molecule.
     */
    ATOM_COUNT,
    /**
     * HBondAcceptorCount, counts hydrogen bond acceptors based on a simplified PHACIR scheme.
     * It includes: Oxygen atoms with formal charge ≤ 0 (excluding: Aromatic ether oxygens and oxygens adjacent to nitrogen).
     * Nitrogen atoms with formal charge ≤ 0 (excluding: Nitrogens adjacent to oxygen).
     */
    H_BOND_ACCEPTOR_COUNT,
    /**
     * HBondDonorCount, counts hydrogen bond donors based on a simplified PHACIR classification.
     * It includes: OH groups where the oxygen has a formal charge ≥ 0 and NH groups where the nitrogen has a formal charge ≥ 0
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
     * BCUT descriptor, calculates Burden matrix modified eigenvalues with different weighting schemes. Returns 6 values:<br>
     * 1. BCUTw-1l, BCUTw-2l ... - nhigh lowest atom weighted BCUTS<br>
     * 2. BCUTw-1h, BCUTw-2h ... - nlow highest atom weighted BCUTS<br>
     * 3. BCUTc-1l, BCUTc-2l ... - nhigh lowest partial charge weighted BCUTS<br>
     * 4. BCUTc-1h, BCUTc-2h ... - nlow highest partial charge weighted BCUTS<br>
     * 5. BCUTp-1l, BCUTp-2l ... - nhigh lowest polarizability weighted BCUTS<br>
     * 6. BCUTp-1h, BCUTp-2h ... - nlow highest polarizability weighted BCUTS<br>
     * Note: No array for one parameter is returned, just the highest and lowest numbers. (Default Parameters: nhigh = 1 and nlow = 1)
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
     * For aromatic bonds counts use AROMATIC_BONDS_COUNT.
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
     * have at least two heavy-atom neighbors. Excluding terminal bonds.
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
     * These indices compare the molecular graph with minimal and maximal molecular graphs. Returns 3 values:<br>
     * 1. Kier1 - First kappa shape index<br>
     * 2. Kier2 - Second kappa shape index<br>
     * 3. Kier3 - Third kappa shape index<br>
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
     * Spiro atom count descriptor, calculates the number of spiro atoms in a molecule.
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
     * Returns 5 values:<br>
     * 1. WTPT1 - molecular ID<br>
     * 2. WTPT2 - molecular ID / number of atoms<br>
     * 3. WTPT3 - sum of path lengths starting from heteroatoms<br>
     * 4. WTPT4 - sum of path lengths starting from oxygens<br>
     * 5. WTPT5 - sum of path lengths starting from nitrogens<br>
     */
    WEIGHTED_PATH,
    /**
     * ZagrebIndex descriptor, calculates the Zagreb index of a molecule.
     * The Zagreb index is the sum of the squares of atom degrees over all heavy atoms,
     * which provides information about the molecular complexity and topological structure.
     */
    ZAGREB_INDEX,
    /**
     * CarbonTypes descriptor, calculates the frequency of occurrence of 9 different types of carbon atoms. Returns 9 values:<br>
     * 1. C1SP1 - triply bound carbon bound to one other carbon<br>
     * 2. C2SP1 - triply bound carbon bound to two other carbons<br>
     * 3. C1SP2 - doubly bound carbon bound to one other carbon<br>
     * 4. C2SP2 - doubly bound carbon bound to two other carbons<br>
     * 5. C3SP2 - doubly bound carbon bound to three other carbons<br>
     * 6. C1SP3 - singly bound carbon bound to one other carbon<br>
     * 7. C2SP3 - singly bound carbon bound to two other carbons<br>
     * 8. C3SP3 - singly bound carbon bound to three other carbons<br>
     * 9. C4SP3 - singly bound carbon bound to four other carbons
     */
    CARBON_TYPES,
    /**
     * ALogP descriptor, calculates Ghose-Crippen LogP values, molar refractivity values
     * and ALogP squared values. Returns 3 values:<br>
     * 1. ALogP (logP value) is the Ghose-Crippen octanol-water partition coefficient.<br>
     * 2. ALogP² is the squared ALogP value.<br>
     * 3. Molar Refractivity (MR) measures the volume occupied by an atom or group of atoms.<br>
     */
    A_LOG_P,
    /**
     * XLogP descriptor, calculates logP based on the atom-type method called XLogP.
     * Requires all hydrogens to be explicit.
     */
    X_LOG_P,
    /**
     * JP_LOG_P descriptor, calculates the octanol-water partition coefficien based on the JPlogP method.
     * Original publication: Junghwan Lee et al. "Estimation of partition coefficients...".
     */
    JP_LOG_P,
    /**
     * APol descriptor, calculates the sum of the atomic polarizabilities (including implicit hydrogens).
     */
    A_POL,
    /**
     * AutocorrelationDescriptorCharge, calculates topological autocorrelation vectors
     * that capture patterns related to charge distribution across the molecular structure.
     * This descriptor correlates atomic partial charges along the molecular topology
     * to characterize charge-related structural patterns in the molecule.
     * Returns 5 values representing charge autocorrelation at different topological distances.
     */
    AUTOCORRELATION_CHARGE,
    /**
     * AutocorrelationDescriptorMass, calculates topological autocorrelation vectors
     * that capture patterns related to atomic mass distribution across the molecular structure.
     * This descriptor correlates atomic masses along the molecular topology
     * to characterize mass-related structural patterns in the molecule.
     * Returns 5 values representing mass autocorrelation at different topological distances.
     */
    AUTOCORRELATION_MASS,
    /**
     * AutocorrelationDescriptorPolarizability, calculates topological autocorrelation vectors
     * that capture patterns related to polarizability distribution across the molecular structure.
     * This descriptor correlates atomic polarizabilities along the molecular topology
     * to characterize polarizability-related structural patterns in the molecule.
     * Returns 5 values representing polarizability autocorrelation at different topological distances.
     * NOTE: Method is not validated in the CDK so not validated in this implementation as well
     */
    AUTOCORRELATION_POLARIZABILITY,
    /**
     * FragmentComplexity descriptor, calculates the complexity of a molecular system.
     * The complexity is defined as [Nilakantan, R. et. al.. Journal of chemical information and modeling. 2006. 46]:
     * C = abs(B^2 - A^2 + A) + H/100
     * where:
     * (C = complexity,
     * A = number of non-hydrogen atoms,
     * B = number of bonds,
     * H = number of heteroatoms,).
     * This provides a measure of structural complexity that correlates with synthetic accessibility.
     */
    FRAGMENT_COMPLEXITY,
    /**
     * ChiChain descriptor, calculates the Kier + Hall chi chain indices of orders 3 through 7.
     * These values characterize a molecular graph based on its chain subgraphs.
     * Returns 10 values:<br>
     * 1. SCH-3 - Simple chain, order 3<br>
     * 2. SCH-4 - Simple chain, order 4<br>
     * 3. SCH-5 - Simple chain, order 5<br>
     * 4. SCH-6 - Simple chain, order 6<br>
     * 5. SCH-7 - Simple chain, order 7<br>
     * 6. VCH-3 - Valence chain, order 3<br>
     * 7. VCH-4 - Valence chain, order 4<br>
     * 8. VCH-5 - Valence chain, order 5<br>
     * 9. VCH-6 - Valence chain, order 6<br>
     * 10. VCH-7 - Valence chain, order 7
     */
    CHI_CHAIN,/**
     * ChiCluster descriptor, calculates Kier + Hall chi cluster indices of orders 3 through 6.
     * These values characterize a molecular graph based on its cluster subgraphs.
     * Returns 8 values:<br>
     * 1. SC-3 - Simple cluster, order 3<br>
     * 2. SC-4 - Simple cluster, order 4<br>
     * 3. SC-5 - Simple cluster, order 5<br>
     * 4. SC-6 - Simple cluster, order 6<br>
     * 5. VC-3 - Valence cluster, order 3<br>
     * 6. VC-4 - Valence cluster, order 4<br>
     * 7. VC-5 - Valence cluster, order 5<br>
     * 8. VC-6 - Valence cluster, order 6
     */
    CHI_CLUSTER,
    /**
     * ChiPathCluster descriptor, calculates Kier + Hall chi path cluster indices of orders 4 through 6.
     * These values characterize a molecular graph based on its path cluster subgraphs.
     * Returns 6 values:<br>
     * 1. SPC-4 - Simple path cluster, order 4<br>
     * 2. SPC-5 - Simple path cluster, order 5<br>
     * 3. SPC-6 - Simple path cluster, order 6<br>
     * 4. VPC-4 - Valence path cluster, order 4<br>
     * 5. VPC-5 - Valence path cluster, order 5<br>
     * 6. VPC-6 - Valence path cluster, order 6
     */
    CHI_PATH_CLUSTER,
    /**
     * ChiPath descriptor, calculates Kier + Hall chi path indices of orders 0 through 7.
     * These values characterize a molecular graph based on its path subgraphs.
     * Returns 16 values:<br>
     * 1.  SP-0 - Simple path, order 0<br>
     * 2.  SP-1 - Simple path, order 1<br>
     * 3.  SP-2 - Simple path, order 2<br>
     * 4.  SP-3 - Simple path, order 3<br>
     * 5.  SP-4 - Simple path, order 4<br>
     * 6.  SP-5 - Simple path, order 5<br>
     * 7.  SP-6 - Simple path, order 6<br>
     * 8.  SP-7 - Simple path, order 7<br>
     * 9.  VP-0 - Valence path, order 0<br>
     * 10. VP-1 - Valence path, order 1<br>
     * 11. VP-2 - Valence path, order 2<br>
     * 12. VP-3 - Valence path, order 3<br>
     * 13. VP-4 - Valence path, order 4<br>
     * 14. VP-5 - Valence path, order 5<br>
     * 15. VP-6 - Valence path, order 6<br>
     * 16. VP-7 - Valence path, order 7<br>
     */
    CHI_PATH,
    /**
     * FractionalPSA descriptor, calculates the ratio of polar surface are to molecular weight.
     * This descriptor provides the polar surface area efficiency, which is the TPSADescriptor value divided by the
     * molecular weight, measure in square Angstroms per Dalton.
     */
    FRACTIONAL_PSA,
    /**
     * LargestPiSystem descriptor, calculates the number of atoms in the largest pi system.
     * This descriptor identifies the largest conjugated pi system within a molecule and
     * returns the count of atoms participating in it.
     */
    LARGEST_PI_SYSTEM,
    /**
     * Descriptor that calculates small ring information.
     * Returns 4 values:<br>
     * 1. nSmallRings - total number of small rings (of size 3 through 9)<br>
     * 2. nAromRings - total number of small aromatic rings<br>
     * 3. nRingBlocks - total number of distinct ring blocks<br>
     * 4. nAromBlocks - total number of aromatically connected components<br>
     * 5. nRings3 - total number of 3-membered rings<br>
     * 6. nRings4 - total number of 4-membered rings<br>
     * 7. nRings5 - total number of 5-membered rings<br>
     * 8. nRings6 - total number of 6-membered rings<br>
     * 9. nRings7 - total number of 7-membered rings<br>
     * 10. nRings8 - total number of 8-membered rings<br>
     * 11. nRings9 - total number of 9-membered rings<br>
     */
    SMALL_RING,
    /**
     * Basic Group Count: Returns the number of basic groups in a molecule.
     */
    BASIC_GROUP_COUNT,
    /**
     * Acidic Group Count: Returns the number of acidic groups in a molecule.
     */
    ACIDIC_GROUP_COUNT,
    /**
     * AminoAcidCount descriptor, calculates the number of each amino acid in a molecule.
     * Returns 20 values, one for each of the 20 standard amino acids:
     * Alanine, Arginine, Asparagine, Aspartic acid, Cysteine, Glutamic acid, Glutamine,
     * Glycine, Histidine, Isoleucine, Leucine, Lysine, Methionine, Phenylalanine,
     * Proline, Serine, Threonine, Tryptophan, Tyrosine, and Valine.
     * This descriptor helps identify and quantify amino acid composition in peptides and proteins.
     */
    AMINO_ACID_COUNT,
    /**
     * Kier-Hall SMARTS descriptor that calculates counts of functional groups and substructures
     * based on the Kier and Hall SMARTS patterns, used for QSAR modeling and molecular characterization.
     * Note: This descriptor provides 79 values representing different molecular fragments.
     */
    KIER_HALL_SMARTS;

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
            // MOLECULAR_WEIGHT has 1 component
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
            BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
            bcutDescriptor.setParameters(new Object[] {1, 1, false}); // nhigh = 1, nlow = 1, checkAromaticity = false
            descriptorToCdkObjectMap.put(BCUT, bcutDescriptor);

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

            // AUTOCORRELATION_CHARGE has 5 components
            descriptorToComponentNumberMap.put(AUTOCORRELATION_CHARGE, 5);
            descriptorToCdkObjectMap.put(AUTOCORRELATION_CHARGE, new AutocorrelationDescriptorCharge());

            // AUTOCORRELATION_MASS has 5 components
            descriptorToComponentNumberMap.put(AUTOCORRELATION_MASS, 5);
            descriptorToCdkObjectMap.put(AUTOCORRELATION_MASS, new AutocorrelationDescriptorMass());

            // AUTOCORRELATION_POLARIZABILITY has 5 components
            descriptorToComponentNumberMap.put(AUTOCORRELATION_POLARIZABILITY, 5);
            descriptorToCdkObjectMap.put(AUTOCORRELATION_POLARIZABILITY, new AutocorrelationDescriptorPolarizability());

            // FRAGMENT_COMPLEXITY has 1 component
            descriptorToComponentNumberMap.put(FRAGMENT_COMPLEXITY, 1);
            descriptorToCdkObjectMap.put(FRAGMENT_COMPLEXITY, new FragmentComplexityDescriptor());

            // CHI_CHAIN has 10 components
            descriptorToComponentNumberMap.put(CHI_CHAIN, 10);
            descriptorToCdkObjectMap.put(CHI_CHAIN, new ChiChainDescriptor());

            // CHI_CLUSTER has 8 components
            descriptorToComponentNumberMap.put(CHI_CLUSTER, 8);
            descriptorToCdkObjectMap.put(CHI_CLUSTER, new ChiClusterDescriptor());

            // CHI_PATH_CLUSTER has 6 components
            descriptorToComponentNumberMap.put(CHI_PATH_CLUSTER, 6);
            descriptorToCdkObjectMap.put(CHI_PATH_CLUSTER, new ChiPathClusterDescriptor());

            // CHI_PATH has 16 components
            descriptorToComponentNumberMap.put(CHI_PATH, 16);
            descriptorToCdkObjectMap.put(CHI_PATH, new ChiPathDescriptor());

            // FRACTIONAL_PSA has 1 component
            descriptorToComponentNumberMap.put(FRACTIONAL_PSA, 1);
            descriptorToCdkObjectMap.put(FRACTIONAL_PSA, new FractionalPSADescriptor());

            // LARGEST_PI_SYSTEM has 1 component
            descriptorToComponentNumberMap.put(LARGEST_PI_SYSTEM, 1);
            LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
            largestPiSystemDescriptor.setParameters(new Object[] {false});
            descriptorToCdkObjectMap.put(LARGEST_PI_SYSTEM, new LargestPiSystemDescriptor());

            // SMALL_RING has 4 components
            descriptorToComponentNumberMap.put(SMALL_RING, 11);
            descriptorToCdkObjectMap.put(SMALL_RING, new SmallRingDescriptor());

            // Basic Group Count has 1 component
            BasicGroupCountDescriptor basicGroupCountDescriptor = new BasicGroupCountDescriptor();
            basicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            descriptorToCdkObjectMap.put(BASIC_GROUP_COUNT, basicGroupCountDescriptor);
            descriptorToComponentNumberMap.put(BASIC_GROUP_COUNT, 1);

            // Acidic Group Count has 1 component
            AcidicGroupCountDescriptor acidicGroupCountDescriptor = new AcidicGroupCountDescriptor();
            acidicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            descriptorToCdkObjectMap.put(ACIDIC_GROUP_COUNT, acidicGroupCountDescriptor);
            descriptorToComponentNumberMap.put(ACIDIC_GROUP_COUNT, 1);

            // AminoAcidCount has 20 components
            descriptorToComponentNumberMap.put(AMINO_ACID_COUNT, 20);
            descriptorToCdkObjectMap.put(AMINO_ACID_COUNT, new AminoAcidCountDescriptor());

            // Kier-Hall SMARTS has 79 components
            descriptorToComponentNumberMap.put(KIER_HALL_SMARTS, 79);
            descriptorToCdkObjectMap.put(KIER_HALL_SMARTS, new KierHallSmartsDescriptor());

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
                case AUTOCORRELATION_CHARGE:
                    setAutocorrelationCharge(anAtomContainer, aVector, aStartIndex);
                    break;
                case AUTOCORRELATION_MASS:
                    setAutocorrelationMass(anAtomContainer, aVector, aStartIndex);
                    break;
                case AUTOCORRELATION_POLARIZABILITY:
                    setAutocorrelationPolarizability(anAtomContainer, aVector, aStartIndex);
                    break;
                case FRAGMENT_COMPLEXITY:
                    setFragmentComplexity(anAtomContainer, aVector, aStartIndex);
                    break;
                case CHI_CHAIN:
                    setChiChain(anAtomContainer, aVector, aStartIndex);
                    break;
                case CHI_CLUSTER:
                    setChiCluster(anAtomContainer, aVector, aStartIndex);
                    break;
                case CHI_PATH_CLUSTER:
                    setChiPathCluster(anAtomContainer, aVector, aStartIndex);
                    break;
                case CHI_PATH:
                    setChiPath(anAtomContainer, aVector, aStartIndex);
                    break;
                case FRACTIONAL_PSA:
                    setFractionalPSA(anAtomContainer, aVector, aStartIndex);
                    break;
                case LARGEST_PI_SYSTEM:
                    setLargestPiSystem(anAtomContainer, aVector, aStartIndex);
                    break;
                case SMALL_RING:
                    setSmallRing(anAtomContainer, aVector, aStartIndex);
                    break;
                case BASIC_GROUP_COUNT:
                    setBasicGroupCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case ACIDIC_GROUP_COUNT:
                    setAcidicGroupCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case AMINO_ACID_COUNT:
                    setAminoAcidCount(anAtomContainer, aVector, aStartIndex);
                    break;
                case KIER_HALL_SMARTS:
                    setKierHallSmarts(anAtomContainer, aVector, aStartIndex);
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
                    BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
                    // Change parameters so that we can use our own setAromaticity method default: checkAromaticity = true
                    bcutDescriptor.setParameters(new Object[] {1, 1, false}); // nhigh = 1, nlow = 1, checkAromaticity = false
                    DoubleArrayResult bcutResult = (DoubleArrayResult) bcutDescriptor.calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        aVector[aStartIndex + i] = (float) bcutResult.get(i);
                    }
                    break;
                case BOND_COUNT_ALL:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new BondCountDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_SPECIFIED:
                    BondCountDescriptor bondCountDesc = new BondCountDescriptor();

                    bondCountDesc.setParameters(new Object[]{"s"});
                    aVector[aStartIndex] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();

                    bondCountDesc.setParameters(new Object[]{"d"});
                    aVector[aStartIndex + 1] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();

                    bondCountDesc.setParameters(new Object[]{"t"});
                    aVector[aStartIndex + 2] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();
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
                    IAtomContainer hybridizationRatioMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    aVector[aStartIndex] = (float) ((DoubleResult) (new HybridizationRatioDescriptor()).calculate(hybridizationRatioMoleculeWithExplicitH).getValue()).doubleValue();
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
                    IAtomContainer aLogPMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    DoubleArrayResult aLogPResult = (DoubleArrayResult) (new ALOGPDescriptor()).calculate(aLogPMoleculeWithExplicitH).getValue();
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
                case AUTOCORRELATION_CHARGE:
                    DoubleArrayResult tmpAutocorrelationChargeResult = (DoubleArrayResult) (new AutocorrelationDescriptorCharge()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpAutocorrelationChargeResult.get(i);
                    }
                    break;
                case AUTOCORRELATION_MASS:
                    DoubleArrayResult tmpAutocorrelationMassResult = (DoubleArrayResult) (new AutocorrelationDescriptorMass()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpAutocorrelationMassResult.get(i);
                    }
                    break;
                case AUTOCORRELATION_POLARIZABILITY:
                    DoubleArrayResult tmpAutocorrelationPolarizabilityResult = (DoubleArrayResult) (new AutocorrelationDescriptorPolarizability()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpAutocorrelationPolarizabilityResult.get(i);
                    }
                    break;
                case FRAGMENT_COMPLEXITY:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new FragmentComplexityDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case CHI_CHAIN:
                    DoubleArrayResult tmpChiChainResult = (DoubleArrayResult) new ChiChainDescriptor().calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 10; i++) {
                        aVector[aStartIndex + i] = (float) tmpChiChainResult.get(i);
                    }
                    break;
                case CHI_CLUSTER:
                    DoubleArrayResult tmpChiClusterResult = (DoubleArrayResult) (new ChiClusterDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 8; i++) {
                        aVector[aStartIndex + i] = (float) tmpChiClusterResult.get(i);
                    }
                    break;
                case CHI_PATH_CLUSTER:
                    DoubleArrayResult chiPathClusterResult = (DoubleArrayResult) new ChiPathClusterDescriptor().calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        aVector[aStartIndex + i] = (float) chiPathClusterResult.get(i);
                    }
                    break;
                case CHI_PATH:
                    DoubleArrayResult tmpArrayResultChiPath = (DoubleArrayResult) (new ChiPathDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 16; i++) {
                        aVector[aStartIndex + i] = (float) tmpArrayResultChiPath.get(i);
                    }
                    break;
                case FRACTIONAL_PSA:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new FractionalPSADescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case LARGEST_PI_SYSTEM:
                    // The largest pi system descriptor is not thread-safe, so we need to copy the molecule otherwise calculation can fail.
                    IAtomContainer tmpMolecule = copyMolecule(anAtomContainer);
                    LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
                    // Change parameters so that we can use our own setAromaticity method
                    largestPiSystemDescriptor.setParameters(new Object[] {false}); // checkAromaticity = false
                    aVector[aStartIndex] = (float) ((IntegerResult) largestPiSystemDescriptor.calculate(tmpMolecule).getValue()).intValue();
                    break;
                case SMALL_RING:
                    IntegerArrayResult smallRingResult = (IntegerArrayResult) (new SmallRingDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 11; i++) {
                        aVector[aStartIndex + i] = (float) smallRingResult.get(i);
                    }
                    break;
                case BASIC_GROUP_COUNT:
                    BasicGroupCountDescriptor basicGroupCountDesc = new BasicGroupCountDescriptor();
                    basicGroupCountDesc.initialise(SilentChemObjectBuilder.getInstance());
                    aVector[aStartIndex] = (float) ((IntegerResult) basicGroupCountDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ACIDIC_GROUP_COUNT:
                    AcidicGroupCountDescriptor acidicGroupCountDesc = new AcidicGroupCountDescriptor();
                    acidicGroupCountDesc.initialise(SilentChemObjectBuilder.getInstance());
                    aVector[aStartIndex] = (float) ((IntegerResult) acidicGroupCountDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case AMINO_ACID_COUNT:
                    IntegerArrayResult aminoAcidCountResult = (IntegerArrayResult) (new AminoAcidCountDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < aminoAcidCountResult.length(); i++) {
                        aVector[aStartIndex + i] = (float) aminoAcidCountResult.get(i);
                    }
                    break;
                case KIER_HALL_SMARTS:
                    IntegerArrayResult kierHallSmartsResult = (IntegerArrayResult) (new KierHallSmartsDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < kierHallSmartsResult.length(); i++) {
                        aVector[aStartIndex + i] = (float) kierHallSmartsResult.get(i);
                    }
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
                    BondCountDescriptor bondCountDesc = (BondCountDescriptor)descriptorToCdkObjectMap.get(BOND_COUNT_SPECIFIED);
                    bondCountDesc.setParameters(new Object[]{"s"});
                    aVector[aStartIndex] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();

                    // Double bonds (d)
                    bondCountDesc.setParameters(new Object[]{"d"});
                    aVector[aStartIndex + 1] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();

                    // Triple bonds (t)
                    bondCountDesc.setParameters(new Object[]{"t"});
                    aVector[aStartIndex + 2] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();
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
                    IAtomContainer hybridizationRatioMoleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(HYBRIDIZATION_RATIO).calculate(hybridizationRatioMoleculeWithExplicitH).getValue()).doubleValue();
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
                    DoubleArrayResult aLogPResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(A_LOG_P).calculate(aLogPMoleculeWithExplicitH).getValue();
                    aVector[aStartIndex] = (float) aLogPResult.get(0);      // ALogP
                    aVector[aStartIndex + 1] = (float) aLogPResult.get(1);  // ALogP squared
                    aVector[aStartIndex + 2] = (float) aLogPResult.get(2);  // Molar Refractivity
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
                case AUTOCORRELATION_CHARGE:
                    DoubleArrayResult tmpAutocorrelationChargeResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(AUTOCORRELATION_CHARGE).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpAutocorrelationChargeResult.get(i);
                    }
                    break;
                case AUTOCORRELATION_MASS:
                    DoubleArrayResult tmpAutocorrelationMassResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(AUTOCORRELATION_MASS).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpAutocorrelationMassResult.get(i);
                    }
                    break;
                case AUTOCORRELATION_POLARIZABILITY:
                    DoubleArrayResult tmpAutocorrelationPolarizabilityResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(AUTOCORRELATION_POLARIZABILITY).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        aVector[aStartIndex + i] = (float) tmpAutocorrelationPolarizabilityResult.get(i);
                    }
                    break;
                case FRAGMENT_COMPLEXITY:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FRAGMENT_COMPLEXITY).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case CHI_CHAIN:
                    DoubleArrayResult tmpChiChainResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_CHAIN).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 10; i++) {
                        aVector[aStartIndex + i] = (float) tmpChiChainResult.get(i);
                    }
                    break;
                case CHI_CLUSTER:
                    DoubleArrayResult tmpChiClusterResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_CLUSTER).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 8; i++) {
                        aVector[aStartIndex + i] = (float) tmpChiClusterResult.get(i);
                    }
                    break;
                case CHI_PATH_CLUSTER:
                    DoubleArrayResult tmpChiPathClusterResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_PATH_CLUSTER).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        aVector[aStartIndex + i] = (float) tmpChiPathClusterResult.get(i);
                    }
                    break;
                case CHI_PATH:
                    DoubleArrayResult tmpArrayResultChiPath = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_PATH).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 16; i++) {
                        aVector[aStartIndex + i] = (float) tmpArrayResultChiPath.get(i);
                    }
                    break;
                case FRACTIONAL_PSA:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FRACTIONAL_PSA).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case LARGEST_PI_SYSTEM:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(LARGEST_PI_SYSTEM).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case SMALL_RING:
                    IntegerArrayResult smallRingResult = (IntegerArrayResult) descriptorToCdkObjectMap.get(SMALL_RING).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 11; i++) {
                        aVector[aStartIndex + i] = (float) smallRingResult.get(i);
                    }
                    break;
                case BASIC_GROUP_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(BASIC_GROUP_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ACIDIC_GROUP_COUNT:
                    aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(ACIDIC_GROUP_COUNT).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case AMINO_ACID_COUNT:
                    IntegerArrayResult aminoAcidCountResult = (IntegerArrayResult) descriptorToCdkObjectMap.get(AMINO_ACID_COUNT).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < aminoAcidCountResult.length(); i++) {
                        aVector[aStartIndex + i] = (float) aminoAcidCountResult.get(i);
                    }
                    break;
                case KIER_HALL_SMARTS:
                    IntegerArrayResult kierHallSmartsResult = (IntegerArrayResult) descriptorToCdkObjectMap.get(KIER_HALL_SMARTS).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < kierHallSmartsResult.length(); i++) {
                        aVector[aStartIndex + i] = (float) kierHallSmartsResult.get(i);
                    }
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
     * Sets molecular weight.
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
     * Sets Wiener number(s).
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
     * Sets atom count.
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
     * Sets hydrogen bond acceptor count.
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
     * Sets hydrogen bond donor count.
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
     * Sets the topological polar surface area (TPSA).
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
     * Sets largest chain size.
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
     * Sets longest aliphatic chain size.
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
     * Sets the Mannhold LogP value (octanol-water partition coefficient).
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
     * Sets BCUT descriptor values.
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
     * Sets bond count.
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
     * Sets specific bond counts (single, double, triple).
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
            BondCountDescriptor bondCountDesc = (BondCountDescriptor)descriptorToCdkObjectMap.get(BOND_COUNT_SPECIFIED);
            bondCountDesc.setParameters(new Object[]{"s"});
            aVector[aStartIndex] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();

            // Double bonds (d)
            bondCountDesc.setParameters(new Object[]{"d"});
            aVector[aStartIndex + 1] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();

            // Triple bonds (t)
            bondCountDesc.setParameters(new Object[]{"t"});
            aVector[aStartIndex + 2] = (float) ((IntegerResult) bondCountDesc.calculate(anAtomContainer).getValue()).intValue();
        } catch (Exception e) {
            for (int i = 0; i < 3; i++) {
                aVector[aStartIndex + i] = Float.NaN;
            }
        }
    }

    /**
     * Sets bond polarizability value.
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
     * Sets the number of aromatic atoms.
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
     * Sets the number of aromatic bonds.
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
     * Sets rotatable bonds count.
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
     * Sets the FMF (Framework Match Fraction) value.
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
     * Sets the FractionalCSP3 value (fraction of sp3 hybridized carbon atoms).
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
     * Sets hybridization ratio (sp3 carbons to sp2 carbons).
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
        try {
            // Create a copy of the molecule with explicit hydrogen atoms
            IAtomContainer moleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);

            // Calculate XLogP using the molecule copy containing explicit hydrogen atoms
            aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(HYBRIDIZATION_RATIO).calculate(moleculeWithExplicitH).getValue()).doubleValue();
        } catch (Exception anException) {
            aVector[aStartIndex] = Float.NaN;
            Descriptor.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setHybridizationRatio: An exception occurred: " + anException.getMessage(),
                    anException
            );
        }
    }

    /**
     * Sets Kappa shape indices.
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
     * Sets Petitjean number value.
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
     * Sets spiro atom count.
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
     * Sets the vertex adjacency information (magnitude).
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
     * Sets weighted path descriptor values.
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
     * Sets Zagreb index value.
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
     * Sets carbon types descriptor values.
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
     * Sets ALogP values (Ghose-Crippen LogP, ALogP squared and molar refractivity).
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
     * Sets XLogP value (prediction of logP based on the atom-type method).
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
     * Sets JP LogP value (octanol-water partition coefficient based on JPlogP method).
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
     * Sets APol value (sum of the atomic polarizabilities).
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

    /**
     * Sets AutocorrelationDescriptorCharge values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAutocorrelationCharge(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(AUTOCORRELATION_CHARGE).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets AutocorrelationDescriptorMass values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAutocorrelationMass(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(AUTOCORRELATION_MASS).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets AutocorrelationDescriptorPolarizability values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAutocorrelationPolarizability(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(AUTOCORRELATION_POLARIZABILITY).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets fragment complexity value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setFragmentComplexity(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FRAGMENT_COMPLEXITY).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets ChiChain descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setChiChain(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_CHAIN).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 10; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets ChiCluster descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setChiCluster(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_CLUSTER).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 8; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets ChiPathCluster descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setChiPathCluster(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_PATH_CLUSTER).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 6; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets ChiPath descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setChiPath(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) descriptorToCdkObjectMap.get(CHI_PATH).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 16; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets FractionalPSA descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     *  @param anAtomContainer Molecule (IS NOT CHANGED)
     *  @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     *  @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setFractionalPSA(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(FRACTIONAL_PSA).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets LargestPiSystem descriptor value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setLargestPiSystem(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(LARGEST_PI_SYSTEM).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets SmallRing descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setSmallRing(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult result = (IntegerArrayResult) descriptorToCdkObjectMap.get(SMALL_RING).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 11; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets basic group count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setBasicGroupCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(BASIC_GROUP_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets acidic group count
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAcidicGroupCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) descriptorToCdkObjectMap.get(ACIDIC_GROUP_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets amino acid count values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setAminoAcidCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult aminoAcidCountResult = (IntegerArrayResult) descriptorToCdkObjectMap.get(AMINO_ACID_COUNT).calculate(anAtomContainer).getValue();
        for (int i = 0; i < aminoAcidCountResult.length(); i++) {
            aVector[aStartIndex + i] = (float) aminoAcidCountResult.get(i);
        }
    }
    /**
     * Sets Kier-Hall SMARTS descriptor values (79 functional group counts).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * Note: Method has to be synchronized due to missing thread-safety of the CDK calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static synchronized void setKierHallSmarts(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult result = (IntegerArrayResult) descriptorToCdkObjectMap.get(KIER_HALL_SMARTS).calculate(anAtomContainer).getValue();
        for (int i = 0; i < result.length(); i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    // Add new descriptor information here!
    //</editor-fold>
    //<editor-fold desc="Public static molecule preparation methods">
    /**
     * Creates a deep copy of the input molecule.
     * Note: This method is used to create a new molecule object
     * without affecting implicit hydrogen atoms.
     *
     * @param aMolecule Source molecule to be copied (NOT MODIFIED)
     * @return New instance of the molecule
     * @throws NullPointerException If the input molecule is null
     * @throws IllegalArgumentException If the input molecule is empty
     * @throws CloneNotSupportedException If the molecule cannot be properly copied
     */
    public static IAtomContainer copyMolecule(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        if (aMolecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Input molecule must not be empty");
        }
        try {
            // Create a new empty atom container with the same properties
            IAtomContainer tmpMoleculeCopy = aMolecule.getBuilder().newInstance(IAtomContainer.class);

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

                // Add atom to new container
                tmpMoleculeCopy.addAtom(tmpNewAtom);
            }

            // Copy bonds
            for (IBond tmpBond : aMolecule.bonds()) {
                IBond tmpNewBond = tmpBond.getBuilder().newInstance(IBond.class);

                // Get atoms for this bond in the new molecule
                IAtom tmpAtom1 = tmpMoleculeCopy.getAtom(aMolecule.indexOf(tmpBond.getBegin()));
                IAtom tmpAtom2 = tmpMoleculeCopy.getAtom(aMolecule.indexOf(tmpBond.getEnd()));

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

                // Add bond to new container
                tmpMoleculeCopy.addBond(tmpNewBond);
            }

            return tmpMoleculeCopy;

        } catch (Exception anException) {
            throw new CloneNotSupportedException("Could not clone molecule: " + anException.getMessage());
        }

    }

    /**
     * Creates a new molecule with all implicit hydrogen atoms made explicit.
     *
     * @param aMolecule Source molecule with implicit hydrogens (NOT MODIFIED)
     * @return New molecule with all hydrogens made explicit
     * @throws NullPointerException If the input molecule is null
     * @throws IllegalArgumentException If the input molecule is empty
     * @throws CloneNotSupportedException If the molecule cannot be properly processed
     */
    public static IAtomContainer createMoleculeWithExplicitHydrogens(
            IAtomContainer aMolecule
    ) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        if (aMolecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Input molecule must not be empty");
        }
        try {
            // Create a deep copy of the molecule first
            IAtomContainer tmpMoleculeCopy = copyMolecule(aMolecule);

            // Add explicit hydrogen atoms
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMoleculeCopy);

            return tmpMoleculeCopy;
        } catch (Exception anException) {
            throw new CloneNotSupportedException("Could not create molecule with explicit hydrogens: " + anException.getMessage());
        }
    }

    /**
     * Uses a specified aromaticity model to modify aMolecule.
     * Note: This method changes the input molecule by applying the specified aromaticity model.
     *
     * @param aMolecule Molecule that will be modified (IS CHANGED)
     * @param anAromaticityModel The aromaticity model that will be used for aromaticity detection
     * @throws NullPointerException If the input molecule or aromaticity model is null
     * @throws IllegalArgumentException If the input molecule is empty
     * @throws Exception If the aromaticity detection fails
     */
    public static void setAromaticity(
            IAtomContainer aMolecule,
            ElectronDonation anAromaticityModel
    ) throws NullPointerException, IllegalArgumentException, Exception {
        if (aMolecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Input molecule must not be empty");
        }
        if (anAromaticityModel == null) {
            throw new NullPointerException("Aromaticity model must not be null");
        }
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);

            // Clears all aromatic flags before applying the aromaticity model.
            Aromaticity.clear(aMolecule);

            Cycles.markRingAtomsAndBonds(aMolecule);
            Aromaticity.apply(anAromaticityModel, aMolecule);
        } catch (Exception anException) {
            throw new Exception("Failed to detect aromaticity: " + anException.getMessage(), anException);
        }
    }

    //</editor-fold>

}

/*
 * CDK-Descriptor-Calculation
 * Copyright (C) 2025 Manuel Schauer, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
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
import org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor;
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
import org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SpiroAtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * Descriptor related calculations based on the CDK for the enrichment of data vectors.
 * <p>
 * Note: There are 5 different "public static boolean setDescriptorsForMolecules...()" methods with different forms of
 * (parallelized) calculation. There is single molecule processing or batch processing for large datastructures.
 * Single and Batch processing can be used with either an IAtomContainer array or with a String Array of SMILES codes.
 * <p>
 * Note: For adding a new descriptor, go to "Add new descriptor information here!"
 *
 * Example usage of the Descriptor class:
 *
 * <pre>{@code
 * // Preprocessing
 * // A: Molecule preprocessing:
 * // 1: Create a molecule
 * IAtomContainer molecule = ...;
 *
 * // 2. IMPORTANT: Perform aromaticity perception before descriptor calculation
 * Descriptor.setAromaticity(molecule);
 *
 * // 3. Create a molecule array
 * IAtomContainer[] molecules = new IAtomContainer[] { molecule };
 *
 * // B: String preprocessing
 * // 1: Create a String array of SMILES codes
 * String[] moleculeSmilesStrings = new String[] {"CCO", "CCC(O)O", ...}
 *
 * // Define parameters
 * // 1: Define descriptor array
 * Descriptor[] descriptors = new Descriptor[] {
 *     Descriptor.MOLECULAR_WEIGHT,
 *     Descriptor.TPSA,
 *     Descriptor.A_LOG_P
 * };
 * // or use a method to get specific descriptors:
 * Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(); // e.g., user-defined selection
 *
 * // 2. Determine total number of components needed
 * int componentCount = Descriptor.getNumberOfComponents(descriptors);
 *
 * // 3. Determine total number of molecules
 * int moleculeCount = molecules.length;
 *
 * // 4. Create data matrix for results with appropriate size
 * float[][] matrix = new float[moleculeCount][componentCount];
 *
 * // 5. Set aStartIndex to 0 for filling the matrix from the beginning
 * int startIndex = 0; // if the matrix is part of a larger data structure, set aStartIndex accordingly
 *
 * // 6. Choose whether to use parallel calculation
 * boolean isParallelCalculation = true; // false for sequential calculation
 *
 * // 7. Create a synchronized List for NaN positions
 * List<int[]> aNanPositionsSynchronized = Collections.synchronizedList(new LinkedList<>());
 *
 * // For Batch Processing:
 * // Define batch size
 * int aBatchSize = 100;
 *
 * // For SMILES Processing
 * // Define anElectronDonation model for aromaticity handling
 * ElectronDonation anElectronDonationModel = Aromaticity.Model.Daylight;
 *
 * // Calculate descriptors (choose one of the following methods)
 *
 * // Option A: High-Performance-Parallel Descriptor Calculation with IAtomContainer array
 * boolean result = Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
 *                 descriptors,
 *                 molecules,
 *                 matrix,
 *                 startIndex,
 *                 isParallelCalculation,
 *                 aNanPositionsSynchronized
 *         );
 * // you can also use Descriptor.setDescriptorsForMoleculesByBatchParallelization(...) -> aBatchSize needs to be specified accordingly
 *
 * // Option B: High-Performance-Parallel Descriptor Calculation with SMILES string array
 * boolean result = Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
 *                descriptors,
 *                moleculeSmilesStrings,
 *                matrix,
 *                startIndex,
 *                anElectronDonationModel,
 *                isParallelCalculation,
 *                aNanPositionsSynchronized
 *        );
 * // you can also use Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(...) -> aBatchSize needs to be specified accordingly
 *
 * }</pre>
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 * @version 1.0.0
 */
public enum Descriptor {
    //<editor-fold desc="Descriptor enumeration and initialization">
    /**
     * Molecular weight, adds up the natural masses (weighted average of all known
     * isotopes of the particular element based on their natural abundances) of every atom
     * in the given molecule, it is NOT the exact mass.
     *
     * @see WeightDescriptor
     */
    MOLECULAR_WEIGHT(true, true, false, false, 1, "Molecular Weight"),
    /**
     * Wiener number, returns Wiener path number and Wiener polarity number.
     * Path number: sum of the distances between any two atoms in the molecule.<br>
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.<br>
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     *
     * @see WienerNumbersDescriptor
     */
    WIENER_NUMBER(true, true, false, false, 2, "Wiener Number"),
    //<editor-fold desc="Basic Bond and Count descriptors">
    /**
     * Atom count, counts the number of all atoms in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT(true, true, false, false, 1, "Atom Count"),
    /**
     * Atom count C, counts the number of all carbon atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_C(true, true, false, false, 1, "Atom Count C"),
    /**
     * Atom count H, counts the number of all hydrogen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_H(true, true, false, false, 1, "Atom Count H"),
    /**
     * Atom count N counts the number of all nitrogen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_N(true, true, false, false, 1, "Atom Count N"),
    /**
     * Atom count O, counts the number of all oxygen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_O(true, true, false, false, 1, "Atom Count O"),
    /**
     * Atom count S, counts the number of all sulfur atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_S(true, true, false, false, 1, "Atom Count S"),
    /**
     * Atom count P, counts the number of all phosphorus atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_P(true, true, false, false, 1, "Atom Count P"),
    /**
     * Atom count F, counts the number of all fluorine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_F(true, true, false, false, 1, "Atom Count F"),
    /**
     * Atom count Br, counts the number of all bromine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_BR(true, true, false, false, 1, "Atom Count Br"),
    /**
     * Atom count Cl, counts the number of all chlorine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_CL(true, true, false, false, 1, "Atom Count Cl"),
    /**
     * Atom count I, counts the number of all iodine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_I(true, true, false, false, 1, "Atom Count I"),
    /**
     * Total bond count, counts the number of all bonds in a molecule, neglecting the order.
     * Double and triple bonds are counted as one bond.
     * Bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_ALL(true, true, false, false, 1, "Bond Count All"),
    /**
     * Bond count single, counts the number of single bonds in a molecule.
     * No bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_SINGLE(true, true, false, false, 1, "Bond Count Single"),
    /**
     * Bond count double, counts the number of double bonds in a molecule.
     * No bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_DOUBLE(true, true, false, false, 1, "Bond Count Double"),
    /**
     * Bond count triple, counts the number of triple bonds in a molecule.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_TRIPLE(true, true, false, false, 1, "Bond Count Triple"),
    /**
     * H bond acceptor count, counts hydrogen bond acceptors based on a simplified PHACIR scheme.
     * It includes: Oxygen atoms with formal charge ≤ 0 (excluding: Aromatic ether oxygens and oxygens adjacent to nitrogen)
     * and nitrogen atoms with formal charge ≤ 0 (excluding: nitrogens adjacent to oxygen).
     *
     * @see HBondAcceptorCountDescriptor
     */
    H_BOND_ACCEPTOR_COUNT(true, true, false, false, 1, "H-Bond Acceptor Count"),
    /**
     * H bond donor count, counts hydrogen bond donors based on a simplified PHACIR classification.
     * It includes: OH groups where the oxygen has a formal charge ≥ 0 and NH groups where the nitrogen has a formal charge ≥ 0.
     *
     * @see HBondDonorCountDescriptor
     */
    H_BOND_DONOR_COUNT(true, true, false, false, 1, "H-Bond Donor Count"),
    /**
     * Aromatic atoms count, counts the number of aromatic atoms in a molecule.
     * Note: Requires that aromatic atoms in the molecule have already been detected and marked.
     *
     * @see AromaticAtomsCountDescriptor
     */
    AROMATIC_ATOMS_COUNT(true, true, false, false, 1, "Aromatic Atoms Count"),
    /**
     * Aromatic bonds count, counts the number of aromatic bonds in a molecule.
     * Note: Requires that aromatic bonds in the molecule have already been detected and marked.
     *
     * @see AromaticBondsCountDescriptor
     */
    AROMATIC_BONDS_COUNT(true, true, false, false, 1, "Aromatic Bonds Count"),
    /**
     * Rotatable bonds count, counts the number of rotatable bonds in a molecule.
     * A rotatable bond is defined as any single non-ring bond, where atoms on both sides
     * have at least two heavy-atom neighbors. Excluding terminal bonds.
     *
     * @see RotatableBondsCountDescriptor
     */
    ROTATABLE_BONDS_COUNT(true, true, false, false, 1, "Rotatable Bonds Count"),
    /**
     * Basic group count, returns the number of basic groups in a molecule.
     *
     * @see BasicGroupCountDescriptor
     */
    BASIC_GROUP_COUNT(true, true, false, false, 1, "Basic Group Count"),
    /**
     * Acidic group count, returns the number of acidic groups in a molecule.
     *
     * @see AcidicGroupCountDescriptor
     */
    ACIDIC_GROUP_COUNT(true, true, false, false, 1, "Acidic Group Count"),
    //</editor-fold>
    /**
     * TPSA descriptor, calculates the topological polar surface area (TPSA) of a molecule.
     * TPSA is the sum of the surface areas of polar atoms (typically oxygen and nitrogen)
     * and their attached hydrogens, based on a topological approximation (2D structure only).
     *
     * @see TPSADescriptor
     */
    TPSA(true, true, false, false, 1, "TPSA"),
    /**
     * Largest chain descriptor, calculates the number of atoms in the longest chain in the molecule.
     * This is a simple topological descriptor that provides a measure of molecular linearity.
     *
     * @see LargestChainDescriptor
     */
    LARGEST_CHAIN(true, true, false, false, 1, "Largest Chain"),
    /**
     * Longest aliphatic chain descriptor, calculates the number of atoms in the longest aliphatic chain.
     * This descriptor provides information about the maximum linear extent of non-aromatic
     * portions of the molecular structure, which relates to molecular shape properties.
     *
     * @see LongestAliphaticChainDescriptor
     */
    LONGEST_ALIPHATIC_CHAIN(true, true, false, false, 1, "Longest Aliphatic Chain"),
    /**
     * BCUT descriptor, calculates Burden matrix modified eigenvalues with different weighting schemes. Returns 6 values:<br>
     * 1. BCUTw-1l, BCUTw-2l ... - nlow lowest atom weighted BCUTS<br>
     * 2. BCUTw-1h, BCUTw-2h ... - nhigh highest atom weighted BCUTS<br>
     * 3. BCUTc-1l, BCUTc-2l ... - nlow lowest partial charge weighted BCUTS<br>
     * 4. BCUTc-1h, BCUTc-2h ... - nhigh highest partial charge weighted BCUTS<br>
     * 5. BCUTp-1l, BCUTp-2l ... - nlow lowest polarizability weighted BCUTS<br>
     * 6. BCUTp-1h, BCUTp-2h ... - nhigh highest polarizability weighted BCUTS<br>
     * Note: No array for one parameter is returned, just the highest and lowest numbers. (Default Parameters: nhigh = 1 and nlow = 1)
     *
     * @see BCUTDescriptor
     */
    BCUT(false, false, false, false, 6, "BCUT"),
    /**
     * Bond polarizability descriptor.
     * The BPolDescriptor calculates the bond polarizability of a molecule.
     * Bond polarizability is a simple sum of polarizability contributions from all bonds, based on bond types and involved atoms.
     * It provides a rough estimate of how easily the electron cloud in a molecule can be distorted,
     * which relates to intermolecular interactions, polarizability and refractive behavior.
     *
     * @see BPolDescriptor
     */
    B_POL(true, true, false, false, 1, "BPol"),
    /**
     * Rule of five descriptor, calculates the number of failures of Lipinski's Rule of Five.
     * The descriptor returns the number of violations (0-4).
     *
     * @see RuleOfFiveDescriptor
     */
    RULE_OF_FIVE(true, true, false, false, 1, "Rule of Five"),
    /**
     * FMF (Framework Match Fraction) descriptor, calculates the ratio of heavy atoms in
     * the framework to the total number of heavy atoms in the molecule.
     * This provides an indication of the proportion of the molecule that is part of the
     * scaffold or core structure, versus the proportion that is in side chains.
     *
     * @see FMFDescriptor
     */
    FMF(true, true, false, false, 1, "FMF"),
    /**
     * Fractional C SP3 descriptor, characterizes the non-flatness of a molecule by calculating
     * the fraction of sp3 hybridized carbon atoms over the total carbon count.
     * This provides information about the three-dimensionality and complexity
     * of a molecule, which relates to drug-likeness properties.
     *
     * @see FractionalCSP3Descriptor
     */
    FRACTIONAL_CSP3(true, true, false, false, 1, "Fractional CSP3"),
    /**
     * Hybridization ratio descriptor, calculates the ratio of sp3 carbons to sp2 carbons.
     * This provides valuable information about the three-dimensionality and flatness
     * of a molecule, which can be useful for predicting drug-like properties and
     * comparing structural characteristics.
     *
     * @see HybridizationRatioDescriptor
     */
    HYBRIDIZATION_RATIO(true, true, false, true, 1, "Hybridization Ratio"),
    /**
     * Kappa shape indices descriptor, calculates Kier and Hall kappa molecular shape indices.
     * These indices compare the molecular graph with minimal and maximal molecular graphs. Returns 3 values:<br>
     * 1. Kier1 - First kappa shape index<br>
     * 2. Kier2 - Second kappa shape index<br>
     * 3. Kier3 - Third kappa shape index<br>
     * Note: Hydrogens are ignored in the calculation.
     *
     * @see KappaShapeIndicesDescriptor
     */
    KAPPA_SHAPE_INDICES(false, true, false, true, 3, "Kappa Shape Indices"),
    /**
     * Petitjean number descriptor, calculates an index characterizing molecular graph topology.
     * This topological descriptor is based on the calculation of the graph eccentricity
     * and provides information about the molecular shape and branching pattern.
     *
     * @see PetitjeanNumberDescriptor
     */
    PETITJEAN_NUMBER(true, true, false, false, 1, "Petitjean Number"),
    /**
     * Spiro atom count descriptor, calculates the number of spiro atoms in a molecule.
     *
     * @see SpiroAtomCountDescriptor
     */
    SPIRO_ATOM_COUNT(true, true, false, false, 1, "Spiro Atom Count"),
    /**
     * VAdjMa descriptor, calculates the Vertex adjacency information (magnitude).
     * This is calculated as 1 + log2 m, where m is the number of heavy-heavy bonds.
     * If m is zero, then zero is returned.
     * This descriptor characterizes molecular complexity in terms of edge connectivity.
     *
     * @see VAdjMaDescriptor
     */
    V_ADJ_MAT(true, true, false, false, 1, "VAdjMa"),
    /**
     * Weighted path descriptor, evaluates the weighted path descriptors for a molecule.
     * Returns 5 values:<br>
     * 1. WTPT1 - molecular ID<br>
     * 2. WTPT2 - molecular ID / number of atoms<br>
     * 3. WTPT3 - sum of path lengths starting from heteroatoms<br>
     * 4. WTPT4 - sum of path lengths starting from oxygens<br>
     * 5. WTPT5 - sum of path lengths starting from nitrogens<br>
     * <p>
     * Note: This descriptor computes all paths which is an NP-hard problem, do not use it for complex molecules.
     * @see WeightedPathDescriptor
     */
    WEIGHTED_PATH(false, true, false, false, 5, "Weighted Path"),
    /**
     * Zagreb index descriptor, calculates the Zagreb index of a molecule.
     * The Zagreb index is the sum of the squares of atom degrees over all heavy atoms,
     * which provides information about the molecular complexity and topological structure.
     *
     * @see ZagrebIndexDescriptor
     */
    ZAGREB_INDEX(true, true, false, false, 1, "Zagreb Index"),
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
     *
     * @see CarbonTypesDescriptor
     */
    CARBON_TYPES(true, true, false, false, 9, "Carbon Types"),
    //<editor-fold desc="LogP descriptors">
    /**
     * ALogP descriptor, calculates Ghose-Crippen LogP values, molar refractivity values
     * and ALogP squared values. Returns 3 values:<br>
     * 1. ALogP (logP value) is the Ghose-Crippen octanol-water partition coefficient.<br>
     * 2. ALogP² is the squared ALogP value.<br>
     * 3. Molar Refractivity (MR) measures the volume occupied by an atom or group of atoms.<br>
     *
     * @see ALOGPDescriptor
     */
    A_LOG_P(true, true, false, true, 3, "ALogP"),
    /**
     * XLogP descriptor, calculates logP based on the atom-type method called XLogP.
     * Requires all hydrogen's to be explicit.
     *
     * @see XLogPDescriptor
     */
    X_LOG_P(true, true, false, true, 1, "XLogP"),
    /**
     * JPlogP descriptor, calculates the octanol-water partition coefficient based on an atom contribution model.
     *
     * @see JPlogPDescriptor
     */
    JP_LOG_P(true, false, false, false, 1, "JPlogP"),
    /**
     * Mannhold LogP descriptor, calculates the octanol-water partition coefficient (logP) using the Mannhold method.
     * LogP describes the hydrophilicity or lipophilicity of a compound and is crucial for
     * predicting solubility, permeability, and bioavailability.
     *
     * @see MannholdLogPDescriptor
     */
    MANNHOLD_LOGP(true, true, false, false, 1, "Mannhold LogP"),
    //</editor-fold>
    /**
     * APol descriptor, calculates the sum of the atomic polarizabilities (including implicit hydrogens).
     *
     * @see APolDescriptor
     */
    A_POL(true, true, false, false, 1, "APol"),
    //<editor-fold desc="Autocorrelation descriptors">
    /**
     * Autocorrelation charge descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to charge distribution across the molecular structure.
     * This descriptor correlates atomic partial charges along the molecular topology
     * to characterize charge-related structural patterns in the molecule.
     * Returns 5 values representing charge autocorrelation at different topological distances.
     *
     * @see AutocorrelationDescriptorCharge
     */
    AUTOCORRELATION_CHARGE(true, false, false, false, 5, "Autocorrelation Charge"),
    /**
     * Autocorrelation mass descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to atomic mass distribution across the molecular structure.
     * This descriptor correlates atomic masses along the molecular topology
     * to characterize mass-related structural patterns in the molecule.
     * Returns 5 values representing mass autocorrelation at different topological distances.
     *
     * @see AutocorrelationDescriptorMass
     */
    AUTOCORRELATION_MASS(true, true, false, false, 5, "Autocorrelation Mass"),
    /**
     * Autocorrelation polarizability descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to polarizability distribution across the molecular structure.
     * This descriptor correlates atomic polarizabilities along the molecular topology
     * to characterize polarizability-related structural patterns in the molecule.
     * Returns 5 values representing polarizability autocorrelation at different topological distances.
     * NOTE: Method is not validated in the CDK so not validated in this implementation as well
     *
     * @see AutocorrelationDescriptorPolarizability
     */
    AUTOCORRELATION_POLARIZABILITY(true, true, false, false, 5, "Autocorrelation Polarizability"),
    //</editor-fold>
    /**
     * Fragment complexity descriptor, calculates the complexity of a molecular system.
     * The complexity is defined as [Nilakantan, R. et al. Journal of chemical information and modeling. 2006. 46]:
     * C = abs(B^2 - A^2 + A) + H/100
     * where:
     * (C = complexity,
     * A = number of non-hydrogen atoms,
     * B = number of bonds,
     * H = number of heteroatoms,).
     * This provides a measure of structural complexity that correlates with synthetic accessibility.
     *
     * @see FragmentComplexityDescriptor
     */
    FRAGMENT_COMPLEXITY(true, true, false, false, 1, "Fragment Complexity"),
    //<editor-fold desc="CHI descriptors">
    /**
     * Chi chain descriptor, calculates the Kier + Hall chi chain indices of orders 3 through 7.
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
     *
     * @see ChiChainDescriptor
     */
    CHI_CHAIN(false, true, false, false, 10, "Chi Chain"),
    /**
     * Chi cluster descriptor, calculates Kier + Hall chi cluster indices of orders 3 through 6.
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
     *
     * @see ChiClusterDescriptor
     */
    CHI_CLUSTER(false, true, false, false, 8, "Chi Cluster"),
    /**
     * Chi path cluster descriptor, calculates Kier + Hall chi path cluster indices of orders 4 through 6.
     * These values characterize a molecular graph based on its path cluster subgraphs.
     * Returns 6 values:<br>
     * 1. SPC-4 - Simple path cluster, order 4<br>
     * 2. SPC-5 - Simple path cluster, order 5<br>
     * 3. SPC-6 - Simple path cluster, order 6<br>
     * 4. VPC-4 - Valence path cluster, order 4<br>
     * 5. VPC-5 - Valence path cluster, order 5<br>
     * 6. VPC-6 - Valence path cluster, order 6
     *
     * @see ChiPathClusterDescriptor
     */
    CHI_PATH_CLUSTER(false, true, false, false, 6, "Chi Path Cluster"),
    /**
     * Chi path descriptor, calculates Kier + Hall chi path indices of orders 0 through 7.
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
     *
     * @see ChiPathDescriptor
     */
    CHI_PATH(false, true, false, false, 16, "Chi Path"),
    //</editor-fold>
    /**
     * Fractional PSA descriptor, calculates the ratio of polar surface area to molecular weight.
     * This descriptor provides the polar surface area efficiency, which is the TPSADescriptor value divided by the
     * molecular weight, measured in square Angstroms per Dalton.
     *
     * @see FractionalPSADescriptor
     */
    FRACTIONAL_PSA(true, true, false, false, 1, "Fractional PSA"),
    /**
     * Largest pi system descriptor, calculates the number of atoms in the largest pi system.
     * This descriptor identifies the largest conjugated pi system within a molecule and
     * returns the count of atoms participating in it.
     *
     * @see LargestPiSystemDescriptor
     */
    LARGEST_PI_SYSTEM(true, true, false, false, 1, "Largest Pi System"),
    /**
     * Descriptor that calculates small ring information.
     * Returns 11 values:<br>
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
     *
     * @see SmallRingDescriptor
     */
    SMALL_RING(true, true, false, false, 11, "Small Ring"),
    /**
     * Amino acid count descriptor, calculates the number of each amino acid in a molecule.
     * Returns 20 values, one for each of the 20 standard amino acids:
     * Alanine, Arginine, Asparagine, Aspartic acid, Cysteine, Glutamic acid, Glutamine,
     * Glycine, Histidine, Isoleucine, Leucine, Lysine, Methionine, Phenylalanine,
     * Proline, Serine, Threonine, Tryptophan, Tyrosine, and Valine.
     * This descriptor helps identify and quantify amino acid composition in peptides and proteins.
     *
     * @see AminoAcidCountDescriptor
     */
    AMINO_ACID_COUNT(false, true, false, false, 20, "Amino Acid Count"),
    /**
     * Kier-Hall SMARTS descriptor that calculates counts of functional groups and substructures
     * based on the Kier and Hall SMARTS patterns, used for QSAR modeling and molecular characterization.
     * Note: This descriptor provides 79 values representing different molecular fragments.
     *
     * @see KierHallSmartsDescriptor
     */
    KIER_HALL_SMARTS(true, true, false, false, 79, "Kier Hall SMARTS"),
    /**
     * Eccentric connectivity index descriptor, calculates a topological descriptor that combines
     * distance and adjacency information.
     * It is defined as the sum of the products of eccentricity and vertex degree for each atom.
     * This index provides information about the distribution of atoms in the molecular structure
     * and helps characterize molecular complexity, branching, and overall shape.
     *
     * @see EccentricConnectivityIndexDescriptor
     */
    ECCENTRIC_CONNECTIVITY_INDEX(true, true, false, false, 1, "Eccentric Connectivity Index"),
    /**
     * MDE descriptor, calculates molecular distance edge descriptors for carbon, oxygen and nitrogen atoms.
     * These descriptors encode information about the connectivity and distance of atoms of specific types
     * and hybridization states in the molecular graph. Returns 19 values representing various molecular
     * distance edge counts between different types of carbon, oxygen, and nitrogen atoms:<br>
     * MDEC-11<br>
     * MDEC-12<br>
     * MDEC-13<br>
     * MDEC-14<br>
     * MDEC-22<br>
     * MDEC-23<br>
     * MDEC-24<br>
     * MDEC-33<br>
     * MDEC-34<br>
     * MDEC-44<br>
     * MDEO-11<br>
     * MDEO-12<br>
     * MDEO-22<br>
     * MDEN-11<br>
     * MDEN-12<br>
     * MDEN-13<br>
     * MDEN-22<br>
     * MDEN-23<br>
     * MDEN-33<br>
     *
     * @see MDEDescriptor
     */
    MDE(true, true, false, false, 19, "MDE"),
    /**
     * VABC descriptor, calculates the volume descriptor using the van der Waals volume calculation approach.
     * This descriptor estimates molecular volume based on atom contributions, considering bond types
     * and atomic properties, providing insights into molecular size and steric properties.
     *
     * @see VABCDescriptor
     */
    VABC(true, true, false, false, 1, "VABC"),
    /**
     * PubChem fingerprinter, generates a 881-bit binary fingerprint based on PubChem's substructure keys.
     * This fingerprint encodes the presence or absence of specific substructural features
     * and is useful for similarity searching and chemical space analysis.
     *
     * @see PubchemFingerprinter
     */
    PUBCHEM_FINGERPRINTER(false, true, true, false, 881, "PubChem Fingerprinter"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 0.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_0(true, true, true, false, 1024, "Circular Fingerprinter ECFP-0"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 0.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_0(true, true, true, false, 1024, "Circular Fingerprinter FCFP-0"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 2.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_2(true, true, true, false, 1024, "Circular Fingerprinter ECFP-2"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 2.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_2(true, true, true, false, 1024, "Circular Fingerprinter FCFP-2"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 4.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_4(true, true, true, false, 1024, "Circular Fingerprinter ECFP-4"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 4.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_4(true, true, true, false, 1024, "Circular Fingerprinter FCFP-4"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 6.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_6(true, true, true, false, 1024, "Circular Fingerprinter ECFP-6"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 6.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_6(true, true, true, false, 1024, "Circular Fingerprinter FCFP-6"),
    /**
     * MACCS fingerprinter, generates a 166-bit binary fingerprint based on the MACCS structural keys.
     *
     * @see MACCSFingerprinter
     */
    MACCS_FINGERPRINTER(true, true, true, false, 166, "MACCS Fingerprinter");

    // Add new descriptor information here!
    //</editor-fold>

    // <editor-fold desc="Descriptor information and constructor">
    /**
     * Indicates whether this descriptor is quickly calculable.
     * Fast descriptors have lower computational complexity and can be calculated efficiently,
     * while slow descriptors require more intensive calculations (e.g., CHI descriptors).
     */
    private final boolean isFast;

    /**
     * Indicates whether this descriptor is safe to calculate.
     * Safe descriptors are reliable and validated, while unsafe descriptors may produce
     * inconsistent results or NaN values under certain conditions.
     */
    private final boolean isSafe;

    /**
     * Indicates whether the enum entry is a fingerprint or not.
     */
    private final boolean isFingerprint;

    /**
     * Indicates whether this descriptor needs explicit hydrogens.
     */
    private final boolean needsExplicitHydrogens;

    /**
     * The number of components calculated by this descriptor.
     */
    private final int descriptorComponentNumber;

    /**
     * The human-readable name of this descriptor for output purposes.
     */
    private final String name;

    /**
     * Constructs a Descriptor with the specified speed, safety characteristics, and component number.
     *
     * @param isFast true if the descriptor is quickly calculable, false if it requires
     *               more intensive computation
     * @param isSafe true if the descriptor is safe and reliable, false if it may produce
     *               inconsistent results or NaN values
     * @param isFingerprint true if the enum entry is a fingerprint, false otherwise
     * @param needsExplicitHydrogens true if the descriptor requires explicit hydrogens, false otherwise
     * @param descriptorComponentNumber the number of components calculated by this descriptor
     * @param name the human-readable name of this descriptor for output purposes
     */
    Descriptor(boolean isFast, boolean isSafe, boolean isFingerprint, boolean needsExplicitHydrogens, int descriptorComponentNumber, String name) {
        this.isFast = isFast;
        this.isSafe = isSafe;
        this.isFingerprint = isFingerprint;
        this.needsExplicitHydrogens = needsExplicitHydrogens;
        this.descriptorComponentNumber = descriptorComponentNumber;
        this.name = name;
    }

    /**
     * Returns whether this descriptor can be calculated quickly. This classification is based on a performance benchmark
     * of 10,000 molecules: descriptors computed in under 1 second for 10,000 molecules are considered fast.
     *
     * @return true if the descriptor has low computational complexity and can be
     *         calculated efficiently, false if it requires intensive computation
     */
    public boolean isFast() {
        return isFast;
    }

    /**
     * Returns whether this descriptor is safe to calculate.
     *
     * @return true if the descriptor is reliable and validated, false if it may
     *         produce inconsistent results or NaN values under certain conditions
     */
    public boolean isSafe() {
        return isSafe;
    }

    /**
     * Returns whether this descriptor calculates a fingerprint.
     *
     * @return true if the enum entry calculates a fingerprint, false otherwise
     */
    public boolean isFingerprint() {
        return isFingerprint;
    }

    /**
     * Returns whether this descriptor needs explicit hydrogens.
     *
     * @return true if the descriptor requires explicit hydrogens, false otherwise
     */
    public boolean needsExplicitHydrogens() {
        return needsExplicitHydrogens;
    }

    /**
     * Returns the number of components calculated by this descriptor.
     *
     * @return the number of components calculated by this descriptor
     */
    public int getDescriptorComponentNumber() {
        return descriptorComponentNumber;
    }

    /**
     * Returns the human-readable name of this descriptor.
     *
     * @return the human-readable name of this descriptor
     */
    public String getName() {
        return name;
    }

    /**
     * Calculates descriptor or fingerprint values for a molecule and stores them in the result vector.
     * <p>
     * This is the core calculation method that dispatches to either fingerprint or molecular descriptor
     * calculation based on the descriptor type. It retrieves pre-initialized CDK descriptor instances
     * from a shared pool to avoid repeated instantiation overhead.
     * <p>
     * <b>Important:</b> This method does NOT perform input validation. All necessary checks
     * must be performed by the calling public methods before invoking this method.
     * <p>
     * <b>Thread Safety:</b> For fingerprints, this method uses a blocking queue pool
     * ({@link #calculateFingerprintFromPool}) because fingerprinter instances are not thread-safe.
     * For molecular descriptors, shared instances from {@link #descriptorToCdkObjectMap} are used.
     * <p>
     * <b>Result Handling:</b> The method handles four CDK result types:
     * <ul>
     *   <li>{@link DoubleResult} - Single double value (e.g., molecular weight)</li>
     *   <li>{@link IntegerResult} - Single integer value (e.g., atom count)</li>
     *   <li>{@link DoubleArrayResult} - Array of doubles (e.g., BCUT returns 6 values)</li>
     *   <li>{@link IntegerArrayResult} - Array of integers (e.g., amino acid counts)</li>
     * </ul>
     *
     * @param anAtomContainer the molecule to calculate descriptors for (IS NOT CHANGED);
     *                        must have aromaticity already perceived if required by the descriptor
     * @param aVector         the result vector to store calculated values (MAY BE CHANGED);
     *                        values are stored starting at aStartIndex
     * @param aStartIndex     the starting index in aVector where results should be written;
     *                        subsequent values are written to aStartIndex + 1, aStartIndex + 2, etc.
     * @throws CDKException if the calculation fails for any reason, wrapping the original exception
     *                      with information about which descriptor failed
     */
    private void calculate(IAtomContainer anAtomContainer, float[] aVector, int aStartIndex) throws CDKException {
        try {
            if (this.isFingerprint()) {
                this.calculateFingerprintFromPool(anAtomContainer, aVector, aStartIndex);
            } else {
                IMolecularDescriptor cdkDescriptor = descriptorToCdkObjectMap.get(this);
                if (cdkDescriptor != null) {

                    IDescriptorResult result = cdkDescriptor.calculate(anAtomContainer).getValue();

                    if (result instanceof DoubleResult) {
                        DoubleResult doubleResult = (DoubleResult) result;
                        aVector[aStartIndex] = (float) doubleResult.doubleValue();
                    } else if (result instanceof IntegerResult) {
                        IntegerResult integerResult = (IntegerResult) result;
                        aVector[aStartIndex] = (float) integerResult.intValue();
                    } else if (result instanceof DoubleArrayResult) {
                        DoubleArrayResult arrayResult = (DoubleArrayResult) result;
                        for (int i = 0; i < this.descriptorComponentNumber; i++) {
                            aVector[aStartIndex + i] = (float) arrayResult.get(i);
                        }
                    } else if (result instanceof IntegerArrayResult) {
                        IntegerArrayResult arrayResult = (IntegerArrayResult) result;
                        for (int i = 0; i < this.descriptorComponentNumber; i++) {
                            aVector[aStartIndex + i] = (float) arrayResult.get(i);
                        }
                    }
                }
            }
        } catch (Exception e) {
            throw new CDKException("Calculation failed for " + this.name(), e);
        }
    }
    //</editor-fold>

    //<editor-fold desc="Private static final LOGGER">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(Descriptor.class.getName());
    //</editor-fold>

    // <editor-fold desc="SMILES parser">
    /**
     * SMILES Parser for SMILES String batch processing
     */
    private static final SmilesParser SMILES_PARSER = new SmilesParser(SilentChemObjectBuilder.getInstance());
    //</editor-fold>

    // <editor-fold desc="CDK descriptor and fingerprinter mappings and static initializer block">
    /**
     * EnumMap that maps a descriptor to an instance of its CDK descriptor class
     */
    private static final EnumMap<Descriptor, IMolecularDescriptor> descriptorToCdkObjectMap = new EnumMap<>(Descriptor.class);
    /**
     * Pool size for fingerprinter instances. Default value is 4 which should be sufficient for regular users.
     * Can be changed via setFingerprintPoolSize() method.
     */
    private static int fingerprintPoolSize = 4;
    /**
     * EnumMap that maps a fingerprint descriptor to its pool of IFingerprinter instances.
     * Pool size is configurable via setFingerprintPoolSize() method.
     * Uses BlockingQueue to ensure thread-safe access to fingerprinter instances.
     */
    private static final EnumMap<Descriptor, BlockingQueue<IFingerprinter>> fingerprintPoolMap = new EnumMap<>(Descriptor.class);

    /*
     * Static initializer block to populate the descriptorToCdkObjectMap.
     */
    static {
        try {
            // MOLECULAR_WEIGHT
            descriptorToCdkObjectMap.put(MOLECULAR_WEIGHT, new WeightDescriptor());

            // WIENER_NUMBER
            descriptorToCdkObjectMap.put(WIENER_NUMBER, new WienerNumbersDescriptor());

            // ATOM_COUNT
            descriptorToCdkObjectMap.put(ATOM_COUNT, new AtomCountDescriptor());

            // ATOM_COUNT_C
            AtomCountDescriptor atomCountCDescriptor = new AtomCountDescriptor();
            atomCountCDescriptor.setParameters(new Object[] {"C"}); // set carbon as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_C, atomCountCDescriptor);

            // ATOM_COUNT_H
            AtomCountDescriptor atomCountHDescriptor = new AtomCountDescriptor();
            atomCountHDescriptor.setParameters(new Object[] {"H"}); // set hydrogen as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_H, atomCountHDescriptor);

            // ATOM_COUNT_N
            AtomCountDescriptor atomCountNDescriptor = new AtomCountDescriptor();
            atomCountNDescriptor.setParameters(new Object[] {"N"}); // set nitrogen as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_N, atomCountNDescriptor);

            // ATOM_COUNT_O
            AtomCountDescriptor atomCountODescriptor = new AtomCountDescriptor();
            atomCountODescriptor.setParameters(new Object[] {"O"}); // set oxygen as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_O, atomCountODescriptor);

            // ATOM_COUNT_S
            AtomCountDescriptor atomCountSDescriptor = new AtomCountDescriptor();
            atomCountSDescriptor.setParameters(new Object[] {"S"}); // set sulfur as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_S, atomCountSDescriptor);

            // ATOM_COUNT_P
            AtomCountDescriptor atomCountPDescriptor = new AtomCountDescriptor();
            atomCountPDescriptor.setParameters(new Object[] {"P"}); // set phosphorus as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_P, atomCountPDescriptor);

            // ATOM_COUNT_F
            AtomCountDescriptor atomCountFDescriptor = new AtomCountDescriptor();
            atomCountFDescriptor.setParameters(new Object[] {"F"}); // set fluorine as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_F, atomCountFDescriptor);

            // ATOM_COUNT_BR
            AtomCountDescriptor atomCountBrDescriptor = new AtomCountDescriptor();
            atomCountBrDescriptor.setParameters(new Object[] {"Br"}); // set bromine as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_BR, atomCountBrDescriptor);

            // ATOM_COUNT_CL
            AtomCountDescriptor atomCountClDescriptor = new AtomCountDescriptor();
            atomCountClDescriptor.setParameters(new Object[] {"Cl"}); // set chlorine as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_CL, atomCountClDescriptor);

            // ATOM_COUNT_I
            AtomCountDescriptor atomCountIDescriptor = new AtomCountDescriptor();
            atomCountIDescriptor.setParameters(new Object[] {"I"}); // set iodine as the atom to count
            descriptorToCdkObjectMap.put(ATOM_COUNT_I, atomCountIDescriptor);

            // H_BOND_ACCEPTOR_COUNT
            descriptorToCdkObjectMap.put(H_BOND_ACCEPTOR_COUNT, new HBondAcceptorCountDescriptor());

            // H_BOND_DONOR_COUNT
            descriptorToCdkObjectMap.put(H_BOND_DONOR_COUNT, new HBondDonorCountDescriptor());

            // TPSA
            descriptorToCdkObjectMap.put(TPSA, new TPSADescriptor());

            // LARGEST_CHAIN
            descriptorToCdkObjectMap.put(LARGEST_CHAIN, new LargestChainDescriptor());

            // LONGEST_ALIPHATIC_CHAIN
            descriptorToCdkObjectMap.put(LONGEST_ALIPHATIC_CHAIN, new LongestAliphaticChainDescriptor());

            // MANNHOLD_LOGP
            descriptorToCdkObjectMap.put(MANNHOLD_LOGP, new MannholdLogPDescriptor());

            // BCUT
            BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
            bcutDescriptor.setParameters(new Object[] {1, 1, false}); // nhigh = 1, nlow = 1, checkAromaticity = false
            descriptorToCdkObjectMap.put(BCUT, bcutDescriptor);

            // BOND_COUNT_ALL
            descriptorToCdkObjectMap.put(BOND_COUNT_ALL, new BondCountDescriptor());

            // BOND_COUNT_SINGLE
            BondCountDescriptor bondCountSingleDescriptor = new BondCountDescriptor();
            bondCountSingleDescriptor.setParameters(new Object[]{"s"});
            descriptorToCdkObjectMap.put(BOND_COUNT_SINGLE, bondCountSingleDescriptor);

            // BOND_COUNT_DOUBLE
            BondCountDescriptor bondCountDoubleDescriptor = new BondCountDescriptor();
            bondCountDoubleDescriptor.setParameters(new Object[]{"d"});
            descriptorToCdkObjectMap.put(BOND_COUNT_DOUBLE, bondCountDoubleDescriptor);

            // BOND_COUNT_TRIPLE
            BondCountDescriptor bondCountTripleDescriptor = new BondCountDescriptor();
            bondCountTripleDescriptor.setParameters(new Object[]{"t"});
            descriptorToCdkObjectMap.put(BOND_COUNT_TRIPLE, bondCountTripleDescriptor);

            // B_POL
            descriptorToCdkObjectMap.put(B_POL, new BPolDescriptor());

            // RULE_OF_FIVE
            descriptorToCdkObjectMap.put(RULE_OF_FIVE, new RuleOfFiveDescriptor());

            // AROMATIC_ATOMS_COUNT
            descriptorToCdkObjectMap.put(AROMATIC_ATOMS_COUNT, new AromaticAtomsCountDescriptor());

            // AROMATIC_BONDS_COUNT
            descriptorToCdkObjectMap.put(AROMATIC_BONDS_COUNT, new AromaticBondsCountDescriptor());

            // ROTATABLE_BONDS_COUNT
            descriptorToCdkObjectMap.put(ROTATABLE_BONDS_COUNT, new RotatableBondsCountDescriptor());

            // FMF
            descriptorToCdkObjectMap.put(FMF, new FMFDescriptor());

            // FRACTIONAL_CSP3
            descriptorToCdkObjectMap.put(FRACTIONAL_CSP3, new FractionalCSP3Descriptor());

            // HYBRIDIZATION_RATIO
            descriptorToCdkObjectMap.put(HYBRIDIZATION_RATIO, new HybridizationRatioDescriptor());

            // KAPPA_SHAPE_INDICES
            descriptorToCdkObjectMap.put(KAPPA_SHAPE_INDICES, new KappaShapeIndicesDescriptor());

            // PETITJEAN_NUMBER
            descriptorToCdkObjectMap.put(PETITJEAN_NUMBER, new PetitjeanNumberDescriptor());

            // SPIRO_ATOM_COUNT
            descriptorToCdkObjectMap.put(SPIRO_ATOM_COUNT, new SpiroAtomCountDescriptor());

            // V_ADJ_MAT
            descriptorToCdkObjectMap.put(V_ADJ_MAT, new VAdjMaDescriptor());

            // WEIGHTED_PATH
            descriptorToCdkObjectMap.put(WEIGHTED_PATH, new WeightedPathDescriptor());

            // ZAGREB_INDEX
            descriptorToCdkObjectMap.put(ZAGREB_INDEX, new ZagrebIndexDescriptor());

            // CARBON_TYPES
            descriptorToCdkObjectMap.put(CARBON_TYPES, new CarbonTypesDescriptor());

            // A_LOG_P
            descriptorToCdkObjectMap.put(A_LOG_P, new ALOGPDescriptor());

            // X_LOG_P
            descriptorToCdkObjectMap.put(X_LOG_P, new XLogPDescriptor());

            // JP_LOG_P
            descriptorToCdkObjectMap.put(JP_LOG_P, new JPlogPDescriptor());

            // A_POL
            descriptorToCdkObjectMap.put(A_POL, new APolDescriptor());

            // AUTOCORRELATION_CHARGE
            descriptorToCdkObjectMap.put(AUTOCORRELATION_CHARGE, new AutocorrelationDescriptorCharge());

            // AUTOCORRELATION_MASS
            descriptorToCdkObjectMap.put(AUTOCORRELATION_MASS, new AutocorrelationDescriptorMass());

            // AUTOCORRELATION_POLARIZABILITY
            descriptorToCdkObjectMap.put(AUTOCORRELATION_POLARIZABILITY, new AutocorrelationDescriptorPolarizability());

            // FRAGMENT_COMPLEXITY
            descriptorToCdkObjectMap.put(FRAGMENT_COMPLEXITY, new FragmentComplexityDescriptor());

            // CHI_CHAIN
            descriptorToCdkObjectMap.put(CHI_CHAIN, new ChiChainDescriptor());

            // CHI_CLUSTER
            descriptorToCdkObjectMap.put(CHI_CLUSTER, new ChiClusterDescriptor());

            // CHI_PATH_CLUSTER
            descriptorToCdkObjectMap.put(CHI_PATH_CLUSTER, new ChiPathClusterDescriptor());

            // CHI_PATH
            descriptorToCdkObjectMap.put(CHI_PATH, new ChiPathDescriptor());

            // FRACTIONAL_PSA
            descriptorToCdkObjectMap.put(FRACTIONAL_PSA, new FractionalPSADescriptor());

            // LARGEST_PI_SYSTEM
            LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
            largestPiSystemDescriptor.setParameters(new Object[] {false}); //do not check aromaticity again
            descriptorToCdkObjectMap.put(LARGEST_PI_SYSTEM, largestPiSystemDescriptor);

            // SMALL_RING
            descriptorToCdkObjectMap.put(SMALL_RING, new SmallRingDescriptor());

            // Basic Group Count
            BasicGroupCountDescriptor basicGroupCountDescriptor = new BasicGroupCountDescriptor();
            basicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            descriptorToCdkObjectMap.put(BASIC_GROUP_COUNT, basicGroupCountDescriptor);

            // Acidic Group Count
            AcidicGroupCountDescriptor acidicGroupCountDescriptor = new AcidicGroupCountDescriptor();
            acidicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            descriptorToCdkObjectMap.put(ACIDIC_GROUP_COUNT, acidicGroupCountDescriptor);

            // AminoAcidCount
            descriptorToCdkObjectMap.put(AMINO_ACID_COUNT, new AminoAcidCountDescriptor());

            // Kier-Hall SMARTS
            descriptorToCdkObjectMap.put(KIER_HALL_SMARTS, new KierHallSmartsDescriptor());

            // ECCENTRIC_CONNECTIVITY_INDEX
            descriptorToCdkObjectMap.put(ECCENTRIC_CONNECTIVITY_INDEX, new EccentricConnectivityIndexDescriptor());

            // MDE
            descriptorToCdkObjectMap.put(MDE, new MDEDescriptor());

            // VABC
            descriptorToCdkObjectMap.put(VABC, new VABCDescriptor());

            // Initialize fingerprint pools with configurable pool size (default: 4)
            initializeFingerprintPools();

            // Add new descriptor information here!

        } catch (Exception anException) {
            throw new RuntimeException("Failed to initialize descriptors, this should never happen. ", anException);
        }

    }
    //</editor-fold>

    //<editor-fold desc="Private static methods for pool management">
    /**
     * Initializes all fingerprint pools with the current fingerprintPoolSize.
     * This method creates new BlockingQueues for each fingerprint type and populates them
     * with fingerprinter instances. Thread-safe for concurrent access.
     */
    private static synchronized void initializeFingerprintPools() {
        // Clear existing pools
        fingerprintPoolMap.clear();

        // PUBCHEM_FINGERPRINTER Pool
        BlockingQueue<IFingerprinter> pubchemPool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            pubchemPool.offer(new PubchemFingerprinter(SilentChemObjectBuilder.getInstance()));
        }
        fingerprintPoolMap.put(PUBCHEM_FINGERPRINTER, pubchemPool);

        // CIRCULAR_FINGERPRINTER_ECFP_0 Pool
        BlockingQueue<IFingerprinter> ecfp0Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            ecfp0Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP0, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_ECFP_0, ecfp0Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_0 Pool
        BlockingQueue<IFingerprinter> fcfp0Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            fcfp0Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP0, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_FCFP_0, fcfp0Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_2 Pool
        BlockingQueue<IFingerprinter> ecfp2Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            ecfp2Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP2, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_ECFP_2, ecfp2Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_2 Pool
        BlockingQueue<IFingerprinter> fcfp2Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            fcfp2Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP2, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_FCFP_2, fcfp2Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_4 Pool
        BlockingQueue<IFingerprinter> ecfp4Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            ecfp4Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_ECFP_4, ecfp4Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_4 Pool
        BlockingQueue<IFingerprinter> fcfp4Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            fcfp4Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP4, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_FCFP_4, fcfp4Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_6 Pool
        BlockingQueue<IFingerprinter> ecfp6Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            ecfp6Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_ECFP_6, ecfp6Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_6 Pool
        BlockingQueue<IFingerprinter> fcfp6Pool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            fcfp6Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP6, 1024));
        }
        fingerprintPoolMap.put(CIRCULAR_FINGERPRINTER_FCFP_6, fcfp6Pool);

        // MACCS_FINGERPRINTER Pool
        BlockingQueue<IFingerprinter> maccsPool = new LinkedBlockingQueue<>(fingerprintPoolSize);
        for (int i = 0; i < fingerprintPoolSize; i++) {
            maccsPool.offer(new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()));
        }
        fingerprintPoolMap.put(MACCS_FINGERPRINTER, maccsPool);
    }
    //</editor-fold>

    //<editor-fold desc="Public static methods">
    /**
     * Returns all available descriptors.
     *
     * @return All available descriptors
     */
    public static Descriptor[] getAllDescriptors(){
        return values();
    }

    /**
     * Returns all available fingerprint descriptors.
     *
     * @return All available fingerprint descriptors
     */
    public static Descriptor[] getAllFingerprints() {
        return Arrays.stream(values())
                .filter(Descriptor::isFingerprint)
                .toArray(Descriptor[]::new);
    }

    /**
     * Returns specified available descriptors.
     * Note: Fast, safe, and non-fingerprint descriptors are automatically included by default.
     *
     * @param isSlowlyCalculableDescriptorInclusion True: Slowly calculable descriptors are included in the result, false: Otherwise.
     * @param isUnsafeDescriptorInclusion True: Unsafe descriptors are included in the result, false: Unsafe descriptors
     *                                     are excluded from the result, this does not mean no NaN's can be produced.
     * @param isFingerprintAsDescriptorInclusion True: Fingerprint descriptors are included in the result, false: Fingerprint descriptors are excluded.
     * @return Specified descriptors
     */
    public static Descriptor[] getSpecifiedDescriptors(
            boolean isSlowlyCalculableDescriptorInclusion,
            boolean isUnsafeDescriptorInclusion,
            boolean isFingerprintAsDescriptorInclusion
    ) {
        // Initialize ArrayList with maximum possible capacity to avoid internal resizing during element addition
        List<Descriptor> result = new ArrayList<>(values().length);

        for (Descriptor descriptor : values()) {
            boolean includeDescriptor = true;

            // Exclude slow descriptors if not requested
            if (!descriptor.isFast() && !isSlowlyCalculableDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude unsafe descriptors if not requested
            if (!descriptor.isSafe() && !isUnsafeDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude fingerprint descriptors if not requested
            if (descriptor.isFingerprint() && !isFingerprintAsDescriptorInclusion) {
                includeDescriptor = false;
            }

            if (includeDescriptor) { // Only fast, safe, and non-fingerprint descriptors remain if all flags are false
                result.add(descriptor);
            }
        }

        return result.toArray(new Descriptor[0]);
    }

    /**
     * Returns sum of number of calculated components of an array of defined descriptors.
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
                tmpTotalNumberOfComponents += tmpDescriptor.getDescriptorComponentNumber();
            }
            return tmpTotalNumberOfComponents;
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.getNumberOfComponents: Number of components could not be evaluated: This should never happen.", anException
            );
            throw new Exception("Descriptor.getNumberOfComponents: An exception occurred: This should never happen.", anException);
        }
    }

    /**
     * Returns the descriptor and component index for a given position in a descriptor array.
     * <p>
     * Example: Given descriptors [MOLECULAR_WEIGHT (1 component), BCUT (6 components), WIENER_NUMBER (2 components)]:
     * <ul>
     *   <li>Index 0 → MOLECULAR_WEIGHT, component 0</li>
     *   <li>Index 1 → BCUT, component 0</li>
     *   <li>Index 6 → BCUT, component 5</li>
     *   <li>Index 7 → WIENER_NUMBER, component 0</li>
     *   <li>Index 8 → WIENER_NUMBER, component 1</li>
     * </ul>
     *
     * @param aDescriptors Array of descriptors (must not be null or empty)
     * @param anIndex The index in the descriptor array (must be &gt;= 0 and &lt; total component count)
     * @return A two-element int array: [descriptor index in aDescriptors, component index within that descriptor]
     * @throws IllegalArgumentException if aDescriptors is null/empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public static int[] getDescriptorAndComponentIndex(Descriptor[] aDescriptors, int anIndex)
            throws IllegalArgumentException {
        // Checks
        if (aDescriptors == null || aDescriptors.length == 0) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: aDescriptors must not be null or empty."
            );
        }
        if (anIndex < 0) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: anIndex must be >= 0."
            );
        }

        // Calculate cumulative component counts and find the target descriptor
        int cumulativeCount = 0;
        for (int i = 0; i < aDescriptors.length; i++) {
            if (aDescriptors[i] == null) {
                throw new IllegalArgumentException(
                        "Descriptor.getDescriptorAndComponentIndex: aDescriptors contains null element at index " + i
                );
            }

            int componentCount = aDescriptors[i].getDescriptorComponentNumber();
            int nextCumulativeCount = cumulativeCount + componentCount;

            // Check if anIndex falls within this descriptor's range
            if (anIndex < nextCumulativeCount) {
                int componentIndex = anIndex - cumulativeCount;
                return new int[] { i, componentIndex };
            }

            cumulativeCount = nextCumulativeCount;
        }

        // If we reach here, anIndex is out of bounds
        throw new IllegalArgumentException(
                "Descriptor.getDescriptorAndComponentIndex: anIndex (" + anIndex +
                        ") exceeds total component count (" + cumulativeCount + ")."
        );
    }

    /**
     * Returns the descriptor and its component name for a given position in a descriptor array.
     * <p>
     * This is a convenience method that extends {@link #getDescriptorAndComponentIndex} by also
     * providing human-readable names for the descriptor and component.
     *
     * @param aDescriptors Array of descriptors (must not be null or empty)
     * @param anIndex The index in the descriptor array (must be &gt;= 0 and &lt; total component count)
     * @return A String array: [descriptor name, component index as string]
     * @throws IllegalArgumentException if aDescriptors is null/empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public static String[] getDescriptorAndComponentInfo(Descriptor[] aDescriptors, int anIndex)
            throws IllegalArgumentException {
        int[] indices = getDescriptorAndComponentIndex(aDescriptors, anIndex);
        Descriptor descriptor = aDescriptors[indices[0]];

        return new String[] {
                descriptor.getName(),
                String.valueOf(indices[1])
        };
    }

    /**
     * Sets the pool size for fingerprinter instances and reinitialized all fingerprint pools.
     * The default pool size is 4, which should be sufficient for regular users.
     * Increasing the pool size can improve performance in highly parallel environments
     * but will increase memory usage.
     * <p>
     * This method is thread-safe and will block until all pools are reinitialized.
     * Any fingerprinter instances currently in use will be returned to the old pools
     * and will eventually be garbage collected.
     *
     * @param aPoolSize The new pool size for fingerprinter instances (must be greater than 0)
     * @throws IllegalArgumentException if aPoolSize is less than or equal to 0
     */
    public static synchronized void setFingerprintPoolSize(int aPoolSize) throws IllegalArgumentException {
        if (aPoolSize <= 0) {
            throw new IllegalArgumentException("Descriptor.setFingerprintPoolSize: aPoolSize must be greater than 0.");
        }
        fingerprintPoolSize = aPoolSize;
        initializeFingerprintPools();
        Descriptor.LOGGER.log(
            Level.INFO,
            () -> "Fingerprint pool size changed to: " + aPoolSize + ". All fingerprint pools have been reinitialized."
        );
    }

    /**
     * Returns the current pool size for fingerprinter instances.
     *
     * @return Current fingerprint pool size
     */
    public static synchronized int getFingerprintPoolSize() {
        return fingerprintPoolSize;
    }

    /**
     * Uses a specified aromaticity model to modify aMolecule.
     * The method performs a complete workflow of:
     * <ol>
     *     <li>Suppresses explicit hydrogen's.</li>
     *     <li>Perceives atom types and configures atoms.</li>
     *     <li>Clears existing aromaticity flags.</li>
     *     <li>Detects rings.</li>
     *     <li>Applies the specified {@link ElectronDonation} model.</li>
     *     <li>Adds implicit hydrogen's.</li>
     * </ol>
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
        //<editor-fold desc="Checks">
        if (aMolecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Input molecule must not be empty");
        }
        if (anAromaticityModel == null) {
            throw new NullPointerException("Aromaticity model must not be null");
        }
        //</editor-fold>
        try {
            AtomContainerManipulator.suppressHydrogens(aMolecule);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
            // Clears all aromatic flags before applying the aromaticity model.
            Aromaticity.clear(aMolecule);
            Cycles.markRingAtomsAndBonds(aMolecule);
            Aromaticity.apply(anAromaticityModel, aMolecule);
            CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(aMolecule.getBuilder());
            adder.addImplicitHydrogens(aMolecule);
        } catch (Exception anException) {
            throw new Exception("Failed to detect aromaticity: " + anException.getMessage(), anException);
        }
    }

    /**
     * Calculates a single descriptor for the given molecule and returns the result as a float array.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for
     * single descriptor calculations. If the calculation fails, all values in the result array
     * will be set to {@link Float#NaN}.
     * <p>
     * <b>Important:</b> The molecule must have aromaticity already perceived if required by the descriptor.
     * Use {@link #setAromaticity(IAtomContainer, ElectronDonation)} before calling this method if needed.
     *
     * @param aDescriptor The descriptor to calculate (must not be null)
     * @param aMolecule The molecule to calculate the descriptor for (must not be null or empty);
     *                  aromaticity must be perceived beforehand if required
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link #getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     * @throws NullPointerException if aDescriptor or aMolecule is null
     * @throws IllegalArgumentException if aMolecule is empty
     */
    public static float[] calculateDescriptor(Descriptor aDescriptor, IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException {
        // Checks
        if (aDescriptor == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: aDescriptor must not be null.");
        }
        if (aMolecule == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: aMolecule must not be null.");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Descriptor.calculateDescriptor: aMolecule must not be empty.");
        }

        float[] result = new float[aDescriptor.descriptorComponentNumber];
        try {
            if (aDescriptor.needsExplicitHydrogens){
                IAtomContainer tmpMoleculeWithExplicitHydrogens = createMoleculeWithExplicitHydrogens(aMolecule);
                aDescriptor.calculate(tmpMoleculeWithExplicitHydrogens, result, 0);
            } else {
                aDescriptor.calculate(aMolecule, result, 0);
            }
        } catch (Exception e) {
            // Fill with NaN on failure and log the error
            Arrays.fill(result, Float.NaN);
            LOGGER.log(Level.WARNING, () ->
                    "Failed to calculate descriptor " + aDescriptor.getName() + ": " + e.getMessage()
            );
        }
        return result;
    }

    /**
     * Calculates a single descriptor for the given molecule with automatic aromaticity perception.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for single descriptor calculations.
     * The method automatically perceives aromaticity using the specified electron donation model before calculating the descriptor.
     * If aromaticity perception or calculation fails, all values in the result array will be set to {@link Float#NaN}.
     *
     * @param aDescriptor The descriptor to calculate (must not be null)
     * @param aMolecule The molecule to calculate the descriptor for (must not be null or empty)
     * @param anElectronDonation The electron donation model to use for aromaticity perception (must not be null)
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link #getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     * @throws NullPointerException if aDescriptor, aMolecule, or anElectronDonation is null
     * @throws IllegalArgumentException if aMolecule is empty
     */
    public static float[] calculateDescriptor(Descriptor aDescriptor, IAtomContainer aMolecule,
                                              ElectronDonation anElectronDonation)
            throws NullPointerException, IllegalArgumentException {
        // Checks
        if (aDescriptor == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: aDescriptor must not be null.");
        }
        if (aMolecule == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: aMolecule must not be null.");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Descriptor.calculateDescriptor: aMolecule must not be empty.");
        }
        if (anElectronDonation == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: anElectronDonation must not be null.");
        }

        float[] result = new float[aDescriptor.descriptorComponentNumber];
        try {
            IAtomContainer tmpMolecule = copyMolecule(aMolecule); // Create a copy to avoid modifying the original molecule
            setAromaticity(tmpMolecule, anElectronDonation);
            if (aDescriptor.needsExplicitHydrogens){
                IAtomContainer tmpMoleculeWithExplicitHydrogens = createMoleculeWithExplicitHydrogens(tmpMolecule);
                aDescriptor.calculate(tmpMoleculeWithExplicitHydrogens, result, 0);
            } else {
                aDescriptor.calculate(tmpMolecule, result, 0);
            }
        } catch (Exception e) {
            Arrays.fill(result, Float.NaN);
            LOGGER.log(Level.WARNING, () ->
                    "Failed to calculate descriptor " + aDescriptor.getName() + ": " + e.getMessage()
            );
        }
        return result;
    }

    /**
     * Calculates a single descriptor for the given SMILES string with automatic parsing and aromaticity perception.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for single descriptor calculations.
     * The method automatically parses the SMILES string into a molecule,
     * perceives aromaticity using the specified electron donation model, and then calculates the descriptor.
     * If parsing, aromaticity perception, or calculation fails, all values in the result array
     * will be set to {@link Float#NaN}.
     *
     * @param aDescriptor The descriptor to calculate (must not be null)
     * @param aSmilesString The SMILES string representing the molecule (must not be null or empty)
     * @param anElectronDonation The electron donation model to use for aromaticity perception (must not be null)
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link #getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     * @throws NullPointerException if aDescriptor, aSmilesString, or anElectronDonation is null
     * @throws IllegalArgumentException if aSmilesString is blank
     */
    public static float[] calculateDescriptor(Descriptor aDescriptor, String aSmilesString,
                                              ElectronDonation anElectronDonation)
            throws NullPointerException, IllegalArgumentException {
        // Checks
        if (aDescriptor == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: aDescriptor must not be null.");
        }
        if (aSmilesString == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: aSmilesString must not be null.");
        }
        if (aSmilesString.isBlank()) {
            throw new IllegalArgumentException("Descriptor.calculateDescriptor: aSmilesString must not be blank.");
        }
        if (anElectronDonation == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: anElectronDonation must not be null.");
        }

        float[] result = new float[aDescriptor.descriptorComponentNumber];
        try {
            IAtomContainer tmpMolecule = SMILES_PARSER.parseSmiles(aSmilesString);
            setAromaticity(tmpMolecule, anElectronDonation);
            if (aDescriptor.needsExplicitHydrogens){
                    IAtomContainer tmpMoleculeWithExplicitHydrogens = createMoleculeWithExplicitHydrogens(tmpMolecule);
                    aDescriptor.calculate(tmpMoleculeWithExplicitHydrogens, result, 0);
            } else {
                    aDescriptor.calculate(tmpMolecule, result, 0);
            }
        } catch (Exception e) {
            Arrays.fill(result, Float.NaN);
            LOGGER.log(Level.WARNING, () ->
                    "Failed to calculate descriptor " + aDescriptor.getName() + " for SMILES '" + aSmilesString + "': " + e.getMessage()
            );
        }
        return result;
    }


    /**
     * Sets calculated descriptor components in vectors (rows) of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecule batches. If parallel computation is used, the atom
     * container batches (molecules) are distributed onto parallel thread, one for each molecule batch, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via Descriptor.setFingerprintPoolSize().
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *                              anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                              column to start filling with descriptors
     * @param aBatchSize Number of molecules to process in each batch
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                              IMPORTANT: For parallel calculations (anIsParallelCalculation=true), this must be thread-safe.
     *                              Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculesByBatchParallelization(
            Descriptor[] aDescriptors,
            IAtomContainer[] anAtomContainerArray,
            float[][] aMatrix,
            int aStartIndex,
            int aBatchSize,
            boolean anIsParallelCalculation,
            List<int[]> aNanPositions
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        final String methodName = "setDescriptorsForMoleculesByMoleculeBatchParallelization";
        Descriptor.validateDescriptors(methodName, aDescriptors);
        Descriptor.validateAtomContainerArray(methodName, anAtomContainerArray);
        Descriptor.validateMatrix(methodName, aDescriptors, aMatrix, anAtomContainerArray, aStartIndex);
        Descriptor.validateNanPositions(methodName, aNanPositions);
        Descriptor.validateBatchSize(methodName, aBatchSize);
        //</editor-fold>

        try {
            int tmpNumberOfMolecules = anAtomContainerArray.length;
            int tmpNumberOfBatches = (int) Math.ceil((double) tmpNumberOfMolecules / aBatchSize);

            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += aDescriptors[i].getDescriptorComponentNumber();
            }

            AtomicBoolean tmpHasNaN = new AtomicBoolean(false);

            if (anIsParallelCalculation) {
                try {
                    IntStream.range(0, tmpNumberOfBatches).parallel().forEach(batchIndex -> {
                        int tmpBatchStart = batchIndex * aBatchSize;
                        int tmpBatchEnd = Math.min(tmpBatchStart + aBatchSize, tmpNumberOfMolecules);

                        for (int i = tmpBatchStart; i < tmpBatchEnd; i++) {
                            try {
                                boolean tmpSuccess = Descriptor.setDescriptorsForSingleMolecule(
                                        aDescriptors,
                                        anAtomContainerArray[i],
                                        aMatrix[i],
                                        tmpStartIndices,
                                        i,
                                        aNanPositions
                                );
                                if (!tmpSuccess) {
                                    tmpHasNaN.set(true);
                                }
                            } catch (Exception anException) {
                                tmpHasNaN.set(true);
                                Descriptor.LOGGER.log(
                                        Level.WARNING,
                                        "Descriptor.setDescriptorsForMoleculesByBatchParallelization: Exception in batch " + batchIndex + ", molecule index: " + i,
                                        anException
                                );
                            }
                        }
                    });
                } catch (Exception anException) {
                    Descriptor.LOGGER.log(
                            Level.WARNING,
                            "Descriptor.setDescriptorsForMoleculesByBatchParallelization: Global exception occurred in descriptor calculation.",
                            anException
                    );
                    return false;
                }
            } else {
                for (int i = 0; i < tmpNumberOfMolecules; i++) {
                    if (!Descriptor.setDescriptorsForSingleMolecule(aDescriptors, anAtomContainerArray[i], aMatrix[i], tmpStartIndices, i, aNanPositions)) {
                        tmpHasNaN.set(true);
                    }
                }
            }

            return !tmpHasNaN.get();
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForMoleculesByBatchParallelization: An exception occurred.", anException);
        }
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecule batches. If parallel computation is used, the smiles string
     * batches (molecules) are distributed onto parallel thread, one for each molecule batch, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via Descriptor.setFingerprintPoolSize().
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param aMoleculeSmilesStringArray Array of molecule smiles strings. Note: aMoleculeSmilesStringArray[i] corresponds to aMatrix[i] data
     *                              vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *                                      aMoleculeSmilesStringArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                                      column to start filling with descriptors
     * @param aBatchSize Number of molecules to process in each batch
     * @param anElectronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     *                                      Aromaticity.Model.Daylight is used as default.
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                                     IMPORTANT: For parallel calculations (anIsParallelCalculation=true), this must be thread-safe.
     *                                     Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
            Descriptor[] aDescriptors,
            String[] aMoleculeSmilesStringArray,
            float[][] aMatrix,
            int aStartIndex,
            int aBatchSize,
            ElectronDonation anElectronDonationModel,
            boolean anIsParallelCalculation,
            List<int[]> aNanPositions
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        final String methodName = "setDescriptorsForMoleculeBySmilesStringsBatchParallelization";
        Descriptor.validateDescriptors(methodName, aDescriptors);
        Descriptor.validateSmilesStringArray(methodName, aMoleculeSmilesStringArray);
        Descriptor.validateMatrix(methodName, aDescriptors, aMatrix, aMoleculeSmilesStringArray, aStartIndex);
        Descriptor.validateNanPositions(methodName, aNanPositions);
        Descriptor.validateBatchSize(methodName, aBatchSize);
        //</editor-fold>

        try {
            int tmpNumberOfMolecules = aMoleculeSmilesStringArray.length;
            int tmpNumberOfBatches = (int) Math.ceil((double) tmpNumberOfMolecules / aBatchSize);

            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += aDescriptors[i].getDescriptorComponentNumber();
            }

            AtomicBoolean tmpHasNaN = new AtomicBoolean(false);

            if (anIsParallelCalculation) {
                try {
                    IntStream.range(0, tmpNumberOfBatches).parallel().forEach(batchIndex -> {
                        int tmpBatchStart = batchIndex * aBatchSize;
                        int tmpBatchEnd = Math.min(tmpBatchStart + aBatchSize, tmpNumberOfMolecules);

                        for (int i = tmpBatchStart; i < tmpBatchEnd; i++) {
                            try {
                                boolean tmpSuccess = Descriptor.setDescriptorsForSingleMoleculeSmilesString(
                                        aDescriptors,
                                        aMoleculeSmilesStringArray[i],
                                        aMatrix[i],
                                        tmpStartIndices,
                                        i,
                                        anElectronDonationModel,
                                        aNanPositions
                                );
                                if (!tmpSuccess) {
                                    tmpHasNaN.set(true);
                                }
                            } catch (Exception anException) {
                                tmpHasNaN.set(true);
                                Descriptor.LOGGER.log(
                                        Level.WARNING,
                                        "Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Exception in batch " + batchIndex + ", molecule index: " + i,
                                        anException
                                );
                            }
                        }
                    });
                } catch (Exception anException) {
                    Descriptor.LOGGER.log(
                            Level.WARNING,
                            "Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Global exception occurred in descriptor calculation.",
                            anException
                    );
                    return false;
                }
            } else {
                for (int i = 0; i < tmpNumberOfMolecules; i++) {
                    if (!Descriptor.setDescriptorsForSingleMoleculeSmilesString(aDescriptors, aMoleculeSmilesStringArray[i], aMatrix[i], tmpStartIndices, i, anElectronDonationModel, aNanPositions)) {
                        tmpHasNaN.set(true);
                    }
                }
            }

            return !tmpHasNaN.get();
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: An exception occurred.", anException);
        }
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecules. If parallel computation is used, the atom
     * containers (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via Descriptor.setFingerprintPoolSize().
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *                              anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                              column to start filling with descriptors
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                              IMPORTANT: For parallel calculations (anIsParallelCalculation=true), this must be thread-safe.
     *                              Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculesByMoleculeParallelization(
            Descriptor[] aDescriptors,
            IAtomContainer[] anAtomContainerArray,
            float[][] aMatrix,
            int aStartIndex,
            boolean anIsParallelCalculation,
            List<int[]> aNanPositions
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        final String methodName = "setDescriptorsForMoleculesByMoleculeParallelization";
        Descriptor.validateDescriptors(methodName, aDescriptors);
        Descriptor.validateAtomContainerArray(methodName, anAtomContainerArray);
        Descriptor.validateMatrix(methodName, aDescriptors, aMatrix, anAtomContainerArray, aStartIndex);
        Descriptor.validateNanPositions(methodName, aNanPositions);
        //</editor-fold>

        try {
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += aDescriptors[i].getDescriptorComponentNumber();
            }

            AtomicBoolean tmpHasNaN = new AtomicBoolean(false);

            if (anIsParallelCalculation) {
                try {
                    // Advise by Oracle: Parallel streams should use the common Fork-join pool
                    IntStream.range(0, anAtomContainerArray.length).parallel().forEach(
                            i ->
                            {
                                try {
                                    boolean tmpSuccess = Descriptor.setDescriptorsForSingleMolecule(
                                                    aDescriptors,
                                                    anAtomContainerArray[i],
                                                    aMatrix[i],
                                                    tmpStartIndices,
                                                    i,
                                                    aNanPositions
                                            );
                                    if (!tmpSuccess) {
                                        tmpHasNaN.set(true);
                                    }
                                } catch (Exception anException) {
                                    tmpHasNaN.set(true);
                                    Descriptor.LOGGER.log(
                                            Level.WARNING,
                                            "Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: One descriptor calculation caused an exception, molecule index: "
                                                    + i
                                                    + ".",
                                            anException
                                    );
                                }
                            }
                    );
                    return !tmpHasNaN.get();
                } catch (Exception anException) {
                    Descriptor.LOGGER.log(
                            Level.WARNING,
                            "Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: Global exception occurred in descriptor calculation: ",
                            anException
                    );
                    return false;
                }
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < anAtomContainerArray.length; i++) {
                    if (!Descriptor.setDescriptorsForSingleMolecule(aDescriptors, anAtomContainerArray[i], aMatrix[i], tmpStartIndices, i, aNanPositions)) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: An exception occurred.", anException);
        }
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecules. If parallel computation is used, the smiles strings
     * (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via Descriptor.setFingerprintPoolSize().
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param aMoleculeSmilesStringArray Array of molecule smiles strings. Note: aMoleculeSmilesStringArray[i] corresponds to aMatrix[i] data
     *                              vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *                              aMoleculeSmilesStringArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                              column to start filling with descriptors
     * @param anElectronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     *                               Aromaticity.Model.Daylight is used as default.
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                              IMPORTANT: For parallel calculations (anIsParallelCalculation=true), this must be thread-safe.
     *                               Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setDescriptorsForMoleculesBySmilesStringParallelization(
            Descriptor[] aDescriptors,
            String[] aMoleculeSmilesStringArray,
            float[][] aMatrix,
            int aStartIndex,
            ElectronDonation anElectronDonationModel,
            boolean anIsParallelCalculation,
            List<int[]> aNanPositions
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        final String methodName = "setDescriptorsForMoleculesBySmilesStringParallelization";
        Descriptor.validateDescriptors(methodName, aDescriptors);
        Descriptor.validateSmilesStringArray(methodName, aMoleculeSmilesStringArray);
        Descriptor.validateMatrix(methodName, aDescriptors, aMatrix, aMoleculeSmilesStringArray, aStartIndex);
        Descriptor.validateNanPositions(methodName, aNanPositions);
        //</editor-fold>

        try {
            int tmpNumberOfMolecules = aMoleculeSmilesStringArray.length;
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += aDescriptors[i].getDescriptorComponentNumber();
            }
            AtomicBoolean tmpHasNaN = new AtomicBoolean(false);
            if (anIsParallelCalculation) {
                try {
                    // Advise by Oracle: Parallel streams should use the common Fork-join pool
                    IntStream.range(0, aMoleculeSmilesStringArray.length).parallel().forEach(
                            i ->
                            {
                                try {
                                    boolean tmpSuccess = Descriptor.setDescriptorsForSingleMoleculeSmilesString(
                                            aDescriptors,
                                            aMoleculeSmilesStringArray[i],
                                            aMatrix[i],
                                            tmpStartIndices,
                                            i,
                                            anElectronDonationModel,
                                            aNanPositions
                                    );
                                    if (!tmpSuccess) {
                                        tmpHasNaN.set(true);
                                    }
                                } catch (Exception anException) {
                                    tmpHasNaN.set(true);
                                    Descriptor.LOGGER.log(
                                            Level.WARNING,
                                            "Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: Exception in molecule index: " + i,
                                            anException
                                    );
                                }

                });
            } catch (Exception anException) {
                Descriptor.LOGGER.log(
                        Level.WARNING,
                        "Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: Global exception occurred in descriptor calculation.",
                        anException
                );
                return false;
            }
        } else {
            for (int i = 0; i < tmpNumberOfMolecules; i++) {
                if (!Descriptor.setDescriptorsForSingleMoleculeSmilesString(aDescriptors, aMoleculeSmilesStringArray[i], aMatrix[i], tmpStartIndices, i, anElectronDonationModel, aNanPositions)) {
                    tmpHasNaN.set(true);
                }
            }
        }

        return !tmpHasNaN.get();
    } catch (Exception anException) {
        throw new Exception("Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: An exception occurred.", anException);
        }
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecules. If parallel computation is used, the atom
     * containers (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: Uses a new descriptor instance for EVERY descriptor calculation which slows down the calculation.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector. (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *               anAtomContainerArray[i], i.e. the molecules define the rows of the matrix. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *      *             column to start filling with descriptors
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (anIsParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
        public static boolean setDescriptorsForMoleculesByMoleculeParallelizationNew(
            Descriptor[] aDescriptors,
            IAtomContainer[] anAtomContainerArray,
            float[][] aMatrix,
            int aStartIndex,
            boolean anIsParallelCalculation,
            List<int[]> aNanPositions
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
            final String methodName = "setDescriptorsForMoleculesByMoleculeParallelizationNew";
            Descriptor.validateDescriptors(methodName, aDescriptors);
            Descriptor.validateAtomContainerArray(methodName, anAtomContainerArray);
            Descriptor.validateMatrix(methodName, aDescriptors, aMatrix, anAtomContainerArray, aStartIndex);
            Descriptor.validateNanPositions(methodName, aNanPositions);
        //</editor-fold>

        try {
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += aDescriptors[i].getDescriptorComponentNumber();
            }

            AtomicBoolean tmpHasNaN = new AtomicBoolean(false);

            if (anIsParallelCalculation) {
                try {
                    // Advise by Oracle: Parallel streams should use the common Fork-join pool
                    IntStream.range(0, anAtomContainerArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                boolean tmpSuccess = Descriptor.setDescriptorsForSingleMoleculeNew(
                                        aDescriptors,
                                        anAtomContainerArray[i],
                                        aMatrix[i],
                                        tmpStartIndices,
                                        i,
                                        aNanPositions
                                    );
                                if (!tmpSuccess) {
                                    tmpHasNaN.set(true);
                                }
                            } catch (Exception anException) {
                                tmpHasNaN.set(true);
                                Descriptor.LOGGER.log(
                                        Level.WARNING,
                                        "Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew: One descriptor calculation caused an exception, molecule index: "
                                                + i
                                                + ".",
                                        anException
                                );
                            }
                        }
                    );
                    return !tmpHasNaN.get();
                } catch (Exception anException) {
                    Descriptor.LOGGER.log(
                            Level.WARNING,
                            "Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew: Global exception occurred in descriptor calculation: ",
                            anException
                    );
                    return false;
                }
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < anAtomContainerArray.length; i++) {
                    if (!Descriptor.setDescriptorsForSingleMoleculeNew(aDescriptors, anAtomContainerArray[i], aMatrix[i], tmpStartIndices, i, aNanPositions)) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForMoleculesByMoleculeParallelizationNew: An exception occurred.", anException);
        }
    }
    //</editor-fold>

    //<editor-fold desc="Private static methods">
    /**
     * Sets calculated descriptor components in aVector (that corresponds to anAtomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in aVector to be filled with calculated components of descriptors
     * @param aMoleculeIndex Index of the current molecule being processed
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setDescriptorsForSingleMolecule (
            Descriptor[] aDescriptors,
            IAtomContainer anAtomContainer,
            float[] aVector,
            int[] aStartIndices,
            int aMoleculeIndex,
            List<int[]> aNanPositions
    ) throws Exception {
        try {
            IAtomContainer tmpMoleculeWithExplicitHydrogens = null;
            for (Descriptor aDescriptor : aDescriptors) {
                if (aDescriptor.needsExplicitHydrogens){
                    tmpMoleculeWithExplicitHydrogens = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    break;
                }
            }
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < aDescriptors.length; i++) {
                IAtomContainer moleculeToUse = aDescriptors[i].needsExplicitHydrogens ? tmpMoleculeWithExplicitHydrogens : anAtomContainer;
                if (!Descriptor.setDescriptor(aDescriptors[i], moleculeToUse, aVector, aStartIndices[i], aMoleculeIndex, aNanPositions)) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForSingleMolecule: An exception occurred: This should never happen.", anException);
        }
    }

    /**
     * Sets calculated descriptor components in aVector (that corresponds to anAtomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param aMoleculeSmilesString Molecule SMILES String (IS NOT CHANGED)
     * @param aVector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in aVector to be filled with calculated components of descriptors
     * @param aMoleculeIndex Index of the current molecule being processed
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @param anElectronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     * Aromaticity.Model.Daylight is used as default.
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setDescriptorsForSingleMoleculeSmilesString (
            Descriptor[] aDescriptors,
            String aMoleculeSmilesString,
            float[] aVector,
            int[] aStartIndices,
            int aMoleculeIndex,
            ElectronDonation anElectronDonationModel,
            List<int[]> aNanPositions
    ) throws Exception {
        try {
            IAtomContainer tmpMolecule = Descriptor.SMILES_PARSER.parseSmiles(aMoleculeSmilesString);
            if (anElectronDonationModel == null) {
                Descriptor.setAromaticity(tmpMolecule, Aromaticity.Model.Daylight);
            } else {
                Descriptor.setAromaticity(tmpMolecule, anElectronDonationModel);
            }

            IAtomContainer tmpMoleculeWithExplicitHydrogens = null;
            for (Descriptor aDescriptor : aDescriptors) {
                if (aDescriptor.needsExplicitHydrogens){
                    tmpMoleculeWithExplicitHydrogens = createMoleculeWithExplicitHydrogens(tmpMolecule);
                    break;
                }
            }
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < aDescriptors.length; i++) {
                IAtomContainer moleculeToUse = aDescriptors[i].needsExplicitHydrogens ? tmpMoleculeWithExplicitHydrogens : tmpMolecule;
                if (!Descriptor.setDescriptor(aDescriptors[i], moleculeToUse, aVector, aStartIndices[i], aMoleculeIndex, aNanPositions)) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForSingleMoleculeSmilesString: An exception occurred: This should never happen.", anException);
        }
    }

    /**
     * Sets calculated descriptor components in aVector (that corresponds to anAtomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Uses a new descriptor instance for EVERY descriptor calculation.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in aVector to be filled with calculated components of descriptor
     * @param aMoleculeIndex Index of the current molecule being processed
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setDescriptorsForSingleMoleculeNew(
        Descriptor[] aDescriptors,
        IAtomContainer anAtomContainer,
        float[] aVector,
        int[] aStartIndices,
        int aMoleculeIndex,
        List<int[]> aNanPositions
    ) throws Exception {
        try {
            IAtomContainer tmpMoleculeWithExplicitHydrogens = null;
            for (Descriptor aDescriptor : aDescriptors) {
                if (aDescriptor.needsExplicitHydrogens){
                    tmpMoleculeWithExplicitHydrogens = createMoleculeWithExplicitHydrogens(anAtomContainer);
                    break;
                }
            }
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < aDescriptors.length; i++) {
                IAtomContainer moleculeToUse = aDescriptors[i].needsExplicitHydrogens ? tmpMoleculeWithExplicitHydrogens : anAtomContainer;
                if (!Descriptor.setDescriptorNew(aDescriptors[i], moleculeToUse, aVector, aStartIndices[i], aMoleculeIndex, aNanPositions)) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            throw new Exception("Descriptor.setDescriptorsForSingleMoleculeNew: An exception occurred: This should never happen.", anException);
        }
    }

    /**
     * Sets component values of aDescriptor for anAtomContainer in aVector beginning with aStartIndex.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule (row of data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @param aMoleculeIndex Index of the current molecule being processed (row index in data matrix)
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setDescriptor(
            Descriptor aDescriptor,
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex,
            int aMoleculeIndex,
            List<int[]> aNanPositions
    ) {
        try {
            aDescriptor.calculate(anAtomContainer, aVector, aStartIndex);

            // Check for NaN values in the calculated result and track them
            int numComponents = aDescriptor.getDescriptorComponentNumber();
            return !Descriptor.checkAndTrackNaNValues(aVector, aStartIndex, numComponents, aMoleculeIndex, aNanPositions);
        } catch (Exception anException) {
            int numComponents = aDescriptor.getDescriptorComponentNumber();
            for (int i = 0; i < numComponents; i++) {
                aVector[aStartIndex + i] = Float.NaN;
                // Track NaN position if aNanPositions is provided
                if (aNanPositions != null) {
                    aNanPositions.add(new int[]{aMoleculeIndex, aStartIndex + i});
                }
            }
            Descriptor.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setDescriptor: An exception occurred while calculating descriptor "
                            + aDescriptor
                            + " for molecule index "
                            + aMoleculeIndex
                            + ".",
                    anException
            );
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
     * @param aVector Vector of molecule (row in data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @param aMoleculeIndex Index of the current molecule being processed
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setDescriptorNew(
            Descriptor aDescriptor,
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex,
            int aMoleculeIndex,
            List<int[]> aNanPositions

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
                case ATOM_COUNT_C:
                    AtomCountDescriptor atomCountCDesc = new AtomCountDescriptor();
                    //set parameter to count carbon atoms
                    atomCountCDesc.setParameters(new Object[]{"C"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountCDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_H:
                    AtomCountDescriptor atomCountHDesc = new AtomCountDescriptor();
                    //set parameter to count hydrogen atoms
                    atomCountHDesc.setParameters(new Object[]{"H"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountHDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_N:
                    AtomCountDescriptor atomCountNDesc = new AtomCountDescriptor();
                    //set parameter to count nitrogen atoms
                    atomCountNDesc.setParameters(new Object[]{"N"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountNDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_O:
                    AtomCountDescriptor atomCountODesc = new AtomCountDescriptor();
                    //set parameter to count oxygen atoms
                    atomCountODesc.setParameters(new Object[]{"O"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountODesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_S:
                    AtomCountDescriptor atomCountSDesc = new AtomCountDescriptor();
                    //set parameter to count sulfur atoms
                    atomCountSDesc.setParameters(new Object[]{"S"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountSDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_P:
                    AtomCountDescriptor atomCountPDesc = new AtomCountDescriptor();
                    //set parameter to count phosphorus atoms
                    atomCountPDesc.setParameters(new Object[]{"P"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountPDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_F:
                    AtomCountDescriptor atomCountFDesc = new AtomCountDescriptor();
                    //set parameter to count fluorine atoms
                    atomCountFDesc.setParameters(new Object[]{"F"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountFDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_BR:
                    AtomCountDescriptor atomCountBrDesc = new AtomCountDescriptor();
                    //set parameter to count bromine atoms
                    atomCountBrDesc.setParameters(new Object[]{"Br"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountBrDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_CL:
                    AtomCountDescriptor atomCountClDesc = new AtomCountDescriptor();
                    //set parameter to count chlorine atoms
                    atomCountClDesc.setParameters(new Object[]{"Cl"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountClDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_I:
                    AtomCountDescriptor atomCountIDesc = new AtomCountDescriptor();
                    //set parameter to count iodine atoms
                    atomCountIDesc.setParameters(new Object[]{"I"});
                    aVector[aStartIndex] = (float) ((IntegerResult) atomCountIDesc.calculate(anAtomContainer).getValue()).intValue();
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
                case BOND_COUNT_SINGLE:
                    BondCountDescriptor bondCountSingleDesc = new BondCountDescriptor();
                    //set parameter to count single bonds
                    bondCountSingleDesc.setParameters(new Object[]{"s"});
                    aVector[aStartIndex] = (float) ((IntegerResult) bondCountSingleDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_DOUBLE:
                    BondCountDescriptor bondCountDoubleDesc = new BondCountDescriptor();
                    //set parameter to count double bonds
                    bondCountDoubleDesc.setParameters(new Object[]{"d"});
                    aVector[aStartIndex] = (float) ((IntegerResult) bondCountDoubleDesc.calculate(anAtomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_TRIPLE:
                    BondCountDescriptor bondCountTripleDesc = new BondCountDescriptor();
                    //set parameter to count triple bonds
                    bondCountTripleDesc.setParameters(new Object[]{"t"});
                    aVector[aStartIndex] = (float) ((IntegerResult) bondCountTripleDesc.calculate(anAtomContainer).getValue()).intValue();
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
                    DoubleArrayResult aLogPResult = (DoubleArrayResult) (new ALOGPDescriptor()).calculate(anAtomContainer).getValue();
                    aVector[aStartIndex] = (float) aLogPResult.get(0);// ALogP
                    aVector[aStartIndex + 1] = (float) aLogPResult.get(1);  // ALogP squared
                    aVector[aStartIndex + 2] = (float) aLogPResult.get(2);  // Molar Refractivity
                    break;
                case X_LOG_P:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new XLogPDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
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
                    LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
                    // Change parameters so that we can use our own setAromaticity method
                    largestPiSystemDescriptor.setParameters(new Object[] {false}); // checkAromaticity = false
                    aVector[aStartIndex] = (float) ((IntegerResult) largestPiSystemDescriptor.calculate(anAtomContainer).getValue()).intValue();
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
                    for (int i = 0; i < 20; i++) {
                        aVector[aStartIndex + i] = (float) aminoAcidCountResult.get(i);
                    }
                    break;
                case KIER_HALL_SMARTS:
                    IntegerArrayResult kierHallSmartsResult = (IntegerArrayResult) (new KierHallSmartsDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 79; i++) {
                        aVector[aStartIndex + i] = (float) kierHallSmartsResult.get(i);
                    }
                    break;
                case ECCENTRIC_CONNECTIVITY_INDEX:
                    aVector[aStartIndex] = (float) ((IntegerResult) (new EccentricConnectivityIndexDescriptor()).calculate(anAtomContainer).getValue()).intValue();
                    break;
                case MDE:
                    DoubleArrayResult aDoubleArrayResult = (DoubleArrayResult) (new MDEDescriptor()).calculate(anAtomContainer).getValue();
                    for (int i = 0; i < 19; i++) {
                        aVector[aStartIndex + i] = (float) aDoubleArrayResult.get(i);
                    }
                    break;
                case VABC:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new VABCDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case PUBCHEM_FINGERPRINTER:
                    try {
                        IBitFingerprint fingerprint = new PubchemFingerprinter(SilentChemObjectBuilder.getInstance()).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < PUBCHEM_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < PUBCHEM_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_0:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP0, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_0:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP0, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_2:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP2, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_2:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP2, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_4:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_4:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP4, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_6:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_6:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP6, 1024).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < 1024; i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;
                case MACCS_FINGERPRINTER:
                    try {
                        IBitFingerprint fingerprint = new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()).getBitFingerprint(anAtomContainer);
                        for (int i = 0; i < MACCS_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            aVector[aStartIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception anException) {
                        for (int i = 0; i < MACCS_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            aVector[aStartIndex + i] = Float.NaN;
                        }
                        LOGGER.log(Level.WARNING, anException.toString(), anException);
                    }
                    break;


                // Add new descriptor information here!
                default:
                    throw new UnsupportedOperationException(aDescriptor + ": This descriptor does not have a routine yet!");
            }
            // Check for NaN values in the calculated result and track them
            int numComponents = aDescriptor.getDescriptorComponentNumber();
            return !checkAndTrackNaNValues(aVector, aStartIndex, numComponents, aMoleculeIndex, aNanPositions);
        } catch (Exception anException) {
            int numComponents = aDescriptor.getDescriptorComponentNumber();
            for (int i = 0; i < numComponents; i++) {
                aVector[aStartIndex + i] = Float.NaN;
                // Track NaN position if aNanPositions is provided
                if (aNanPositions != null) {
                    aNanPositions.add(new int[]{aMoleculeIndex, aStartIndex + i});
                }
            }
            Descriptor.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setDescriptorNew: An exception occurred while calculating descriptor "
                            + aDescriptor
                            + " for molecule index "
                            + aMoleculeIndex
                            + ".",
                    anException
            );
            return false;
        }
    }

    //<editor-fold desc="Helper methods">

    /**
     * Sets fingerprint component values in aVector.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule (row in data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private void calculateFingerprintFromPool(IAtomContainer anAtomContainer, float[] aVector, int aStartIndex) throws CDKException, InterruptedException {
        IFingerprinter fingerprinter = null;
        try {
            BlockingQueue<IFingerprinter> pool = fingerprintPoolMap.get(this);
            fingerprinter = pool.take();

            IBitFingerprint bitFingerprint = fingerprinter.getBitFingerprint(anAtomContainer);
            for (int i = 0; i < this.getDescriptorComponentNumber(); i++) {
                aVector[aStartIndex + i] = bitFingerprint.get(i) ? 1.0f : 0.0f;
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new InterruptedException("Thread interrupted during fingerprint calculation.");
        } catch (Exception anException) {
            throw new CDKException("Failed to calculate: " + this.name(), anException);
        } finally {
            fingerprintPoolMap.get(this).offer(fingerprinter);
        }
    }

    /**
     * Helper method to check for NaN values in a calculated descriptor result and track their positions.
     * This method is thread-safe when used with a synchronized LinkedList.
     *
     * @param aVector The vector containing calculated descriptor values
     * @param aStartIndex The start index in the vector for this descriptor
     * @param numComponents The number of components for this descriptor
     * @param aMoleculeIndex The index of the current molecule
     * @param aNanPositions List to track NaN positions (can be null if NaN positions should not be tracked);
     *                      must be thread-safe for parallel access; will be filled with int[] {[moleculeIndex, componentIndex]}
     *                      pairs of NaN values
     * @return true if any NaN values were found, false otherwise
     */
    private static boolean checkAndTrackNaNValues(float[] aVector, int aStartIndex, int numComponents,
                                          int aMoleculeIndex, List<int[]> aNanPositions) {
        boolean foundNaN = false;
        for (int i = 0; i < numComponents; i++) {
            if (Float.isNaN(aVector[aStartIndex + i])) {
                foundNaN = true;
                // Track NaN position if aNanPositions is provided
                // Note: aNanPositions must be thread-safe for parallel access
                if (aNanPositions != null) {
                    // Synchronize the add operation to ensure thread safety
                    aNanPositions.add(new int[]{aMoleculeIndex, aStartIndex + i});
                }
            }
        }
        return foundNaN;
    }

    //<editor-fold desc="Validation Methods">
    /**
     * Validates descriptor array input.
     *
     * @param methodName the calling method name for error messages
     * @param aDescriptorArray the descriptor array to validate
     * @throws NullPointerException if descriptor array is null
     * @throws IllegalArgumentException if descriptor array is empty
     */
    private static void validateDescriptors(String methodName, Descriptor[] aDescriptorArray){
        if (aDescriptorArray == null) {
            throw new NullPointerException(methodName + ": Given descriptor array is null");
        }
        if (aDescriptorArray.length == 0) {
            throw new IllegalArgumentException(methodName + ": Given descriptor array is empty");
        }
        for (Descriptor tmpDescriptor : aDescriptorArray) {
            if (tmpDescriptor == null) {
                throw new IllegalArgumentException(methodName + ": A single descriptor in aDescriptors is null.");
            }
        }
    }

    /**
     * Validates an atom container array input and its molecules.
     *
     * @param methodName the calling method name for error messages
     * @param anAtomContainerArray the atom container to validate
     * @throws NullPointerException if molecule collection is null
     * @throws IllegalArgumentException if molecule collection is empty
     */
    private static void validateAtomContainerArray(String methodName, IAtomContainer[] anAtomContainerArray) {
        if (anAtomContainerArray == null || anAtomContainerArray.length == 0) {
            throw new IllegalArgumentException(methodName + ": anAtomContainerArray is null or has length 0.");
        }
        for (IAtomContainer tmpMolecule : anAtomContainerArray) {
            if (tmpMolecule == null || tmpMolecule.isEmpty()) {
                throw new IllegalArgumentException(methodName + ": A single molecule in anAtomContainerArray is null or empty.");
            }
        }
    }

    /**
     * Validates a string array input and its SMILES strings.
     *
     * @param methodName the calling method name for error messages
     * @param aMoleculeSmilesStringArray the molecule smiles string array to validate
     * @throws NullPointerException if molecule collection is null
     * @throws IllegalArgumentException if molecule collection is empty
     */
    private static void validateSmilesStringArray(String methodName, String[] aMoleculeSmilesStringArray) {
        if (aMoleculeSmilesStringArray == null || aMoleculeSmilesStringArray.length == 0) {
            throw new IllegalArgumentException(methodName + ": anAtomContainerArray is null or has length 0.");
        }
        for (String tmpMolecule : aMoleculeSmilesStringArray) {
            if (tmpMolecule == null || tmpMolecule.isEmpty()) {
                throw new IllegalArgumentException(methodName + ": A single molecule in anAtomContainerArray is null or empty.");
            }
        }
    }

    /**
     * Validates NaN positions list.
     *
     * @param methodName the calling method name for error messages
     * @param aNanPositions the NaN positions list (can be null)
     */
    private static void validateNanPositions(String methodName, List<int[]> aNanPositions) {
        if (aNanPositions == null) {
            throw new IllegalArgumentException(methodName + ": aNaNPositions is null.");
        }
    }

    /**
     * Validates batch size.
     *
     * @param methodName the calling method name for error messages
     * @param aBatchSize a batch size (>= 0)
     */
    private static void validateBatchSize(String methodName, int aBatchSize) {
        if (aBatchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": aBatchSize must be greater than 0.");
        }
    }

    /**
     * Validates matrix dimensions and start index.
     *
     * @param aDescriptors the descriptors to validate
     * @param methodName the calling method name for error messages
     * @param aMatrix the matrix to validate
     * @param anAtomContainerArray the atom container array to validate
     * @param aStartIndex the start index in the matrix
     * @throws NullPointerException if matrix is null
     * @throws IllegalArgumentException if matrix dimensions are invalid
     * @throws ArrayIndexOutOfBoundsException if start index is out of bounds
     */
    private static void validateMatrix(String methodName, Descriptor[] aDescriptors, float[][] aMatrix, IAtomContainer[] anAtomContainerArray, int aStartIndex) throws Exception {
        if (aMatrix == null || aMatrix.length == 0) {
            throw new IllegalArgumentException(methodName + ": aMatrix is null or has length 0.");
        }
        if (aMatrix.length != anAtomContainerArray.length) {
            throw new IllegalArgumentException(methodName + ": aMatrix and anAtomContainerArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                throw new IllegalArgumentException(methodName + ": A vector in aMatrix is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                throw new IllegalArgumentException(methodName + ": aStartIndex is greater than or equal to vector length.");
            }
            int tmpNumberOfComponents = Descriptor.getNumberOfComponents(aDescriptors);
            if (aStartIndex + tmpNumberOfComponents > tmpVector.length) {
                throw new IllegalArgumentException(methodName + ": Not enough space in vector for descriptors.");
            }
        }
    }

    /**
     * Validates matrix dimensions and start index.
     *
     * @param aDescriptors the descriptors to validate
     * @param methodName the calling method name for error messages
     * @param aMatrix the matrix to validate
     * @param aMoleculeSmilesStringArray the smiles string array to validate
     * @param aStartIndex the start index in the matrix
     * @throws NullPointerException if matrix is null
     * @throws IllegalArgumentException if matrix dimensions are invalid
     * @throws ArrayIndexOutOfBoundsException if start index is out of bounds
     */
    private static void validateMatrix(String methodName, Descriptor[] aDescriptors, float[][] aMatrix, String[] aMoleculeSmilesStringArray, int aStartIndex) throws Exception {
        if (aMatrix == null || aMatrix.length == 0) {
            throw new IllegalArgumentException(methodName + ": aMatrix is null or has length 0.");
        }
        if (aMatrix.length != aMoleculeSmilesStringArray.length) {
            throw new IllegalArgumentException(methodName + ": aMatrix and anAtomContainerArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                throw new IllegalArgumentException(methodName + ": A vector in aMatrix is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                throw new IllegalArgumentException(methodName + ": aStartIndex is greater than or equal to vector length.");
            }
            int tmpNumberOfComponents = Descriptor.getNumberOfComponents(aDescriptors);
            if (aStartIndex + tmpNumberOfComponents > tmpVector.length) {
                throw new IllegalArgumentException(methodName + ": Not enough space in vector for descriptors.");
            }
        }
    }
    //</editor-fold>

    //<editor-fold desc="Package privat static molecule processing methods">
    /**
     * Creates a deep copy of the input molecule.
     * Note: This method is used to create a new molecule object
     * without affecting implicit hydrogen atoms.
     * Note: If necessary, atom types must be perceived and configured manually after creation.
     *
     * @param aMolecule Source molecule to be copied (NOT MODIFIED)
     * @return New instance of the molecule
     * @throws NullPointerException If the input molecule is null
     * @throws IllegalArgumentException If the input molecule is empty
     * @throws CloneNotSupportedException If the molecule cannot be properly copied
     */
    protected static IAtomContainer copyMolecule(IAtomContainer aMolecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Checks">
        if (aMolecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Input molecule must not be empty");
        }
        //</editor-fold>
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
     * @return New molecule (copy of the parameter) with all hydrogens made explicit
     * @throws NullPointerException If the input molecule is null
     * @throws IllegalArgumentException If the input molecule is empty
     * @throws CloneNotSupportedException If the molecule cannot be properly processed
     */
    protected static IAtomContainer createMoleculeWithExplicitHydrogens(
            IAtomContainer aMolecule
    ) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        //<editor-fold desc="Checks">
        if (aMolecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (aMolecule.isEmpty()) {
            throw new IllegalArgumentException("Input molecule must not be empty");
        }
        //</editor-fold>
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
    //</editor-fold>
    //</editor-fold>
    //</editor-fold>
}

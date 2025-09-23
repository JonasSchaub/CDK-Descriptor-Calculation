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

import static de.unijena.cheminf.clustering.desccalc.Descriptor.ACIDIC_GROUP_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.AMINO_ACID_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.AROMATIC_ATOMS_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.AROMATIC_BONDS_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.ATOM_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.ATOM_COUNT_ORGANIC_SUBSET;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.AUTOCORRELATION_CHARGE;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.AUTOCORRELATION_MASS;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.AUTOCORRELATION_POLARIZABILITY;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.A_LOG_P;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.A_POL;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.BASIC_GROUP_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.BCUT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.BOND_COUNT_ALL;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.BOND_COUNT_SPECIFIED;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.B_POL;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.CARBON_TYPES;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.CHI_CHAIN;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.CHI_CLUSTER;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.CHI_PATH;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.CHI_PATH_CLUSTER;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.ECCENTRIC_CONNECTIVITY_INDEX;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.FMF;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.FRACTIONAL_CSP3;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.FRACTIONAL_PSA;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.FRAGMENT_COMPLEXITY;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.HYBRIDIZATION_RATIO;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.H_BOND_ACCEPTOR_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.H_BOND_DONOR_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.JP_LOG_P;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.KAPPA_SHAPE_INDICES;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.KIER_HALL_SMARTS;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.LARGEST_CHAIN;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.LARGEST_PI_SYSTEM;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.LONGEST_ALIPHATIC_CHAIN;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.MANNHOLD_LOGP;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.MDE;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.MOLECULAR_WEIGHT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.PETITJEAN_NUMBER;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.ROTATABLE_BONDS_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.RULE_OF_FIVE;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.SMALL_RING;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.SPIRO_ATOM_COUNT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.TPSA;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.VABC;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.V_ADJ_MAT;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.WEIGHTED_PATH;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.WIENER_NUMBER;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.X_LOG_P;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.ZAGREB_INDEX;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.createMoleculeWithExplicitHydrogens;
import static de.unijena.cheminf.clustering.desccalc.Descriptor.getDescriptorToCdkObjectMap;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * This class provides descriptor calculation methods WITHOUT thread-safety guarantees!
 * CRITICAL WARNING: THREAD-SAFETY VIOLATIONS WILL CAUSE DATA CORRUPTION!
 * This class is designed ONLY for high-performance single-threaded calculations where thread safety
 * is not required or is handled externally. The methods in this class directly access shared CDK
 * descriptor instances without ANY synchronization, making them faster but EXTREMELY DANGEROUS
 * in multi-threaded environments.
 *
 * If you need thread-safe descriptor calculations, use the regular Descriptor class instead!
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
@Deprecated
public class DescriptorCalculator {

    /**
     * Private constructor to prevent instantiation as this is a utility class
     * with only static methods.
     */
    private DescriptorCalculator() {

    }

    //<editor-fold desc="Private static final LOGGER">
    /**
     * Logger of this class.
     */
    private static final Logger LOGGER = Logger.getLogger(DescriptorCalculator.class.getName());
    //</editor-fold>

    //<editor-fold desc="Private static final descriptor instances for descriptors with parameters">
    /**
     * Elements considered in the organic subset (C, H, N, O, S, P, F, Br, Cl, I) for the AtomCountDescriptor
     */
    private static final String[] ORGANIC_SUBSET_ELEMENTS = {
            "C", "H", "N", "O", "S", "P", "F", "Br", "Cl", "I"
    };

    /**
     * Bond types for BondCountDescriptor (single, double, triple bonds)
     */
    private static final String[] BOND_TYPES = {
            "s", "d", "t"  // single, double, triple
    };

    /**
     * Pre-created AtomCountDescriptor instances for each element in the organic subset to avoid
     * repeated instantiation during calculations.
     */
    private static final AtomCountDescriptor[] ORGANIC_SUBSET_DESCRIPTORS;
    /**
     * Pre-created BondCountDescriptor instances for each bond type to avoid
     * repeated instantiation during calculations.
     */
    private static final BondCountDescriptor[] BOND_COUNT_DESCRIPTORS;

    // Initialize BOND_COUNT_DESCRIPTORS in the static block (add to existing block)
    static {
        try {
            ORGANIC_SUBSET_DESCRIPTORS = createOrganicSubsetDescriptors();
            BOND_COUNT_DESCRIPTORS = createBondCountDescriptors();
        } catch (CDKException e) {
            throw new RuntimeException(e);
        }
    }

    // Method to create AtomCountDescriptor instances for the organic subset elements
    private static AtomCountDescriptor[] createOrganicSubsetDescriptors() throws CDKException {
        AtomCountDescriptor[] descriptors = new AtomCountDescriptor[ORGANIC_SUBSET_ELEMENTS.length];
        for (int i = 0; i < ORGANIC_SUBSET_ELEMENTS.length; i++) {
            AtomCountDescriptor desc = new AtomCountDescriptor();
            try {
                desc.setParameters(new Object[]{ORGANIC_SUBSET_ELEMENTS[i]});
            } catch (CDKException e) {
                throw new RuntimeException(e);
            }
            descriptors[i] = desc;
        }
        return descriptors;
    }

    // Method to create BondCountDescriptor instances for the bond types
    private static BondCountDescriptor[] createBondCountDescriptors() throws CDKException {
        BondCountDescriptor[] descriptors = new BondCountDescriptor[BOND_TYPES.length];
        for (int i = 0; i < BOND_TYPES.length; i++) {
            BondCountDescriptor desc = new BondCountDescriptor();
            try {
                desc.setParameters(new Object[]{BOND_TYPES[i]});
            } catch (CDKException e) {
                throw new RuntimeException(e);
            }
            descriptors[i] = desc;
        }
        return descriptors;
    }

    //</editor-fold>

    /**
     * Sets calculated descriptor components in vectors (rows) of a aMatrix (that corresponds to anAtomContainerArray)
     * beginning with aStartIndex by (optional) parallelization of molecules. If parallel computation is used, the atom
     * containers (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     *
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anAtomContainerArray Array of molecules. Note: anAtomContainerArray[i] corresponds to aMatrix[i] data
     *                              vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: Data vector aMatrix[i] corresponds to molecule
     *               anAtomContainerArray[i]. (MAY BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param anIsParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param aNanPositions List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (anIsParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
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
        if (aDescriptors == null || aDescriptors.length == 0) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aDescriptors is null or has length 0."
            );
            throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aDescriptor is null or has length 0.");
        }
        for (Descriptor tmpDescriptor : aDescriptors) {
            if (tmpDescriptor == null) {
                DescriptorCalculator.LOGGER.log(
                        Level.SEVERE,
                        "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: A descriptor in aDescriptors is null."
                );
                throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: A descriptor in aDescriptors is null.");
            }
        }
        if (anAtomContainerArray == null || anAtomContainerArray.length == 0) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: anAtomContainerArray is null or has length 0."
            );
            throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: anAtomContainerArray is null or has length 0.");
        }
        for (IAtomContainer tmpMolecule : anAtomContainerArray) {
            if (tmpMolecule == null || tmpMolecule.isEmpty()) {
                DescriptorCalculator.LOGGER.log(
                        Level.SEVERE,
                        "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: A molecule in anAtomContainerArray is null or empty."
                );
                throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: A molecule in anAtomContainerArray is null or empty.");
            }
        }
        if (aMatrix == null || aMatrix.length == 0) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix is null or has length 0."
            );
            throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix is null or has length 0.");
        }
        if (aMatrix.length != anAtomContainerArray.length) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix and anAtomContainerArray must have the same length."
            );
            throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aMatrix and anAtomContainerArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                DescriptorCalculator.LOGGER.log(
                        Level.SEVERE,
                        "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: A vector is null or has length 0."
                );
                throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: A vector is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                DescriptorCalculator.LOGGER.log(
                        Level.SEVERE,
                        "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
            }
            try {
                if (aStartIndex + Descriptor.getNumberOfComponents(aDescriptors) > tmpVector.length) {
                    DescriptorCalculator.LOGGER.log(
                            Level.SEVERE,
                            "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal."
                    );
                    throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.");
                }
            } catch (Exception anException) {
                DescriptorCalculator.LOGGER.log(
                        Level.SEVERE,
                        "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.", anException
                );
                throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aStartIndex is illegal.", anException);
            }
        }
        if (aNanPositions == null) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aNaNPositions is null."
            );
            throw new IllegalArgumentException("DescriptorCalculator.setCalculatedDescriptorComponentsByMoleculeParallelization: aNaNPositions is null.");
        }
        //</editor-fold>

        try {
            int[] tmpStartIndices = new int[aDescriptors.length];
            for (int i = 0; i < aDescriptors.length; i++) {
                tmpStartIndices[i] = aStartIndex;
                aStartIndex += aDescriptors[i].getDescriptorComponentNumber();
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
                                            DescriptorCalculator.setDescriptorsForSingleMolecule(
                                                    aDescriptors,
                                                    anAtomContainerArray[i],
                                                    aMatrix[i],
                                                    tmpStartIndices,
                                                    i,
                                                    aNanPositions
                                            );
                                } catch (Exception anException) {
                                    tmpIsDescriptorCalculations[i] = false;
                                    DescriptorCalculator.LOGGER.log(
                                            Level.WARNING,
                                            "DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelizationUnsafe: One descriptor calculation caused an exception, molecule index: "
                                                    + i
                                                    + ".",
                                            anException
                                    );
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
                    DescriptorCalculator.LOGGER.log(
                            Level.WARNING,
                            "DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelizationUnsafe: Global exception occurred in descriptor calculation: ",
                            anException
                    );
                    return false;
                }
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < anAtomContainerArray.length; i++) {
                    if (!DescriptorCalculator.setDescriptorsForSingleMolecule(aDescriptors, anAtomContainerArray[i], aMatrix[i], tmpStartIndices, i, aNanPositions)) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelizationUnsafe: An exception occurred: This should never happen.", anException
            );
            throw anException;
        }
    }
    /**
     * Sets calculated descriptor components in aVector (that corresponds to anAtomContainer, a row in the data matrix)
     * at aStartIndices.
     * WARNING: This method is NOT thread-safe and calls unsafe descriptor calculation methods.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
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
            boolean tmpIsSuccessful = true;
            for (int i = 0; i < aDescriptors.length; i++) {
                if (!DescriptorCalculator.setDescriptor(aDescriptors[i], anAtomContainer, aVector, aStartIndices[i], aMoleculeIndex, aNanPositions)) {
                    tmpIsSuccessful = false;
                }
            }
            return tmpIsSuccessful;
        } catch (Exception anException) {
            DescriptorCalculator.LOGGER.log(
                    Level.SEVERE,
                    "DescriptorCalculator.setCalculatedDescriptorComponents: An exception occurred: This should never happen.", anException
            );
            throw new Exception("DescriptorCalculator.setCalculatedDescriptorComponents: An exception occurred: This should never happen.", anException);
        }
    }

    /**
     * Sets component values of aDescriptor for anAtomContainer in aVector beginning with aStartIndex.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
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
                case ATOM_COUNT_ORGANIC_SUBSET:
                    setAtomCountOrganicSubset(anAtomContainer, aVector, aStartIndex);
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
                case ECCENTRIC_CONNECTIVITY_INDEX:
                    setEccentricConnectivityIndex(anAtomContainer, aVector, aStartIndex);
                    break;
                case MDE:
                    setMDE(anAtomContainer, aVector, aStartIndex);
                    break;
                case VABC:
                    setVABC(anAtomContainer, aVector, aStartIndex);
                    break;

                // Add new descriptor information here!
                default:
                    throw new UnsupportedOperationException(aDescriptor + ": This descriptor does not have a routine yet!");
            }
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
            DescriptorCalculator.LOGGER.log(
                    Level.WARNING,
                    "DescriptorCalculator.setDescriptor: An exception occurred while calculating descriptor "
                            + aDescriptor
                            + " for molecule index "
                            + aMoleculeIndex
                            + ".",
                    anException
            );
            return false;
        }
    }

    //<editor-fold desc="Private static UNSAFE descriptor calculation methods">
    /**
     * Sets molecular weight.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setMolecularWeight(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(MOLECULAR_WEIGHT).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets Wiener number(s).
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setWienerNumber(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult tmpResult = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(WIENER_NUMBER).calculate(anAtomContainer).getValue();
        aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
        aVector[aStartIndex + 1] = (float) tmpResult.get(1); //Wiener polarity number
    }

    /**
     * Sets atom count.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAtomCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(ATOM_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets atom count for organic subset (C, H, N, O, S, P, F, Br, Cl, I).
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAtomCountOrganicSubset(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            for (int i = 0; i < ORGANIC_SUBSET_DESCRIPTORS.length; i++) {
                aVector[aStartIndex + i] = (float) ((IntegerResult) ORGANIC_SUBSET_DESCRIPTORS[i].calculate(anAtomContainer).getValue()).intValue();
            }
        } catch (Exception e) {
            for (int i = 0; i < 10; i++) {
                aVector[aStartIndex + i] = Float.NaN;
            }
            DescriptorCalculator.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setAtomCountOrganicSubset: An exception occurred while calculating atom counts for molecule index "
                            + aStartIndex
                            + ".",
                    e
            );
        }
    }

    /**
     * Sets hydrogen bond acceptor count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setHBondAcceptorCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(H_BOND_ACCEPTOR_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets hydrogen bond donor count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setHBondDonorCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(H_BOND_DONOR_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the topological polar surface area (TPSA).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setTPSA(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(TPSA).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets largest chain size.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setLargestChain(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(LARGEST_CHAIN).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets longest aliphatic chain size.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setLongestAliphaticChain(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(LONGEST_ALIPHATIC_CHAIN).calculate(anAtomContainer).getValue()).intValue();
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
    private static void setMannholdLogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(MANNHOLD_LOGP).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets BCUT descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setBCUT(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(BCUT).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 6; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets bond count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setBondCountAll(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(BOND_COUNT_ALL).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets specific bond counts (single, double, triple).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setBondCountSpecified(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            // Single bonds (s)
            aVector[aStartIndex] = (float) ((IntegerResult) BOND_COUNT_DESCRIPTORS[0].calculate(anAtomContainer).getValue()).intValue();
            // Double bonds (d)
            aVector[aStartIndex + 1] = (float) ((IntegerResult) BOND_COUNT_DESCRIPTORS[1].calculate(anAtomContainer).getValue()).intValue();
            // Triple bonds (t)
            aVector[aStartIndex + 2] = (float) ((IntegerResult) BOND_COUNT_DESCRIPTORS[2].calculate(anAtomContainer).getValue()).intValue();
        } catch (Exception e) {
            for (int i = 0; i < 3; i++) {
                aVector[aStartIndex + i] = Float.NaN;
            }
            DescriptorCalculator.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setBondCountSpecified: An exception occurred while calculating bond counts for molecule index "
                            + aStartIndex
                            + ".",
                    e
            );
        }
    }

    /**
     * Sets bond polarizability value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setBPol(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(B_POL).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets Lipinski's Rule of Five violations count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setRuleOfFive(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(RULE_OF_FIVE).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the number of aromatic atoms.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAromaticAtomsCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(AROMATIC_ATOMS_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the number of aromatic bonds.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAromaticBondsCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(AROMATIC_BONDS_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets rotatable bonds count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setRotatableBondsCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(ROTATABLE_BONDS_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the FMF (Framework Match Fraction) value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setFMF(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(FMF).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets the FractionalCSP3 value (fraction of sp3 hybridized carbon atoms).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setFractionalCSP3(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(FRACTIONAL_CSP3).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets hybridization ratio (sp3 carbons to sp2 carbons).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setHybridizationRatio(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            // Create a copy of the molecule with explicit hydrogen atoms
            IAtomContainer moleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
            // Calculate XLogP using the molecule copy containing explicit hydrogen atoms
            aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(HYBRIDIZATION_RATIO).calculate(moleculeWithExplicitH).getValue()).doubleValue();
        } catch (Exception anException) {
            aVector[aStartIndex] = Float.NaN;
            DescriptorCalculator.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setHybridizationRatio: An exception occurred while calculating hybridization ratio for molecule index "
                            + aStartIndex + ".",
                    anException
            );
        }
    }

    /**
     * Sets Kappa shape indices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setKappaShapeIndices(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(KAPPA_SHAPE_INDICES).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 3; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets Petitjean number value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setPetitjeanNumber(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(PETITJEAN_NUMBER).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets spiro atom count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setSpiroAtomCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(SPIRO_ATOM_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets the vertex adjacency information (magnitude).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setVAdjMat(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(V_ADJ_MAT).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets weighted path descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setWeightedPath(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(WEIGHTED_PATH).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets Zagreb index value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setZagrebIndex(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(ZAGREB_INDEX).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets carbon types descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setCarbonTypesDescriptor(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult result = (IntegerArrayResult) getDescriptorToCdkObjectMap().get(CARBON_TYPES).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 9; i++) {
            aVector[aStartIndex + i] = result.get(i);
        }
    }

    /**
     * Sets ALogP values (Ghose-Crippen LogP, ALogP squared and molar refractivity).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setALogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            // Create a copy of the molecule with explicit hydrogen atoms
            IAtomContainer moleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
            // Calculate ALogP using the molecule copy containing explicit hydrogen atoms
            DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(A_LOG_P).calculate(moleculeWithExplicitH).getValue();
            aVector[aStartIndex] = (float) result.get(0);      // ALogP
            aVector[aStartIndex + 1] = (float) result.get(1);  // ALogP squared
            aVector[aStartIndex + 2] = (float) result.get(2);  // Molar Refractivity
        } catch (Exception anException) {
            aVector[aStartIndex] = Float.NaN;
            aVector[aStartIndex + 1] = Float.NaN;
            aVector[aStartIndex + 2] = Float.NaN;
            DescriptorCalculator.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setALogP: An exception occurred while calculating ALogP for molecule index "
                            + aStartIndex
                            + ".",
                    anException
            );
        }
    }

    /**
     * Sets XLogP value (prediction of logP based on the atom-type method).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     * Note: XLogP requires explicit hydrogens for correct calculation
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setXLogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        try {
            // Create a copy of the molecule with explicit hydrogen atoms
            IAtomContainer moleculeWithExplicitH = createMoleculeWithExplicitHydrogens(anAtomContainer);
            // Calculate XLogP using the molecule copy containing explicit hydrogen atoms
            aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(X_LOG_P).calculate(moleculeWithExplicitH).getValue()).doubleValue();
        } catch (Exception anException) {
            aVector[aStartIndex] = Float.NaN;
            DescriptorCalculator.LOGGER.log(
                    Level.WARNING,
                    "Descriptor.setXLogP: An exception occurred while calculating XLogP for molecule index "
                            + aStartIndex
                            + ".",
                    anException
            );
        }
    }
    /**
     * Sets JP LogP value (octanol-water partition coefficient based on JPlogP method).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setJPLogP(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(JP_LOG_P).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets APol value (sum of the atomic polarizabilities).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAPol(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(A_POL).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets AutocorrelationDescriptorCharge values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAutocorrelationCharge(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(AUTOCORRELATION_CHARGE).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets AutocorrelationDescriptorMass values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAutocorrelationMass(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(AUTOCORRELATION_MASS).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets AutocorrelationDescriptorPolarizability values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAutocorrelationPolarizability(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(AUTOCORRELATION_POLARIZABILITY).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 5; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets fragment complexity value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setFragmentComplexity(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(FRAGMENT_COMPLEXITY).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets ChiChain descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setChiChain(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(CHI_CHAIN).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 10; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets ChiCluster descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setChiCluster(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(CHI_CLUSTER).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 8; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets ChiPathCluster descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setChiPathCluster(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(CHI_PATH_CLUSTER).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 6; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets ChiPath descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setChiPath(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult result = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(CHI_PATH).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 16; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets FractionalPSA descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     *  @param anAtomContainer Molecule (IS NOT CHANGED)
     *  @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     *  @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setFractionalPSA(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(FRACTIONAL_PSA).calculate(anAtomContainer).getValue()).doubleValue();
    }

    /**
     * Sets LargestPiSystem descriptor value.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setLargestPiSystem(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(LARGEST_PI_SYSTEM).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets SmallRing descriptor values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setSmallRing(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult result = (IntegerArrayResult) getDescriptorToCdkObjectMap().get(SMALL_RING).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 11; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets basic group count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setBasicGroupCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(BASIC_GROUP_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets acidic group count.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAcidicGroupCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(ACIDIC_GROUP_COUNT).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets amino acid count values.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setAminoAcidCount(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult aminoAcidCountResult = (IntegerArrayResult) getDescriptorToCdkObjectMap().get(AMINO_ACID_COUNT).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 20; i++) {
            aVector[aStartIndex + i] = (float) aminoAcidCountResult.get(i);
        }
    }
    /**
     * Sets Kier-Hall SMARTS descriptor values (79 functional group counts).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setKierHallSmarts(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        IntegerArrayResult result = (IntegerArrayResult) getDescriptorToCdkObjectMap().get(KIER_HALL_SMARTS).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 79; i++) {
            aVector[aStartIndex + i] = (float) result.get(i);
        }
    }

    /**
     * Sets the EccentricConnectivityIndex descriptor value in aVector beginning with aStartIndex.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector to be filled with descriptor value (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector
     */
    private static void setEccentricConnectivityIndex(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((IntegerResult) getDescriptorToCdkObjectMap().get(ECCENTRIC_CONNECTIVITY_INDEX).calculate(anAtomContainer).getValue()).intValue();
    }

    /**
     * Sets MDE descriptor values (molecular distance edge between atoms of specific types).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setMDE(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        DoubleArrayResult tmpResult = (DoubleArrayResult) getDescriptorToCdkObjectMap().get(MDE).calculate(anAtomContainer).getValue();
        for (int i = 0; i < 19; i++) {
            aVector[aStartIndex + i] = (float) tmpResult.get(i);
        }
    }

    /**
     * Sets VABC value (volume descriptor based on atom contributions).
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     * WARNING: This method is NOT thread-safe and accesses shared CDK descriptor instances without synchronization.
     * Use only in single-threaded contexts or when thread safety is ensured externally.
     *
     * @param anAtomContainer Molecule (IS NOT CHANGED)
     * @param aVector Vector of molecule to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     */
    private static void setVABC(
            IAtomContainer anAtomContainer,
            float[] aVector,
            int aStartIndex
    ) {
        aVector[aStartIndex] = (float) ((DoubleResult) getDescriptorToCdkObjectMap().get(VABC).calculate(anAtomContainer).getValue()).doubleValue();
    }

    // Add new descriptor information here!
    //</editor-fold>

}

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

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;

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
    MOLECULER_WEIGHT,
    /**
     * Wiener number, returns Wiener path number and Wiener polarity number.
     * Path number: sum of the distances between any two atoms in the molecule.
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     */
    WIENER_NUMBER;

    // Add new descriptor information here!

    /**
     * EnumMap that maps a descriptor to its number of calculated components
     */
    private static final EnumMap<Descriptor, Integer> descriptorToComponentNumberMap = new EnumMap<>(Descriptor.class);
    /**
     * EnumMap that maps a descriptor to an instance of its CDK descriptor class
     */
    private static final EnumMap<Descriptor, IMolecularDescriptor> descriptorToCdkObjectMap = new EnumMap<>(Descriptor.class);
    static {
        // MOLECULER_WEIGHT has 1 component
        descriptorToComponentNumberMap.put(MOLECULER_WEIGHT, 1);
        descriptorToCdkObjectMap.put(MOLECULER_WEIGHT, new WeightDescriptor());

        // WIENER_NUMBER has 2 components, Wiener path number and Wiener polarity number
        descriptorToComponentNumberMap.put(WIENER_NUMBER, 2);
        descriptorToCdkObjectMap.put(WIENER_NUMBER, new WienerNumbersDescriptor());

        // Add new descriptor information here!

    }
    //</editor-fold>

    //<editor-fold desc="Private static final LOGGER">
    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(Descriptor.class.getName());
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
                case MOLECULER_WEIGHT:
                    setMolecularWeight(anAtomContainer, aVector, aStartIndex);
                    break;
                case WIENER_NUMBER:
                    setWienerNumber(anAtomContainer, aVector, aStartIndex);
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
                case MOLECULER_WEIGHT:
                    aVector[aStartIndex] = (float) ((DoubleResult) (new WeightDescriptor()).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case WIENER_NUMBER:
                    DoubleArrayResult tmpResult = (DoubleArrayResult) (new WienerNumbersDescriptor()).calculate(anAtomContainer).getValue();
                    aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
                    aVector[aStartIndex + 1] = (float) tmpResult.get(1); //Wiener polarity number
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
                case MOLECULER_WEIGHT:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MOLECULER_WEIGHT).calculate(anAtomContainer).getValue()).doubleValue();
                    break;
                case WIENER_NUMBER:
                    DoubleArrayResult tmpResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(WIENER_NUMBER).calculate(anAtomContainer).getValue();
                    aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
                    aVector[aStartIndex + 1] = (float) tmpResult.get(1); //Wiener polarity number
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
        aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MOLECULER_WEIGHT).calculate(anAtomContainer).getValue()).doubleValue();
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

    // Add new descriptor information here!
    //</editor-fold>

}

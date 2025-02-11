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

import java.util.EnumMap;
import java.util.logging.Level;
import java.util.logging.Logger;

// TODO Jonas: Implement IAtomContainer based methods instead of/in addition to SMILES based ones.
// TODO: Define useful descriptor profiles

/**
 * Descriptor related calculations based on the CDK
 * for the enrichment of data vectors for clustering
 * Note: For adding a new descriptor goto "Add new descriptor here!"
 *
 * @author Jonas Schaub, Achim Zielesny
 */
public enum Descriptor {

    //<editor-fold desc="Descriptor enumeration and initialization">
    /*
     * Molecular weight
     */
    MOLECULER_WEIGHT,
    /*
     * Wiener number
     */
    WIENER_NUMBER;
    // Add new descriptor here!

    /*
     * EnumMap that maps a descriptor to its number of calculated components
     */
    private static final EnumMap<Descriptor, Integer> descriptorToComponentNumberMap = new EnumMap<>(Descriptor.class);
    /*
     * EnumMap that maps a descriptor to its CDK object
     */
    private static final EnumMap<Descriptor, Object> descriptorToCdkObjectMap = new EnumMap<>(Descriptor.class);
    static {
        // MOLECULER_WEIGHT has 1 component
        descriptorToComponentNumberMap.put(MOLECULER_WEIGHT, 1);
        // descriptorToCdkObjectMap.put(MOLECULER_WEIGHT, <new CdkObjectForMolWeightCalculation>);

        // WIENER_INDEX has 1 component
        descriptorToComponentNumberMap.put(WIENER_NUMBER, 1);
        // descriptorToCdkObjectMap.put(WIENER_NUMBER, <new CdkObjectForWienerNumberCalculation>);

        // Add new descriptor here!
    }
    //</editor-fold>

    //<editor-fold desc="Private static final LOGGER">
    /*
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
     * Returns sum of number of calculated components of an array of defined descriptors
     *
     * @param aDescriptors Array of  descriptors
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
                throw new IllegalArgumentException("Descriptor.getNumberOfComponents: A descriptor in aDescriptors is null.");
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
     * Sets calculated descriptor components in vectors of a aMatrix (that corresponds to aSmilesArray) beginning with
     * aStartIndex.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param aSmilesArray Array of SMILES strings of molecules. Note: aSmilesArray[i] corresponds to aMatrix[i]. (IS
     *                     NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: aMatrix[i] corresponds to aSmilesArray[i]. (MAY
     *                BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors
     * @param aNumberOfConcurrentCalculationThreads Number of concurrent calculation threads. If zero, then the
     *                                              all calculations are performed one after another (sequentially).
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setCalculatedDescriptorComponents (
        Descriptor[] aDescriptors,
        String[] aSmilesArray,
        float[][] aMatrix,
        int aStartIndex,
        int aNumberOfConcurrentCalculationThreads
    ) throws IllegalArgumentException, Exception {
        //<editor-fold desc="Checks">
        if (aDescriptors == null || aDescriptors.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aDescriptors is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aDescriptor is null or has length 0.");
        }
        for (Descriptor tmpDescriptor : aDescriptors) {
            if (tmpDescriptor == null) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponents: A descriptor in aDescriptors is null."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: A descriptor in aDescriptors is null.");
            }
        }
        if (aSmilesArray == null || aSmilesArray.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aSmilesArray is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aSmilesArray is null or has length 0.");
        }
        for (String tmpSmiles : aSmilesArray) {
            if (tmpSmiles == null || tmpSmiles.isEmpty()) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponents: A SMILES in aSmilesArray is null or empty."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: A SMILES in aSmilesArray is null or empty.");
            }
        }
        if (aMatrix == null || aMatrix.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aMatrix is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aMatrix is null or has length 0.");
        }
        if (aMatrix.length != aSmilesArray.length) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aMatrix and aSmilesArray must have the same length."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aMatrix and aSmilesArray must have the same length.");
        }
        for (float[] tmpVector : aMatrix) {
            if (tmpVector == null || tmpVector.length == 0) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponents: A vector is null or has length 0."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: A vector is null or has length 0.");
            }
            if (aStartIndex >= tmpVector.length) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponents: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aStartIndex is illegal.");
            }
            try {
                if (aStartIndex + Descriptor.getNumberOfComponents(aDescriptors) > tmpVector.length) {
                    Descriptor.LOGGER.log(
                        Level.SEVERE,
                        "Descriptor.setCalculatedDescriptorComponents: aStartIndex is illegal."
                    );
                    throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aStartIndex is illegal.");
                }
            } catch (Exception anException) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponents: aStartIndex is illegal."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aStartIndex is illegal.");
            }
        }
        if (aNumberOfConcurrentCalculationThreads < 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aNumberOfConcurrentCalculationThreads is illegal."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aNumberOfConcurrentCalculationThreads is illegal.");
        }
        //</editor-fold>

        try {
            if (aNumberOfConcurrentCalculationThreads > 0) {
                // TODO: Implement parallelized calculation
                // TODO Jonas: Has single static CDK instance thread-safe descriptor calculations?
                return false;
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < aSmilesArray.length; i++) {
                    if (!Descriptor.setCalculatedDescriptorComponents(aDescriptors, aSmilesArray[i], aMatrix[i], aStartIndex)) {
                        tmpIsSuccessful = false;
                    }
                }
                return tmpIsSuccessful;
            }
        } catch (Exception anException) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: An exception occurred: This should never happen."
            );
            throw anException;
        }
    }
    //</editor-fold>

    //<editor-fold desc="Private static methods">
    /**
     * Sets calculated descriptor components in aVector (that corresponds to aSmiles) beginning with aStartIndex
     * Note: No checks are NOT performed. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated
     * @param aSmiles SMILES of molecule
     * @param aVector Component vector of molecule
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setCalculatedDescriptorComponents (
            Descriptor[] aDescriptors,
            String aSmiles,
            float[] aVector,
            int aStartIndex
    ) throws Exception {
        try {
            boolean tmpIsSuccessful = true;
            for (Descriptor tmpDescriptor : aDescriptors) {
                if (!Descriptor.setComponentValues(tmpDescriptor, aSmiles, aVector, aStartIndex)) {
                    tmpIsSuccessful = false;
                }
                aStartIndex += descriptorToComponentNumberMap.get(tmpDescriptor);
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
     * Sets component values of aDescriptor for aSmiles in aVector beginning with aStartIndex.
     * Note: No checks are NOT performed. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated
     * @param aSmiles SMILES of molecule
     * @param aVector Vector of molecule to be filled with calculated components of descriptors
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setComponentValues(
        Descriptor aDescriptor,
        String aSmiles,
        float[] aVector,
        int aStartIndex
    ) {
        try {
            switch (aDescriptor) {
                case MOLECULER_WEIGHT:
                    // TODO Jonas: Calculate components of descriptor for aSmiles with CDK and set aVector
                    // aVector[aStartIndex] = (float) ((<CdkMolWeightObject>) descriptorToCdkObjectMap.get(MOLECULER_WEIGHT)).getValue(aSmiles);
                    // Demo code:
                    aVector[aStartIndex] = 500.0f;
                    break;
                case WIENER_NUMBER:
                    // TODO Jonas: Calculate components of descriptor for aSmiles with CDK and set aVector
                    // aVector[aStartIndex] = (float) ((<CdkWienerNumberObject>) descriptorToCdkObjectMap.get(WIENER_NUMBER)).getValue(aSmiles);
                    // Demo code:
                    aVector[aStartIndex] = 70.0f;
                    break;
                // Add new descriptor here!
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

}

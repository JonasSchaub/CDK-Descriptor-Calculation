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

// TODO Jonas: Implement IAtomContainer based methods instead of/in addition to SMILES based ones. -> I would make the user take care of converting from SMILES to AC
// TODO: Define useful descriptor profiles

/**
 * Descriptor related calculations based on the CDK
 * for the enrichment of data vectors for clustering.
 * Note: For adding a new descriptor goto "Add new descriptor here!"
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
     * Wiener number, returns Wiener path nr and Wiener polarity nr.
     * Path number: sum of the distances between any two atoms in the molecule.
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     */
    WIENER_NUMBER;
    // Add new descriptor here!

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

        // WIENER_INDEX has 2 components, Wiener path nr. and Wiener polarity nr.
        descriptorToComponentNumberMap.put(WIENER_NUMBER, 2);
        descriptorToCdkObjectMap.put(WIENER_NUMBER, new WienerNumbersDescriptor());

        // Add new descriptor here!
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
     * Returns sum of number of calculated components of an array of defined descriptors
     *
     * @param aDescriptors Array of descriptors
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
     * Sets calculated descriptor components in vectors of a aMatrix (that corresponds to anIAtomContainerArray) beginning with
     * aStartIndex.
     *
     * @param aDescriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param anIAtomContainerArray Array of molecules. Note: anIAtomContainerArray[i] corresponds to aMatrix[i]. (IS
     *                     NOT CHANGED)
     * @param aMatrix Matrix of component vectors of molecules. Note: aMatrix[i] corresponds to anIAtomContainerArray[i]. (MAY
     *                BE CHANGED)
     * @param aStartIndex Start index in a vector to be filled with calculated components of descriptors
     * @param aNumberOfConcurrentCalculationThreads Number of concurrent calculation threads. If zero, then
     *                                              all calculations are performed one after another (sequentially). TODO: if 1, it would also be sequential, wouldn't it?
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a descriptor
     * calculation is NaN
     * @throws IllegalArgumentException Thrown if an argument is illegal
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    public static boolean setCalculatedDescriptorComponents (
        Descriptor[] aDescriptors,
        IAtomContainer[] anIAtomContainerArray,
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
        if (anIAtomContainerArray == null || anIAtomContainerArray.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: anIAtomContainerArray is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: anIAtomContainerArray is null or has length 0.");
        }
        for (IAtomContainer tmpMol : anIAtomContainerArray) {
            if (tmpMol == null || tmpMol.isEmpty()) {
                Descriptor.LOGGER.log(
                    Level.SEVERE,
                    "Descriptor.setCalculatedDescriptorComponents: A molecule in anIAtomContainerArray is null or empty."
                );
                throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: A molecule in anIAtomContainerArray is null or empty.");
            }
        }
        if (aMatrix == null || aMatrix.length == 0) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aMatrix is null or has length 0."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aMatrix is null or has length 0.");
        }
        if (aMatrix.length != anIAtomContainerArray.length) {
            Descriptor.LOGGER.log(
                Level.SEVERE,
                "Descriptor.setCalculatedDescriptorComponents: aMatrix and anIAtomContainerArray must have the same length."
            );
            throw new IllegalArgumentException("Descriptor.setCalculatedDescriptorComponents: aMatrix and anIAtomContainerArray must have the same length.");
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
                // TODO Jonas: Has single static CDK instance thread-safe descriptor calculations? -> don't know, we have to test whether it's thread-safe
                return false;
            } else {
                boolean tmpIsSuccessful = true;
                for (int i = 0; i < anIAtomContainerArray.length; i++) {
                    if (!Descriptor.setCalculatedDescriptorComponents(aDescriptors, anIAtomContainerArray[i], aMatrix[i], aStartIndex)) {
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
     * Sets calculated descriptor components in aVector (that corresponds to anIAtomContainer) beginning with aStartIndex
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptors Array of descriptors to be calculated
     * @param anIAtomContainer molecule
     * @param aVector Component vector of molecule
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     * @throws Exception Thrown if fatal error occurs (this should never happen)
     */
    private static boolean setCalculatedDescriptorComponents (
            Descriptor[] aDescriptors,
            IAtomContainer anIAtomContainer,
            float[] aVector,
            int aStartIndex
    ) throws Exception {
        try {
            boolean tmpIsSuccessful = true;
            for (Descriptor tmpDescriptor : aDescriptors) {
                if (!Descriptor.setComponentValues(tmpDescriptor, anIAtomContainer, aVector, aStartIndex)) {
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
     * Sets component values of aDescriptor for anIAtomContainer in aVector beginning with aStartIndex.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param aDescriptor Descriptor to be calculated
     * @param anIAtomContainer molecule
     * @param aVector Vector of molecule to be filled with calculated components of descriptors
     * @param aStartIndex Start index in aVector to be filled with calculated components of descriptors
     * @return True: Operation was successful, false: Operation failed, i.e. at least one component in a
     * descriptor calculation is NaN
     */
    private static boolean setComponentValues(
        Descriptor aDescriptor,
        IAtomContainer anIAtomContainer,
        float[] aVector,
        int aStartIndex
    ) {
        try {
            switch (aDescriptor) {
                case MOLECULER_WEIGHT:
                    aVector[aStartIndex] = (float) ((DoubleResult) descriptorToCdkObjectMap.get(MOLECULER_WEIGHT).calculate(anIAtomContainer).getValue()).doubleValue();
                    break;
                case WIENER_NUMBER:
                    DoubleArrayResult tmpResult = (DoubleArrayResult) descriptorToCdkObjectMap.get(WIENER_NUMBER).calculate(anIAtomContainer).getValue();
                    aVector[aStartIndex] = (float) tmpResult.get(0); //Wiener path number
                    aVector[aStartIndex+1] = (float) tmpResult.get(1); //Wiener polarity number
                    break;
                // Add new descriptor here!
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

}

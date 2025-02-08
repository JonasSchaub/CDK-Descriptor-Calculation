package de.unijena.cheminf.clustering.desccalc;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test class for descriptor
 */
class DescriptorTest {

    /**
     * Test method for descriptor MOLECULAR_WEIGHT
     */
    @Test
    public void test_MOLECULAR_WEIGHT() {
        // Ethanol
        String tmpSmiles = "CC(=O)O";
        String[] tmpSmilesArray = new String[] {tmpSmiles};
        float[][] tmpMatrix = new float[][]
            {
                {0f}
            };
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT};
        int tmpNumberOfConcurrentCalculationThreads = 0;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(Descriptor.setCalculatedDescriptorComponents(tmpDescriptors, tmpSmilesArray, tmpMatrix, tmpStartIndex, tmpNumberOfConcurrentCalculationThreads));
            Assertions.assertEquals(500.0f, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test method for descriptor WIENER_INDEX
     */
    @Test
    public void test_WIENER_INDEX() {
        // Ethanol
        String tmpSmiles = "CC(=O)O";
        String[] tmpSmilesArray = new String[] {tmpSmiles};
        float[][] tmpMatrix = new float[][]
                {
                        {0f}
                };
        int tmpStartIndex = 0;
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.WIENER_INDEX};
        int tmpNumberOfConcurrentCalculationThreads = 0;

        try {
            Assertions.assertEquals(1, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(Descriptor.setCalculatedDescriptorComponents(tmpDescriptors, tmpSmilesArray, tmpMatrix, tmpStartIndex, tmpNumberOfConcurrentCalculationThreads));
            Assertions.assertEquals(70.0f, tmpMatrix[0][0]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

    /**
     * Test multiple descriptors
     */
    @Test
    public void test_MultipleDescriptors() {
        // Ethanol
        String tmpSmiles = "CC(=O)O";
        String[] tmpSmilesArray = new String[] {tmpSmiles};
        int tmpStartIndex = 0;
        float[][] tmpMatrix = new float[][]
                {
                        {0f, 0f}
                };
        Descriptor[] tmpDescriptors = new Descriptor[] {Descriptor.MOLECULER_WEIGHT, Descriptor.WIENER_INDEX};
        int tmpNumberOfConcurrentCalculationThreads = 0;

        try {
            Assertions.assertEquals(2, Descriptor.getNumberOfComponents(tmpDescriptors));
            Assertions.assertTrue(Descriptor.setCalculatedDescriptorComponents(tmpDescriptors, tmpSmilesArray, tmpMatrix, tmpStartIndex, tmpNumberOfConcurrentCalculationThreads));
            Assertions.assertEquals(500.0f, tmpMatrix[0][0]);
            Assertions.assertEquals(70.0f, tmpMatrix[0][1]);
        } catch (Exception anException) {
            Assertions.fail();
        }
    }

}
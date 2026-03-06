package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.io.SeclibParser;
import java.util.ArrayList;
import java.util.List;

public class TestSetCreator {

    public static void main(String[] args) throws Exception {

        List<Protein> allProteins = SeclibParser.parseSeclib("cb513.db");

        //split: train: 400 , rest ist predict
        List<Protein> trainingSet = allProteins.subList(0, 400);
        List<Protein> potentialTestSet = allProteins.subList(400, allProteins.size());

        List<Protein> cleanTestSet = new ArrayList<>();


        //test protein und train protein vergleichen
        //exclude proteins with high sequence similarity to
        //the training set from your test set
        for (Protein testProt : potentialTestSet) {
            boolean isTooSimilar = false;

            for (Protein trainProt : trainingSet) {

                double identity = calculateAlignmentIdentity(testProt.getSequence(), trainProt.getSequence());

                if (identity > 0.25) { // 25% Threshold
                    isTooSimilar = true;
                    break;
                }
            }

            if (!isTooSimilar) {
                cleanTestSet.add(testProt);
            }
        }
        
        
    }

    //Alignemnt?
    private static double calculateAlignmentIdentity(String sequence, String sequence1) {
        return 0;
    }

}

package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.io.SeclibParser;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class TestSetCreator {

    public static void main(String[] args) throws Exception {

        String trainingPath = "C:/Users/gross/IdeaProjects/progprakt-17/tools/CB513DSSP.db";
        String inputPath = "C:/Users/gross/IdeaProjects/progprakt-17/tools/casp_all.seclib";
        String outputPath = "clean_test_casp.db";

        List<Protein> trainingSet = SeclibParser.parseSeclib(trainingPath);
        List<Protein> potentialTestSet = SeclibParser.parseSeclib(inputPath);
        List<Protein> cleanTestSet = new ArrayList<>();
        //grenze und mindestlänge, 0.25 was zu nuedrig
        double identityThreshold = 0.3;
        int minLength = 18;

        int processedCount = 0;

        for (Protein testProt : potentialTestSet) {
            processedCount++;

            // längen filter
            if (testProt.getSequence().length() < minLength) {
                continue;
            }

            boolean isTooSimilar = false;
            String seqTest = testProt.getSequence();

            //Homologie checken gg training set
            for (Protein trainProt : trainingSet) {
                double identity = calculateAlignmentIdentity(seqTest, trainProt.getSequence());

                if (identity > identityThreshold) {
                    isTooSimilar = true;
                    break; // Zu ähnlich -> Abbruch der inneren Schleife
                }
            }

            // wenn kein protein im training zu ähnlich war, ab ins test set
            if (!isTooSimilar) {
                cleanTestSet.add(testProt);
            }
        }


        //Ergebnisse speichern
        try (PrintWriter pw = new PrintWriter(outputPath)) {
            for (Protein p : cleanTestSet) {
                pw.println("> " + p.getId());
                pw.println("AS " + p.getSequence());
                pw.println("SS " + p.getStructure());
                pw.println();
            }
        }
    }

   //NW
    private static double calculateAlignmentIdentity(String seq1, String seq2) {
        int n = seq1.length();
        int m = seq2.length();
        int GAP = -1;

        // DP-Matrix
        int[][] F = new int[n + 1][m + 1];
        for (int i = 0; i <= n; i++) F[i][0] = i * GAP;
        for (int j = 0; j <= m; j++) F[0][j] = j * GAP;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int match = F[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? 1 : -1);
                int del = F[i - 1][j] + GAP;
                int ins = F[i][j - 1] + GAP;
                F[i][j] = Math.max(match, Math.max(del, ins));
            }
        }

        // traceback berechnung matrix
        int matches = 0;
        int alignedLength = 0;
        int i = n, j = m;

        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && F[i][j] == F[i - 1][j - 1] + (seq1.charAt(i - 1) == seq2.charAt(j - 1) ? 1 : -1)) {
                if (seq1.charAt(i - 1) == seq2.charAt(j - 1)) {
                    matches++;
                }
                i--; j--;
            } else if (i > 0 && F[i][j] == F[i - 1][j] + GAP) {
                i--;
            } else {
                j--;
            }
            alignedLength++;
        }

        // identität = identische reste / länge des alignments
        return (double) matches / alignedLength;
    }
}
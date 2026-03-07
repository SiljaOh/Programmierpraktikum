package de.lmu.ifi.bio.gor.cli;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.io.SeclibParser;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;


public class TestSetCreator {


    public static void main(String[] args) throws Exception {


        List<Protein> allProteins = SeclibParser.parseSeclib("cb513.db");


        //80/20 split
        List<Protein> trainingSet    = allProteins.subList(0, 400);

        List<Protein> potentialTestSet = allProteins.subList(400, allProteins.size());

        List<Protein> cleanTestSet = new ArrayList<>();

        for (Protein testProt : potentialTestSet) {
            boolean isTooSimilar = false;

            for (Protein trainProt : trainingSet) {
                double identity = calculateAlignmentIdentity(
                        testProt.getSequence(), trainProt.getSequence()
                );

                System.out.println(testProt.getId() + " vs " + trainProt.getId() + " = " + identity); // NEU
                //test proteine die mit irgendeinem Train protein eine Identität größer 25%
                // haben werden ausgeschlossen
                if (identity > 0.30) {
                    isTooSimilar = true;
                    break;
                }
            }

            if (!isTooSimilar) {
                cleanTestSet.add(testProt);
            }

        }
        try (PrintWriter pw = new PrintWriter("clean_test.db")) {
            for (Protein p : cleanTestSet) {
                pw.println("> " + p.getId());
                pw.println("AS " + p.getSequence());
                pw.println("SS " + p.getStructure());
                pw.println();
            }

        }
        try (PrintWriter pw = new PrintWriter("train_400.db")) {
            for (Protein p : trainingSet) {
                pw.println("> " + p.getId());
                pw.println("AS " + p.getSequence());
                pw.println("SS " + p.getStructure());
                pw.println();
            }
        }
        System.out.println("Trainingset: " + trainingSet.size());
        System.out.println("clean test set size: " + cleanTestSet.size());

    }

    private static double calculateAlignmentIdentity(String seq1, String seq2) {
        int n = seq1.length();
        int m = seq2.length();
        int GAP = -1;

        // Score-Matrix aufbauen
        int[][] F = new int[n + 1][m + 1];
        for (int i = 0; i <= n; i++) F[i][0] = i * GAP;
        for (int j = 0; j <= m; j++) F[0][j] = j * GAP;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int match = F[i-1][j-1] + (seq1.charAt(i-1) == seq2.charAt(j-1) ? 1 : -1);
                int del   = F[i-1][j] + GAP;
                int ins   = F[i][j-1] + GAP;
                F[i][j]   = Math.max(match, Math.max(del, ins));
            }
        }

        // Traceback
        StringBuilder a1 = new StringBuilder();
        StringBuilder a2 = new StringBuilder();
        int i = n, j = m;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0 && F[i][j] == F[i-1][j-1] + (seq1.charAt(i-1) == seq2.charAt(j-1) ? 1 : -1)) {
                a1.append(seq1.charAt(i-1));
                a2.append(seq2.charAt(j-1));
                i--; j--;
            } else if (i > 0 && F[i][j] == F[i-1][j] + GAP) {
                a1.append(seq1.charAt(i-1));
                a2.append('-');
                i--;
            } else {
                a1.append('-');
                a2.append(seq2.charAt(j-1));
                j--;
            }
        }

        String aln1 = a1.reverse().toString();
        String aln2 = a2.reverse().toString();

        int matches = 0;
        for (int k = 0; k < aln1.length(); k++) {
            if (aln1.charAt(k) != '-' && aln1.charAt(k) == aln2.charAt(k)) matches++;
        }
        return (double) matches / aln1.length();
    }
}

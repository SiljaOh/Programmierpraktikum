package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.io.SeclibParser;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class TestSetCreator {

    public static void main(String[] args) throws Exception {

        // 1. Pfade definieren (Passe diese ggf. an deine Ordnerstruktur an)
        String trainingPath = "C:/Users/gross/IdeaProjects/progprakt-17/tools/CB513DSSP.db";
        String inputPath = "C:/Users/gross/IdeaProjects/progprakt-17/tools/casp_all.seclib";
        String outputPath = "clean_test_casp.db";

        // 2. Daten laden
        List<Protein> trainingSet = SeclibParser.parseSeclib(trainingPath);
        List<Protein> potentialTestSet = SeclibParser.parseSeclib(inputPath);
        List<Protein> cleanTestSet = new ArrayList<>();

        double identityThreshold = 0.3; // 25% Grenze
        int minLength = 18;             // Mindestlänge für GOR-Kompatibilität

        System.out.println("Starte Filterung...");
        System.out.println("Trainings-Set (CB513) Größe: " + trainingSet.size());
        System.out.println("Potential Test-Set Größe: " + potentialTestSet.size());
        System.out.println("---------------------------------------------------");

        int processedCount = 0;

        // 3. Haupt-Filter-Schleife
        for (Protein testProt : potentialTestSet) {
            processedCount++;

            // Fortschrittsanzeige alle 10 Proteine
            if (processedCount % 10 == 0 || processedCount == potentialTestSet.size()) {
                System.out.print("\rVerarbeite Protein " + processedCount + " von " + potentialTestSet.size() + "...");
            }

            // A. Längen-Filter: Alles unter 30 Aminosäuren fliegt sofort raus
            if (testProt.getSequence().length() < minLength) {
                continue;
            }

            boolean isTooSimilar = false;
            String seqTest = testProt.getSequence();

            // B. Homologie-Check gegen jedes Protein im Trainings-Set
            for (Protein trainProt : trainingSet) {
                double identity = calculateAlignmentIdentity(seqTest, trainProt.getSequence());

                if (identity > identityThreshold) {
                    isTooSimilar = true;
                    break; // Zu ähnlich -> Abbruch der inneren Schleife
                }
            }

            // Wenn kein Protein im Training zu ähnlich war, ab ins Test-Set
            if (!isTooSimilar) {
                cleanTestSet.add(testProt);
            }
        }

        System.out.println("\n---------------------------------------------------");
        System.out.println("Filterung abgeschlossen!");

        // 4. Ergebnis speichern
        try (PrintWriter pw = new PrintWriter(outputPath)) {
            for (Protein p : cleanTestSet) {
                pw.println("> " + p.getId());
                pw.println("AS " + p.getSequence());
                pw.println("SS " + p.getStructure());
                pw.println();
            }
        }

        System.out.println("Sauberes Test-Set gespeichert in: " + outputPath);
        System.out.println("Finale Test-Set Größe: " + cleanTestSet.size());
    }

    /**
     * Berechnet die globale Sequenzidentität mittels Needleman-Wunsch Alignment.
     */
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

        // Traceback zur Berechnung der Matches
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

        // Identität = Identische Reste / Länge des Alignments
        return (double) matches / alignedLength;
    }
}
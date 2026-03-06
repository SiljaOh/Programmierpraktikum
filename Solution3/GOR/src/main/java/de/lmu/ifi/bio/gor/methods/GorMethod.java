package de.lmu.ifi.bio.gor.methods;
import de.lmu.ifi.bio.gor.model.Alignment;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.model.GorMatrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public abstract class GorMethod {

    protected GorMatrix matrix;
    public GorMatrix getMatrix() {
        return this.matrix;
    }

    protected boolean useProbability = false;

    public void setUseProbability(boolean useProbability) {
        this.useProbability = useProbability;
    }

    public static List<Protein> readFasta(String filePath) throws IOException {
        List<Protein> proteins = new ArrayList<Protein>();
        String id = null;
        StringBuilder sequence = new StringBuilder();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {

                    if (id != null && sequence.length() > 0) {
                        proteins.add(new Protein(id, sequence.toString(), " ".repeat(sequence.length())));
                        sequence.setLength(0); // StringBuilder leeren für das nächste Protein
                    }

                    id = line.substring(1).trim();
                } else {
                    sequence.append(line.replaceAll("\\s+", ""));
                }
            }
            if (id != null && sequence.length() > 0) {
                proteins.add(new Protein(id, sequence.toString(), " ".repeat(sequence.length())));
            }
        }

        if (proteins.isEmpty()) {
            throw new IOException("Datei '" + filePath + "' enthält keinen gültigen FASTA-Header (startet nicht mit '>').");
        }
        return proteins;
    }
    public static String whichMethod(String modelPath) throws IOException {
        // Header um Methode zu bestimmen
        try (BufferedReader reader = new BufferedReader(new FileReader(modelPath))) {
            String line;

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                // Header finden mehr nicht
                if (line.startsWith("//")) {
                    //wie beim train gor rausfinden
                    switch (line) {
                        case "// Matrix4D":
                            return "GOR3";
                        case "// Matrix6D":
                            return "GOR4";
                        default:
                            return "GOR1";
                    }
                }
            }
        }
        // kein Header einfach gor1
        return "GOR1";
    }

    public abstract void train(List<Protein> data);
    public abstract String predict(String sequence);

    private String getPaddedWindow(String seq, int center) {
        int m = 8;
        StringBuilder sb = new StringBuilder();
        for (int j = center - m; j <= center + m; j++) {
            if (j < 0 || j >= seq.length()) sb.append("-");
            else sb.append(seq.charAt(j));
        }
        return sb.toString();
    }
    public abstract void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result);

    protected String slidingWindow(Protein protein, boolean isTraining) {
        String sequence = protein.getSequence();
        String structure = protein.getStructure();
        StringBuilder result = new StringBuilder();

        int windowSize = 17;
        int m = windowSize / 2;

        //prediction fix oder so
        if (!isTraining && sequence.length() < windowSize) {
            return "C".repeat(sequence.length());
        }

        for (int i = m; i < sequence.length() - m; i++) {

            // Fenster ausschneiden (17 Zeichen lang)
            // substring(start, end) -> end ist exklusiv, daher + m + 1
            String window = sequence.substring(i - m, i + m + 1);

            char target = ' ';
            if (isTraining) {
                target = structure.charAt(i);
            }

            handleWindow(window, i, target, isTraining, result);
        }

        // predict: ignorierte Ränder mit Platzhaltern auffüllen
        if (!isTraining) {
            return "-".repeat(m) + result.toString() + "-".repeat(m);
        }

        return result.toString();
    }

    protected abstract double[] calculateRawScores(String window);

    public double[][] getScores(String sequence) {
        int len = sequence.length();
        double[][] scores = new double[len][3];
        int m = 8;

        for (int i = 0; i < len; i++) {

            String window = getPaddedWindow(sequence, i);

            double[] p = calculateRawScores(window);
            scores[i][0] = p[0]; // Score für H
            scores[i][1] = p[1]; // Score für E
            scores[i][2] = p[2]; // Score für C
        }
        return scores;
    }

    //GOR 5
    public String predictGor5(Alignment alignment) {
        String querySeq = alignment.getReferenceSeq();
        //länge mit gaps
        int queryLen = querySeq.length();
        //Summen arrays für H E C Scores pro Position
        double[] sumH = new double[queryLen];
        double[] sumE = new double[queryLen];
        double[] sumC = new double[queryLen];

        //Sequenz plus alle homologen
        List<String> allSeqs = new ArrayList<>();
        allSeqs.add(querySeq);
        allSeqs.addAll(alignment.getHomologues());

        //zu kleine sequenzen filtern
        for (String seq : allSeqs) {
            String seqWithoutGaps = seq.replace("-", "");
            if (seqWithoutGaps.length() < 17) continue;
            //[position][0/1/2]  hier score
            //für jede Position der gap-freien Sequenz drei Scores.
            double[][] readySeq = getScores(seqWithoutGaps);

            //wie posi ohne gaps
            //seq:       H  P  -  T  F  -  N  K
            //seqPos:    0  1     2  3     4  5
            int seqPos = 0;
            for (int i = 0; i < queryLen; i++) {
                //gaps in allen normalen seq skippen
                if (seq.charAt(i) == '-') continue;
                //haupt seq zählen
                if (querySeq.charAt(i) == '-') {
                    seqPos++;
                    continue; }

                //Score dieser Sequenz an dieser Position zur Summe addieren
                sumH[i] += readySeq[seqPos][0];
                sumE[i] += readySeq[seqPos][1];
                sumC[i] += readySeq[seqPos][2];
                seqPos++;
            }
        }

        StringBuilder result = new StringBuilder();
        for (int i = 0; i < queryLen; i++) {
            //gaps in query überspringen
            if (querySeq.charAt(i) == '-') continue;
            //für jede posi in sequenz höchster score
            if (sumH[i] >= sumE[i] && sumH[i] >= sumC[i]) {
                result.append('H');
            }
            else if (sumE[i] >= sumH[i] && sumE[i] >= sumC[i]) {
                result.append('E');
            }
            else result.append('C');
        }
        return result.toString();
    }


    //BONUS: postprocess, noch mehr ausprobieren!
    //mind 4 H, 2 E, C nur zwischne E und H
    public String postProcess(String prediction) {
        int minH = 4;
        int minE = 2;
        int n = prediction.length();
        char[] processed = prediction.toCharArray();

        // 1 mal vor 1 mal nach c korrektur
        for (int pass = 0; pass < 2; pass++) {

            // einzelnes C zwischen gleichen Nachbarn korrigieren
            // z.B. HHHCHHH -> HHHHHHH, aber HHHCEEE bleibt
            if (pass == 1) {
                for (int i = 1; i < n - 1; i++) {
                    if (processed[i] == 'C'
                            && processed[i - 1] == processed[i + 1]
                            && processed[i - 1] != 'C') {
                        processed[i] = processed[i - 1];
                    }
                }
            }

            for (int i = 0; i < n; i++) {
                char current = processed[i];

                // wv gleiche buchstaben links nebeneinander
                int left = 0;
                for (int j = i - 1; j >= 0; j--) {
                    if (processed[j] == current) left++;
                    else break;
                }

                // wv gleiche buchstaben rechts nebeneinander
                int right = 0;
                for (int j = i + 1; j < n; j++) {
                    if (processed[j] == current) right++;
                    else break;
                }

                // gesamt window anschauen
                int sameState = left + 1 + right;

                if (current == 'H' && sameState < minH) processed[i] = 'C';
                if (current == 'E' && sameState < minE) processed[i] = 'C';
            }
        }

        return new String(processed);
    }
}

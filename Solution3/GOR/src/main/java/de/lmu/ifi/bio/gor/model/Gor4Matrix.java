package de.lmu.ifi.bio.gor.model;

import java.io.*;

public class Gor4Matrix extends Gor3Matrix {

    //6D-Matrix
    //3*20*20*17*20*17 = 6.936.000 parameters
    //use all possible pairs in the window and sum their probabilities
    //counts[Zustand][Zentrum AS][erster Nachbar][erste nachbar posi][zweiter nachbar][zweiter nachbar posi]
    protected final int[][][][][][] counts4 = new int[3][20][20][17][20][17];
    protected final double[][][][][][] scores4 = new double[3][20][20][17][20][17];

     int[][] totalCountsSX = new int[3][20];


    //train
    public void addCount4(char state, char centralAA, char neighborOneAA, int neighborOnePos, char neighborTwoAA, int neighborTwoPos) {
        int sIdx = getStateIndex(state);
        int zIdx = getAaIndex(centralAA);
        int a1Idx = getAaIndex(neighborOneAA);
        int a2Idx = getAaIndex(neighborTwoAA);

        //Wie oft kamen Zustand S, Zentrum-AS, und das Paar (Rₖ an k, Rₗ an l) zusammen im Training vor counten
        //ür jede Kombination aus: Zustand S (H, E, C), Zentrum-AS Rⱼ, Paar (Rₖ an Position k, Rₗ an Position l)
        if (sIdx != -1 && zIdx != -1 && a1Idx != -1 && a2Idx != -1) {
            counts4[sIdx][zIdx][a1Idx][neighborOnePos][a2Idx][neighborTwoPos]++;

            //wie oft kam Zentrum-AS z in Zustand s vor, Kombi genau 1 mal pro fenster! -> erstes paar überhaupt
            if (neighborOnePos == 0 && neighborTwoPos == 1) {
                totalCountsSX[sIdx][zIdx]++;
            }
        }
    }

    //train
    @Override
    public void save(String path) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {
            pw.println("// Matrix6D");

            String aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
            char[] states = {'C', 'E', 'H'};
            //dim 1
            for (char s : states) {
                //dim2
                for (int z = 0; z < 20; z++) {
                    char centralAA = aminoAcids.charAt(z);
                    //dim 3
                    for (int n1 = 0; n1 < 20; n1++) {
                        //dim4
                        for (int k = 0; k < 17; k++) {
                            char neighbourOne = aminoAcids.charAt(n1);
                            //k=0 -8 ist links k=8 (mitte) ist dann 0
                            int posi = k - 8;
                            pw.println("=" + s + "," + centralAA + "," + neighbourOne + "," + posi + "=");
                            //alle mögichen nachbar AS in buchstabe
                            //dim5
                            for (int n2 = 0; n2 < 20; n2++) {
                                char neighbourTwo = aminoAcids.charAt(n2);
                                //anfang zeile
                                pw.print(neighbourTwo + "\t");
                                //dim6
                                for (int l = 0; l < 17; l++) {
                                    //korrekten counnt für jede zeile pro mögliche as
                                    int sidx = getStateIndex(s);
                                    pw.print(counts4[sidx][z][n1][k][n2][l] + "\t");

                                }
                                //neue zeile nue as
                                pw.println();
                            }
                            //neuer header -> neuer möglicher erster nachbar
                            pw.println();
                        }
                    }
                }
            }
            //GOR3 hinten dran hängen!!
            pw.println();
            pw.println("// Matrix4D");
            pw.println();

            char[] states3 = {'C', 'E', 'H'};
            for (int x = 0; x < aminoAcids.length(); x++) {
                char centralAA = aminoAcids.charAt(x);
                for (char s : states3) {
                    pw.println("=" + centralAA + "," + s + "=");
                    for (int y = 0; y < aminoAcids.length(); y++) {
                        char neighborAA = aminoAcids.charAt(y);
                        pw.print(neighborAA + "\t");
                        for (int m = 0; m < 17; m++) {
                            pw.print(counts3[getStateIndex(s)][m][x][y] + "\t");
                        }
                        pw.println();
                    }
                    pw.println();
                }
            }
        }
    }

    //predict)
    public void calculateTotalsFromCounts4() {
        // 1. Reset der geerbten Gor3-Felder (für P(S,X) und P(S))
        for (int s = 0; s < 3; s++) {
            java.util.Arrays.fill(totalCountsSX[s], 0);
            totalCountsPerState3[s] = 0;
        }

        // 2. Totals aus dem 4D-Teil (GOR III) extrahieren
        // Da loadModel4 die counts3 befüllt hat, nutzen wir diese für die Basis-Statistik
        for (int s = 0; s < 3; s++) {
            for (int z = 0; z < 20; z++) {
                // Wir schauen uns die Zentrum-Position (8) im counts3 an
                for (int a = 0; a < 20; a++) {
                    int c = counts3[s][8][z][a];
                    if (c > 0) {
                        totalCountsSX[s][z] += c;
                        totalCountsPerState3[s] += c;
                    }
                }
            }
        }
    }

    //predict
    //header: Zustand, Zentrum, FixerNachbar, Offset
     //=H,Y,M-8= heisst Y ist im Zentrum und an pos -8 ist ein M
    public void loadModel4(String filePath) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;

        // Hilfsvariablen für 6D
        int sIdx = -1, zIdx = -1, n1Idx = -1, p1Idx = -1;
        // Hilfsvariablen für 4D
        int gor3Zustand = -1, gor3Zentrum = -1;

        boolean isMatrix4DSection = false;

        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;

            // Sektions-Erkennung
            if (line.startsWith("// Matrix4D")) {
                isMatrix4DSection = true;
                continue;
            }

            if (line.startsWith("=")) {
                String content = line.replace("=", "");
                String[] parts = content.split(",");

                if (!isMatrix4DSection && parts.length == 4) {
                    // 6D Header parsen: Zustand, Zentrum, Nachbar1, Pos1
                    sIdx = getStateIndex(parts[0].charAt(0));
                    zIdx = getAaIndex(parts[1].charAt(0));
                    n1Idx = getAaIndex(parts[2].charAt(0));
                    p1Idx = Integer.parseInt(parts[3]) + 8;
                } else if (isMatrix4DSection && parts.length == 2) {
                    // 4D Header parsen: Zentrum, Zustand
                    gor3Zentrum = getAaIndex(parts[0].charAt(0));
                    gor3Zustand = getStateIndex(parts[1].charAt(0));
                }
            } else {
                String[] values = line.split("\\s+");
                if (values.length == 18) {
                    int n2Idx = getAaIndex(values[0].charAt(0));
                    if (n2Idx == -1) continue;

                    if (!isMatrix4DSection && sIdx != -1) {
                        // In counts4 speichern
                        for (int p2Idx = 0; p2Idx < 17; p2Idx++) {
                            counts4[sIdx][zIdx][n1Idx][p1Idx][n2Idx][p2Idx] = Integer.parseInt(values[p2Idx + 1]);
                        }
                    } else if (isMatrix4DSection && gor3Zustand != -1) {
                        // In counts3 speichern (geerbt von Gor3Matrix)
                        for (int pIdx = 0; pIdx < 17; pIdx++) {
                            counts3[gor3Zustand][pIdx][gor3Zentrum][n2Idx] = Integer.parseInt(values[pIdx + 1]);
                        }
                    }
                }
            }
        }
        reader.close();

        // Wichtig: Da die Datei nun beide Daten enthält, musst du NICHTS mehr "retten".
        // Du berechnest einfach beide Scores.
        calculateTotalsFromCounts4();
        calculateScores4();
    }

    //predict
    //formel:
    //I = (2/(2m+1)) * sum_k sum_l>k log( f(S, a_j+k, a_j+l) / f(¬S, a_j+k, a_j+l) ) 6D Array
    // - ((2m-1)/(2m+1)) * sum_k log( f(S, a_j+k) / f(¬S, a_j+k) ) 4D Array aus Gor3
    public void calculateScores4() {

        double pseudoCount = 0.5;

        for (int s = 0; s < 3; s++) {
            for (int z = 0; z < 20; z++) { // ZentrumAS
                for (int n1 = 0; n1 < 20; n1++) { //Nachbar1AS an Pos k
                    for (int p1 = 0; p1 < 17; p1++) {//PosK
                        for (int n2 = 0; n2 < 20; n2++) { //Nachbar2AS an Pos l
                            for (int p2 = 0; p2 < 17; p2++) {//PosL

                                // f(S, a_k, a_l)
                                double fS = counts4[s][z][n1][p1][n2][p2] + pseudoCount;

                                // f(¬S, a_k, a_l)
                                double fNotS = 2* pseudoCount;
                                for (int o = 0; o < 3; o++) {
                                    if (o != s) {
                                        fNotS += counts4[o][z][n1][p1][n2][p2];
                                    }
                                }

                                // log( f(S) / f(¬S) )
                                scores4[s][z][n1][p1][n2][p2] = Math.log(fS / fNotS);
                            }
                        }
                    }
                }
            }
        }
        calculateScores3();
    }


    //predict
    // predict
    public char getBestState4(String window) {
        int m = 8;
        int N = 17;
        double[] totals = new double[3];

        // PERFORMANCE-BOOSTER:
        // Indizes des gesamten Fensters vorab einmal berechnen
        int[] windowIndices = new int[N];
        for (int i = 0; i < N; i++) {
            windowIndices[i] = getAaIndex(window.charAt(i));
        }

        // Das Zentrum aus dem vorbereiteten Array holen
        int centerIdx = windowIndices[m];
        if (centerIdx == -1) centerIdx = 0; // Fallback für unbekannte AS

        for (int s = 0; s < 3; s++) {
            double sum6D = 0;
            double sum4D = 0;

            for (int k = 0; k < N; k++) {
                int aaK = windowIndices[k]; // Schnellzugriff aufs Array
                if (aaK == -1) continue;

                // Der GOR 3 Teil (Minus-Term)
                sum4D += scores3[s][k][centerIdx][aaK];

                // Der GOR 4 Teil (6D Doppelsumme l > k)
                for (int l = k + 1; l < N; l++) {
                    int aaL = windowIndices[l]; // Schnellzugriff aufs Array
                    if (aaL == -1) continue;

                    sum6D += scores4[s][centerIdx][aaK][k][aaL][l];
                }
            }
            // Formel: (2/17) * sum6D - (15/17) * sum4D
            totals[s] = (2.0 / N) * sum6D - ((2.0 * m - 1.0) / N) * sum4D;
        }

        // Winner-Logik
        if (totals[0] >= totals[1] && totals[0] >= totals[2]) return 'H';
        if (totals[1] >= totals[0] && totals[1] >= totals[2]) return 'E';
        return 'C';
    }

}

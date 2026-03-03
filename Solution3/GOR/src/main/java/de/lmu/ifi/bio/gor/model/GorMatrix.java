package de.lmu.ifi.bio.gor.model;

import java.io.*;

public class GorMatrix {

    //3 observing states H,E,C
    //window : paddding + zentrum + padding = 8+1+8 = 17
    //20 amino acids
    //Idee für jede AS deren Struktur man predicten werden müssen 3 * 17 * 20 = 1020 Parameter geschätzt werden
    protected final int[][][] counts = new int[3][17][20];
    protected final double[][][] scores = new double[3][17][20];

    int[] totalCountsperState = new int[3];
    int totalAA = 0;


    // Hilfsmethode um 20 Zahlen in Matrix die Aminosäuren zuzuordnen
    //zb. A ist Index 0 in Matrix
    public int getAaIndex(char aa) {
        String aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
        return aminoAcids.indexOf(aa);
    }

    // Hilfsmethode um 3 Zahlen in Matrix den 3 Zuständen zuzuordnen
    //zb. H ist Index 0 in Matrix
    public int getStateIndex(char ss) {
        return switch (ss) {
            case 'H' -> 0;
            case 'E' -> 1;
            default -> 2; //Alles andere ist Coil
        };
    }

    //Matrix befüllen
    public void addCount(char state, int windowPos, char aa){
        int stateIndex = getStateIndex(state);
        int aaIndex = getAaIndex(aa);

        if (aaIndex == -1) {
            return;
        }

        counts[stateIndex][windowPos][aaIndex]++;

        //windowPos muss 8 sein sonst zählen wir zu oft
        if(windowPos == 8){
            totalCountsperState[stateIndex]++;
            totalAA++;
        }
    }

    public int getCount(char state, int windowPos, char aa) {
        return counts[getStateIndex(state)][windowPos][getAaIndex(aa)];
    }

    //predict
    public void calculateScores() {

        int[] aminoAcidTotals = new int[20];
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                for (int a = 0; a < 20; a++) {

                    aminoAcidTotals[a] += counts[s][p][a];
                }
            }
        }
        //H E C
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                for (int a = 0; a < 20; a++) { //für alle 20 Aminosäuren

                    double countSR = counts[s][p][a];
                    double countS = totalCountsperState[s];
                    double countR = aminoAcidTotals[a];

                    if (countSR == 0  || countS == 0  || countR == 0) {
                        // Strafe für 0
                        scores[s][p][a] = -5.0;
                    } else {
                        // Score = ln( P(Struktur, Aminosäure) / (P(Struktur) * P(Aminosäure)) )
                        scores[s][p][a] = Math.log((countSR * totalAA) / (countS * countR));
                    }
                }
            }
        }
    }



    public int getTotalCountsPerState(char state) {
        return totalCountsperState[getStateIndex(state)];
    }

    public int getTotalAA() {
        return totalAA;
    }


    // hier zeilen AS und spalten window
    public void loadModel(String filePath) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;
        //0 H, 1 E, 2 C
        int stateIdx = -1;

        //gegeben lesen
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            // wenn Zeile mit "=" beginnt, haben wir einen neuen Zustand
            if (line.startsWith("=")) {
                //zwischen = zustand ss finden
                char stateChar = line.charAt(1);
                //als zahl speichern
                stateIdx = getStateIndex(stateChar);
            }
            //zeile AS mit 17 Positionswete
            else if (stateIdx != -1) {
                String[] values = line.split("\\s+");
                // Hier prüfen, ob wir genug Daten haben, um zu parsen
                if (values.length == 18) {
                    int aa = getAaIndex(values[0].charAt(0));
                    for (int pos = 0; pos < 17; pos++) {
                        counts[stateIdx][pos][aa] = Integer.parseInt(values[pos + 1]);
                    }
                }
            }
        }
        calculateTotalsFromCounts();
        calculateScores();
    }

    public char getBestState(String window) {
        //O H E
        double[] totals = new double[3];

        //SUmme
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                char aa = window.charAt(p);
                int aaIdx = getAaIndex(aa);
                //alle von 1 Zustand addieren
                totals[s] += scores[s][p][aaIdx];
            }
        }

        // welcher Zustand hat den höchsten Score
        if (totals[0] >= totals[1] && totals[0] >= totals[2]) {
            return 'H';
        }
        if (totals[1] >= totals[0] && totals[1] >= totals[2]) {
            return 'E';
        }
        return 'C';
    }

    public void calculateTotalsFromCounts() {
        totalCountsperState = new int[3];
        totalAA = 0;
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                for (int a = 0; a < 20; a++) {
                    if (p == 8) { // Nur im Zentrum zählen wir für die Statistik
                        totalCountsperState[s] += counts[s][p][a];
                        totalAA += counts[s][p][a];
                    }
                }
            }
        }
    }

    public void save(String path) throws IOException {
        PrintWriter pw = new PrintWriter(new FileWriter(path));

        pw.println("// Matrix3D");
        pw.println(); // Leere Zeile nach Header

        char[] states = {'C', 'E', 'H'};

        String aminoAcids = "ACDEFGHIKLMNPQRSTVWY";

        for (char s : states) {
            pw.println("=" + s + "=");
            pw.println("\t");

            for (int a = 0; a < aminoAcids.length(); a++) {
                char aa = aminoAcids.charAt(a);
                pw.print(aa + "\t"); //aminosäure + Tab

                for (int p = 0; p < 17; p++) {
                    //counts!
                    pw.print(this.counts[getStateIndex(s)][p][getAaIndex(aa)] + "\t");
                }
                pw.println();
            }
            pw.println();
            pw.println();
        }
        pw.close();
    }

}
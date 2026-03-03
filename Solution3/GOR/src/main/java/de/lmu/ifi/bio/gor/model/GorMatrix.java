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

    //train
    // Hilfsmethode um 20 Zahlen in Matrix die Aminosäuren zuzuordnen
    //zb. A ist Index 0 in Matrix
    public int getAaIndex(char aa) {
        String aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
        return aminoAcids.indexOf(aa);
    }

    //train
    // Hilfsmethode um 3 Zahlen in Matrix den 3 Zuständen zuzuordnen
    //zb. H ist Index 0 in Matrix
    public int getStateIndex(char ss) {
        return switch (ss) {
            case 'H' -> 0;
            case 'E' -> 1;
            default -> 2; //Alles andere ist Coil
        };
    }
    //train
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

    //train
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

    //train und predict
    public void calculateTotalsFromCounts() {
        //Wahrscheinlichkeit für Helix = totalCountsperState[0]/totalAA
        //für jeden state wie oft AS, für P(S)
        totalCountsperState = new int[3];
        //alle AS die ma nim zentrum sieht
        totalAA = 0;
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                for (int a = 0; a < 20; a++) {
                    // Nur im Zentrum AS
                    if (p == 8) {
                        //wie oft as  an posi 8 wenn zustand s
                        totalCountsperState[s] += counts[s][p][a];
                        //wert egal welcher zustand addiert
                        totalAA += counts[s][p][a];
                    }
                }
            }
        }
    }

    //predict
    // hier zeilen AS und spalten window
    public void loadModel(String filePath) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;
        //0 H, 1 E, 2 C
        int stateIdx = -1;

        //gegebene file lesen
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            // wenn Zeile mit "=" beginnt, haben wir einen neuen Zustand
            if (line.startsWith("=")) {
                //zwischen = zustand ss finden, zb H
                char stateChar = line.charAt(1);
                //als zahl speichern (0 1 oder 2)
                stateIdx = getStateIndex(stateChar);
            }
            //zeile AS mit 17 Positionswete
            else if (stateIdx != -1) {
                //leerzeichen etc weg
                String[] values = line.split("\\s+");
                //AS plus 17er Window
                if (values.length == 18) {
                    //AS als zahl speichern
                    int aa = getAaIndex(values[0].charAt(0));
                    for (int pos = 0; pos < 17; pos++) {
                        counts[stateIdx][pos][aa] = Integer.parseInt(values[pos + 1]);
                    }
                }
            }
        }
        calculateTotalsFromCounts();
        calculateScores1();
    }

    //predict VERSCHIEBENNN ODER ???
    public void calculateScores1() {

        for (int s = 0; s < 3; s++) {

            //Formel 2. Teil immer gleich, log( f_s̄ / f_s ), wie viel zustand s
            int fS = totalCountsperState[s];
            //wie viele nicht zustand s
            int f_S = totalAA - fS;

            for (int p = 0; p < 17; p++) {
                for (int a = 0; a < 20; a++) {

                    //1. Teil Formel
                    //wie oft as a an posi p wenn  zustand s
                    double fSA    = counts[s][p][a];
                    //wie oft a an p wenn nicht zustand
                    double f_SA = 0;
                    for (int o = 0; o < 3; o++) {
                        if (o != s){
                            f_SA += counts[o][p][a];
                        }
                    }

                    // Formel = I = log( f_s,a / f_s̄,a ) + log( f_s̄ / f_s )
                    scores[s][p][a] = Math.log(fSA / f_SA) + Math.log((double) f_S / fS);
                }
            }
        }
    }

    

    //predict
    public char getBestState(String window) {
        double[] totals = new double[3];
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                char aa = window.charAt(p);
                int aaIdx = getAaIndex(aa);
                if (aaIdx == -1) continue;
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
}
package de.lmu.ifi.bio.gor.model;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Gor3Matrix extends GorMatrix {

    //3 observing states H,E,C
    //window: paddding + zentrum + padding = 8+1+8 = 17
    //20 aminoacids : Zentrum
    //20 aminoacids : Nachbar
    //Idee für jede AS deren Struktur man predicten werden müssen 3 * 17 * 20 *20 = 20400 Parameter geschätzt werden
    protected final int[][][][] counts3 = new int[3][17][20][20];
    protected final double[][][][] scores3 = new double[3][17][20][20];

    //gor3 hat zusätzlcih; Wie oft tritt Zentrum X in Zustand S auf?
    //P(S,X)
    protected final int[][] totalCountsSX = new int[3][20];

    private final int[] totalCountsPerState3 = new int[3];

    //train
    // Eigene addCount Methode für GOR3
    public void addCount3(char state, int windowPos, char centralAA, char neighborAA) {
        int sIdx = getStateIndex(state);
        int xIdx = getAaIndex(centralAA); // Residuum X an Position j
        int yIdx = getAaIndex(neighborAA); // Residuum Y an Position j+m

        if (sIdx != -1 && xIdx != -1 && yIdx != -1) {
            counts3[sIdx][windowPos][xIdx][yIdx]++;

            //nur einmal pro fenster zentrum zählen
            if (windowPos == 8) {
                totalCountsSX[sIdx][xIdx]++;
            }
        }
    }

    //train
    @Override
    public void save(String path) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {
            pw.println("// Matrix4D");

            String aminoAcids = "ACDEFGHIKLMNPQRSTVWY";
            char[] states = {'C', 'E', 'H'};

            for (int x = 0; x < aminoAcids.length(); x++) {
                char centralAA = aminoAcids.charAt(x);
                for (char s : states) {
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

    //train (and predict)
    protected void calculateTotals3() {
        java.util.Arrays.fill(totalCountsPerState3, 0);
        for (int s = 0; s < 3; s++) {
            for (int x = 0; x < 20; x++) {
                // Wir nutzen die Werte am Zentrum (Index 8),
                // um die Gesamtzahl der Zustände S für Zentrum X zu erhalten
                for (int y = 0; y < 20; y++) {
                    totalCountsPerState3[s] += counts3[s][8][x][y];
                }
            }
        }
    }

    //predict
    public void load3(String path) throws IOException {}

    //predict
    public void calculateScores3(){}

    //predict
    public char getBestState3(String window) {
        return window.charAt(0);
    }



}

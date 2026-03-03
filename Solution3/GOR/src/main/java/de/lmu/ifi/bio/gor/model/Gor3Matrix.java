package de.lmu.ifi.bio.gor.model;

import java.io.*;

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
    //s = 3 , x = 20
    protected final int[][] totalCountsSX = new int[3][20];

    int[] totalCountsPerState3 = new int[3];

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
    public void calculateTotalsFromCounts3() {
        //Wahrscheinlichkeit für Helix = totalCountsperState[0]/totalAA
        //für jeden state wie oft AS, für P(S)
        totalCountsPerState3 = new int[3];

        for (int s = 0; s < 3; s++) {
            //zentrum as
            for (int z = 0; z < 20; z++) {
                //nachbar as
                for (int a = 0; a < 20; a++) {
                    totalCountsPerState3[s] += counts3[s][8][z][a];
                    //wie oft kam Zentrum-AS z in Zustand s vor, egal welcher Nachbar
                    totalCountsSX[s][a] += counts3[s][8][z][a];
                }
            }
        }
    }

    //predict
    //Header: Zentrum AS und Zustand
    //Zeile: Nachbar AS + 17 Werte -> für jede 20 zentrum AS eine GOR! MATRIX
    public void loadModel3(String filePath) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;
        //0 H, 1 E, 2 C
        int stateIdx = -1;
        //zentrum as 20
        int zentrumAS = -1;

        //gegebene file lesen
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            // wenn Zeile mit "=" beginnt, haben wir einen neuen Header mit AS und Zustand
            if (line.startsWith("=")) {
                char stateChar = line.charAt(3);
                char centralAA = line.charAt(1);
                stateIdx = getStateIndex(stateChar);
                zentrumAS = getAaIndex(centralAA);
            }
            //zeile nachbar AS mit 17 Positionswete
            else if (stateIdx != -1 && zentrumAS != -1) {
                //leerzeichen etc weg
                String[] values = line.split("\\s+");
                //AS plus 17er Window
                if (values.length == 18) {
                    //AS als zahl speichern
                    int aa = getAaIndex(values[0].charAt(0));
                    if (aa == -1) continue;
                    for (int pos = 0; pos < 17; pos++) {
                        counts3[stateIdx][pos] [zentrumAS][aa] = Integer.parseInt(values[pos + 1]);
                    }
                }
            }
        }
        calculateTotalsFromCounts3();
        calculateScores3();
    }


    //predict
    //I = log( f(S, a_j, a_{j+k}) / f(¬S, a_j, a_{j+k}) ) + log( f(¬S, a_j) / f(S, a_j) )
    public void calculateScores3() {
        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                //Zentrum AS a_j
                for (int c = 0; c < 20; c++) {
                    // Nachbar-AS a_{j+k}
                    for (int a = 0; a < 20; a++) {

                        // log( f(S, a_j, a_{j+k}) / f(¬S, a_j, a_{j+k}) )
                        //wie oft zentrum c nachbar a paar an posi p zusand s
                        double fSA = counts3[s][p][c][a];
                        //wie oft gleiches par anderer zustand
                        double f_SA = 0;
                        for (int o = 0; o < 3; o++) {
                            if (o != s) {
                                f_SA += counts3[o][p][c][a];
                            }
                        }

                        // log( f(¬S, a_j) / f(S, a_j) )
                        //zentrum as c alleine in zustand
                        double fSX = totalCountsSX[s][c];
                        //wie oft in anderem
                        double f_SX = 0;
                        for (int o = 0; o < 3; o++) {
                            if (o != s) f_SX += totalCountsSX[o][c];  // f(¬S, a_j)
                        }
                        //nich bwi 1
                        if (fSA == 0 || f_SA == 0 || fSX == 0 || f_SX == 0) {
                            scores3[s][p][c][a] = 0.0;
                            continue;
                        }

                        scores3[s][p][c][a] = Math.log(fSA / f_SA) + Math.log(f_SX / fSX);
                    }
                }
            }
        }
    }


    public char getBestState3(String window) {
        // Zentrum-AS bestimmen posi 8 inn zahl
        int centerAS = getAaIndex(window.charAt(8));
        //unbekannt im zentrum
        if (centerAS == -1) centerAS = 0;
        //ein gesamtscore pro zustnd
        double[] totals = new double[3];

        for (int s = 0; s < 3; s++) {
            for (int p = 0; p < 17; p++) {
                //as aktuelle posi (=Nachbar)
                char aa = window.charAt(p);
                //nachbar in zahl
                int aaIdx = getAaIndex(aa);
                if (aaIdx == -1) continue;
                //Score für diesen Zustand, diese Position, diese Zentrum-AS und diese Nachbar-AS
                totals[s] += scores3[s][p][centerAS][aaIdx];
            }
        }

        if (totals[0] >= totals[1] && totals[0] >= totals[2]) return 'H';
        if (totals[1] >= totals[0] && totals[1] >= totals[2]) return 'E';
        return 'C';
    }

}

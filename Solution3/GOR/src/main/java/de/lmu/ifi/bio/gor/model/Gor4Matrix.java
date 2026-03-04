package de.lmu.ifi.bio.gor.model;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Gor4Matrix extends Gor3Matrix {

    //6D-Matrix
    //3*20*20*17*20*17 = 6.936.000 parameters
    //use all possible pairs in the window and sum their probabilities
    protected final int[][][][][][] counts4 = new int[3][20][20][17][20][17];
    protected final double[][][][][][] scores4 = new double[3][20][20][17][20][17];


    //train
    public void addCount4(char state, int windowPos, char centralAA, char neighborAA) {
    }

    //train
    @Override
    public void save(String path) throws IOException {}

    //train (and predict)
    public void calculateTotalsFromCounts4() {}

    //predict
    //header: Zustand, Zentrum, FixerNachbar, Offset
    // =H,Y,M-8= heisst Y ist im Zentrum und an pos -8 ist ein M
    public void loadModel4(String filePath) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filePath));
        String line;

        int sIdx = -1;
        int zIdx = -1;
        int n1Idx = -1;
        int p1Idx = -1;

        //gegebene file lesen
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty()) continue;
            // wenn Zeile mit "=" beginnt, haben wir einen neuen Header mit AS und Zustand
            if (line.startsWith("=")) {
                //header parsen
                String content = line.replace("=", "");
                String[] parts = content.split(",");

                if (parts.length == 4) {
                    sIdx = getStateIndex(parts[0].charAt(0));
                    zIdx = getAaIndex(parts[1].charAt(0));
                    n1Idx = getAaIndex(parts[2].charAt(0));
                    // Wichtig: Offset -8 bis 8 auf Index 0 bis 16 mappen
                    p1Idx = Integer.parseInt(parts[3]) + 8;
                }
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
        calculateTotalsFromCounts4();
        calculateScores4();
    }

    //predict
    public void calculateScores4() {

    }

    //predict
    public char getBestState4(String window) {
        return 0;
    }



}

package de.lmu.ifi.bio.gor.model;

import java.io.File;

public class GorMatrix {

    //3 observing states H,E,C
    //17 window: paddding + zentrum + padding = 8+1+8 = 17
    //20 amino acids
    //Idee für jede AS deren Struktur man predicten werden müssen 3 * 17 * 20 = 1020 Parameter geschätzt werden
    private final int[][][] matrix = new int[3][17][20];


    int[] totalCountsperState = new int[3];
    int totalAA;


    // Hilfsmethode um 20 Zahlen in Matrix die Aminosäuren zuzuordnen
    //zb. A ist Index 0 in Matrix
    public int getAaIndex(char aa) {
        String aminoAcids = "ARNDCQEGHILKMFPSTWYV";
        int index = aminoAcids.indexOf(aa);
        return (index == -1) ? 19 : index; //Fallback für unbekannte AS (z.B. 'X')
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
        matrix[stateIndex][windowPos][aaIndex]++;

        //windowPos muss 8 sein sonst zählen wir zu oft
        if(windowPos == 8){
            totalCountsperState[stateIndex]++;
            totalAA++;
        }
    }

    public void save(File file){} //matrix speichern???


}

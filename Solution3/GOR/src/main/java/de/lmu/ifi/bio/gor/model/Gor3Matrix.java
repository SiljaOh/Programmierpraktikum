package de.lmu.ifi.bio.gor.model;

public class Gor3Matrix extends GorMatrix {


    //3 observing states H,E,C
    //window: paddding + zentrum + padding = 8+1+8 = 17
    //20 amino acids
    //Idee für jede AS deren Struktur man predicten werden müssen 3 * 17 * 20 = 1020 Parameter geschätzt werden
    protected final int[][][][] countsGor3 = new int[3][17][20][20];
    protected final int[][][][] scoresGor3 = new int[3][17][20][20];



}

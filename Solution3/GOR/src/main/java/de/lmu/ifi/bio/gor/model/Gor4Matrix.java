package de.lmu.ifi.bio.gor.model;

import java.io.IOException;

public class Gor4Matrix extends Gor3Matrix {

    //6D-Matrix
    //3*20*20*17*20*17 = 6.936.000 parameters
    //use all possible pairs in the window and sum their probabilities
    protected final int[][][][][][] counts3 = new int[3][20][20][17][20][17];
    protected final double[][][][][][] scores3 = new double[3][20][20][17][20][17];


    //train
    public void addCount4(char state, int windowPos, char centralAA, char neighborAA) {
    }

    //train
    @Override
    public void save(String path) throws IOException {}

    //train (and predict)
    public void calculateTotalsFromCounts4() {}

    //predict
    public void loadModel4(String filePath) throws IOException {
    }

    //predict
    public void calculateScores4() {}

    //predict
    public char getBestState4(String window) {
        return 0;
    }



}

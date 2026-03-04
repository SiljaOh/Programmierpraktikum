package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.GorMatrix;
import de.lmu.ifi.bio.gor.model.Protein;


import java.util.List;

//KONSTANTEN

public class Gor1 extends GorMethod{

    public Gor1(){
        this.matrix=new GorMatrix();
    }

    @Override
    //calculate probabilities of observing state i = {H, E, C}
    // at the central window position when finding residue X at
    //position +-8 from central position in training set
    //aim: 3 matrices with 17*20 entries each
    public void train(List<Protein> data) {
        for(Protein p : data){
            this.slidingWindow(p,true);
        }
    }

    @Override
    public String predict(String sequence) {
        Protein protein = new Protein("ID", sequence, " ".repeat(sequence.length()));
        return slidingWindow(protein, false);
    }

    @Override
    public void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result) {
        if(isTraining) {
            for (int w = 0; w < window.length(); w++) {
                char aa = window.charAt(w);
                this.matrix.addCount(target, w, aa);
            }
        }
        else {
            char bestState = this.matrix.getBestState(window);
            result.append(bestState);
        }
    }

    //GorV:
    @Override
    protected double[] calculateRawScores(String window) {
        return this.matrix.getBestStateGorV(window);
    }

}
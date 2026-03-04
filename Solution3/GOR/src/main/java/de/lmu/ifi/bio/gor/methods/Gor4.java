package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.Gor4Matrix;
import de.lmu.ifi.bio.gor.model.Protein;

import java.util.List;

public class Gor4 extends GorMethod {


    public Gor4() {
        this.matrix=new Gor4Matrix();
    }

    public Gor4Matrix getMatrix4() {
        return (Gor4Matrix) this.matrix;
    }

    @Override
    public void train(List<Protein> data) {
        for(Protein p : data){
            slidingWindow(p,true);
        }
    }

    @Override
    public String predict(String sequence) {
        Protein protein = new Protein("ID", sequence, " ".repeat(sequence.length()));
        return slidingWindow(protein, false);
    }

    @Override
    public void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result) {
        if (!isTraining){
            char bestState = getMatrix4().getBestState4(window);
            result.append(bestState);
            return;
        }

        Gor4Matrix gor4Matrix = (Gor4Matrix) this.matrix;

        char centralAA = window.charAt(8);

        //Ohne l > k zählt man jedes Paar doppelt
        for(int k=0; k < 17; k++){
            //wid hier mit k+1 sichergestellt
            for (int l = k+1; l < 17; l++){
                char neighbourOneAA = window.charAt(k);
                char neighbourTwoAA = window.charAt(l);
                gor4Matrix.addCount4(target, centralAA, neighbourOneAA, k, neighbourTwoAA, l );
            }
        }

        //GOR3 noch befüllen für 2. teil
        for (int m = 0; m < 17; m++) {
            char neighbourAA = window.charAt(m);
            gor4Matrix.addCount3(target, m, centralAA, neighbourAA);
        }
    }
    //GorV:
    @Override
    protected double[] calculateRawScores(String window) {
        return this.matrix.getBestStateGorV(window);
    }
}

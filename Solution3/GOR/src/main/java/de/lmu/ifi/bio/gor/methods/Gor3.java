package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.Gor3Matrix;
import de.lmu.ifi.bio.gor.model.Protein;

import java.util.List;

public class Gor3 extends GorMethod{

    private final Gor3Matrix gor3Matrix = new Gor3Matrix();

    public Gor3() {
        this.matrix = gor3Matrix;
    }

    public Gor3Matrix getMatrix3() {
        return gor3Matrix;
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
            char bestState = getMatrix3().getBestState3(window);
            result.append(bestState);
            return;
        }
        char centralAA = window.charAt(8);

        //für das gesamte zentrum:
        for (int m = 0; m < 17; m++) {
            //Nachbar Y an Position j+m
            char neighborAA = window.charAt(m);
            //Paar infos (Zentrum, Nachbar, Position)
            gor3Matrix.addCount3(target, m, centralAA, neighborAA);
        }
    }

    //Gor5:
    @Override
    protected double[] calculateRawScores(String window, boolean probabilities) {
        return this.matrix.getBestStateGorV(window, probabilities);
    }
}

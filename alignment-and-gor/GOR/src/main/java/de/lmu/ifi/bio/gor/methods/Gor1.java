package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.GorMatrix;
import de.lmu.ifi.bio.gor.model.Protein;
import java.util.List;

public class Gor1 extends GorMethod {
    public Gor1() {
        this.matrix = new GorMatrix();
    }

    public void train(List<Protein> data) {
        for(Protein p : data) {
            this.slidingWindow(p, true);
        }

    }

    public String predict(String sequence) {
        Protein protein = new Protein("ID", sequence, " ".repeat(sequence.length()));
        return this.slidingWindow(protein, false);
    }

    public void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result) {
        if (isTraining) {
            for(int w = 0; w < window.length(); ++w) {
                char aa = window.charAt(w);
                this.matrix.addCount(target, w, aa);
            }
        } else {
            char bestState = this.matrix.getBestState(window);
            result.append(bestState);
        }

    }

    protected double[] calculateRawScores(String window) {
        return this.matrix.getBestStateGorV(window);
    }
}

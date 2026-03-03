package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.Gor3Matrix;
import de.lmu.ifi.bio.gor.model.Protein;

import java.util.List;

public class Gor3 extends GorMethod{

    public Gor3() {
        this.matrix=new Gor3Matrix();
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
        //matrix auf Gor3Matrix casten! (Gor3Matrix erbt von GorMatrix)
        Gor3Matrix gor3Matrix = (Gor3Matrix) this.matrix;

        char centralAA = window.charAt(8); //Anker X an pos j

        //für das gesamte zentrum:
        for (int m = 0; m < 17; m++) {
            char neighborAA = window.charAt(m); //Nachbar Y an Position j+m

            //Paar-Information (Zentrum, Nachbar, Position)
            gor3Matrix.addCount3(target, m, centralAA, neighborAA);
        }
    }
}

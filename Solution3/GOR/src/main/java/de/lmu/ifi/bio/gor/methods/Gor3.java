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
    }

    @Override
    public String predict(String sequence) {
        return "";
    }

    @Override
    public void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result) {

    }
}

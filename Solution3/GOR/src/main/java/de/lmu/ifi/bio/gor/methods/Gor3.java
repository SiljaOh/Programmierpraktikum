package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.Protein;

import java.util.List;

public class Gor3 extends GorMethod{


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

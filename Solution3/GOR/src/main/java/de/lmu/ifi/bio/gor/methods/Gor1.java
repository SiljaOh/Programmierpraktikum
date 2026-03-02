package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.Protein;

import java.util.List;

public class Gor1 extends GorMethod{
    @Override
    //calculate probabilities of observing state i = {H, E, C}
    // at the central window position when finding residue X at
    //position +-8 from central position in training set
    //aim: 3 matrices with 17*20 entries each
    public void train(List<Protein> data) {

    }

    @Override
    //calculate 3 probabilities for each residue to be predicted
    //for being helix H, sheet E or coil C
    //aim: predict state with max probability for each residue
    public String predict(String sequence) {
        return "";
    }

    @Override
    public void handleWindow(String window, int position, char targetStrucutre) {

    }
}

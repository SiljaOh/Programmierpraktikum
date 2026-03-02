package de.lmu.ifi.bio.gor.methods;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.model.GorMatrix;

import java.util.List;

abstract class GorMethod {
    //Vererbungsklasse für Gor1,3,4,5
    //alle Methoden die mehrmals benötigt werden

    protected GorMatrix matrix = new GorMatrix();

    public abstract void train(List<Protein> data);

    public abstract String predict(String sequence);

    public abstract void handleWindow(String window, int position, char targetStrucutre);


    private String applyPadding(String sequence) {
        String padding = "XXXXXXXX"; //17er Fenster; 8 links, Zentrum, 8 rechts
        return padding + sequence + padding;
    }

    public void slidingWindow(String sequence) {
        String paddedSequence = applyPadding(sequence);
        for(int i = 0; i < sequence.length(); i++) {
            String window = paddedSequence.substring(i, i + 17); //Window size = 17
            //handleWindow(window,i,structureAtPosition);
        }

    }









}

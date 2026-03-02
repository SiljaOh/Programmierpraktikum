package de.lmu.ifi.bio.gor.methods;
import de.lmu.ifi.bio.gor.model.Protein;

import java.util.List;

abstract class GorMethod {
    //Vererbungsklasse für Gor1,3,4,5

    public abstract void train(List<Protein> data);

    public abstract String predict(String sequence);









}

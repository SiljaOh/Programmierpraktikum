package de.lmu.ifi.bio.gor.model;

public class Protein {
    private String id; //PDB
    private String sequence; // Die AS Zeile
    private String structure; // Die SS Zeile

    public Protein(String id, String as, String ss) {
        this.id = id.trim();
        this.sequence = as.trim();
        this.structure = ss.trim();
    }

    public String getId() {
        return id;
    }

    public String getSequence() {
        return sequence;
    }

    public String getStructure() {
        return structure;
    }
}


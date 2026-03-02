package de.lmu.ifi.bio.gor.model;

public class Protein {
    private String id;
    private String sequence; // Die AS Zeile
    private String structure; // Die SS Zeile


    public String getId() {
        return id;
    }
    public void setId(String id) {
        this.id = id;
    }
    public String getSequence() {
        return sequence;
    }
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    public String getStructure() {
        return structure;
    }
    public void setStructure(String structure) {
        this.structure = structure;
    }
}


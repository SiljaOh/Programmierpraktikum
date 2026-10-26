package de.lmu.ifi.bio.gor.model;

import java.util.ArrayList;
import java.util.List;

public class Alignment {

    private final String id;
    private final String referenceSeq;
    private final String referenceSS;
    private final List<String> homologues = new ArrayList<>();

    public Alignment(String id, String referenceSeq, String referenceSS) {
        this.id = id;
        this.referenceSeq = referenceSeq;
        this.referenceSS = referenceSS;
    }

    public String getId() {
        return id;
    }
    public String getReferenceSeq() {
        return referenceSeq;
    }
    public String getReferenceSS() {
        return referenceSS;
    }
    public List<String> getHomologues() {
        return homologues;
    }
    public void addHomologues(String seq) {
        homologues.add(seq);
    }
}

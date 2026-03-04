package de.lmu.ifi.bio.gor.model;

import java.util.ArrayList;
import java.util.List;

public class Alignment {

    private String id; //> line
    private String referenceSeq; //AS line
    private String referenceSS; //SS line

    private final List<String> homologues = new ArrayList<String>();

    public Alignment(String id, String referenceSeq, String referenceSS) {
        this.id = id;
        this.referenceSeq = referenceSeq;
        this.referenceSS = referenceSS;
    }

    public String getId() {
        return id;
    }
    public void setId(String id) {
        this.id = id;
    }

    public String getReferenceSeq() {
        return referenceSeq;
    }
    public void setReferenceSeq(String referenceSeq) {
        this.referenceSeq = referenceSeq;
    }

    public String getReferenceSS() {
        return referenceSS;
    }
    public void setReferenceSS(String referenceSS) {
        this.referenceSS = referenceSS;
    }

    public List<String> getHomologues() {
        return homologues;
    }
    public void addHomologues(String seq) {
        this.homologues.add(seq);
    }

    //hier drin lassen?
    public String getWindow(String seq, int centerPos) {
        int m = 8;
        StringBuilder sb = new StringBuilder();
        for (int i = centerPos - m; i <= centerPos + m; i++) {
            if (i < 0 || i >= seq.length()) {
                sb.append("-"); // Padding
            } else {
                sb.append(seq.charAt(i));
            }
        }
        return sb.toString();
    }

}

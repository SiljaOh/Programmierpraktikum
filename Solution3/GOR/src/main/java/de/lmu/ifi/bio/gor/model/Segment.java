package de.lmu.ifi.bio.gor.model;

public class Segment {

    public char type;
    public int start;
    public int end;

    public Segment(char type, int start, int end) {
        this.type = type;   // 'this' sorgt dafür, dass das Feld der Klasse gefüllt wird
        this.start = start;
        this.end = end;
    }

    public int length() {
        return end - start + 1;
    }
}

package de.lmu.ifi.bio.gor.methods;

import de.lmu.ifi.bio.gor.model.Segment;

import java.util.ArrayList;
import java.util.List;

public class Evaluator {

    //Q3 Total Accuracy: Anteil der AS deren Struktur korrekt vorhergesagt wurde
    //Problem: ignoriert zusammenhängende Segmente
    //Vergleiche zb. HHHHHHHH und HHHHCHHH __> sehr hohe Q3 Obwohl Helix unterbrochen wird
    //Q3 = (Anzahl identischer Reste as Pos i) / (Gesamtzahl der Reste) * 100
    public double calculateQ3(String predicition, String reference) {

        if (predicition.isEmpty()) return 0.0;

        int correct = 0;
        for (int i = 0; i < reference.length(); i++) {
            if (predicition.charAt(i) == '-') continue;
            if (reference.charAt(i) == '-') continue;
            if (predicition.charAt(i) == reference.charAt(i)) correct++;
        }

        int total = (int) predicition.chars().filter(c -> c != '-').count();
        return (double) correct / total * 100.0;
    }

    //Qobs für qH , qE , qC (Q3 pro Zustand)
    //Wie oft kommt Zustand in SS vor (N) und wie oft in PS (n)
    // qH = n/N * 100
    public double calculateQobs(String prediction, String reference, char state) {

        int totalState = 0; //N
        int correctState = 0;//n
        for (int i = 0; i < reference.length(); i++) {
            if (prediction.charAt(i) == '-') continue;
            if (reference.charAt(i) == state) {
                totalState++;
                if (prediction.charAt(i) == state) {
                    correctState++;
                }
            }
        }
        if (totalState == 0) return -1.0; //Abfangen für - im Ouput
        return (double) correctState / totalState * 100.0;
    }

    //------------SOV Hilfsmethoden---------

    //teilt Strings dynamisch ein
    //"HHHHHCCCCCEEE" wird zu 3 Segment Objekten
    //Objekt 1: type: 'H', start: 0, end: 4 (Länge 5)
    //Objekt 2: type: 'C', start: 5, end: 9 (Länge 5)
    //Objekt 3: type: 'E', start: 10, end: 12 (Länge 3)
    public static List<Segment> getSegments(String s) {
        List<Segment> segments = new ArrayList<>();

        if (s == null || s.isEmpty()) return segments;

        int start = 0;
        char currentType = s.charAt(0);

        for (int i = 1; i < s.length(); i++) {
            if (s.charAt(i) != currentType) {
                segments.add(new Segment(currentType, start, i - 1));
                start = i;
                currentType = s.charAt(i);
            }
        }
        segments.add(new Segment(currentType, start, s.length() - 1));
        return segments;
    }

    //finde die Segmente die mit dem gewolltem char state übereinstimmen
    //für H wird aus [H,H,C,C,E,E] --> [H,H]
    private List<Segment> filterByType(List<Segment> list, char state) {
        List<Segment> filtered = new ArrayList<>();
        for (Segment s : list) {
            if (s.type == state) {
                filtered.add(s);
            }
        }
        return filtered;
    }

    //Überschneiden sich 2 Segment Objekte?
    private boolean isOverlapping(Segment s1, Segment s2) {
        return s1.start <= s2.end && s2.start <= s1.end;
    }


    //minov(s1,s2) = minimum overlap, Schnittmenge
    //Welche AS haben in beiden Segmenten den gleichen Strukturtyp an der gleichen Position
    private double getMinov(Segment s1, Segment s2) {
        return Math.min(s1.end, s2.end) - Math.max(s1.start, s2.start) + 1;
    }

    //maxov(s1,s2) = maximum overlap, Gesamtausdehnung
    //Berechne gesamte Länge, die beide Segmente abdecken
    private double getMaxov(Segment s1, Segment s2) {
        return Math.max(s1.end, s2.end) - Math.min(s1.start, s2.start) + 1;
    }

    //delta als realisitscher puffer bei minov um kleine ungenauigkiten nicht so hart
    // zu bestrafen wie zb. ein Segment das halbiert wurde (das muss härter bestraft werden)
    //kleine Verschiebungen werden dann mehr toleriert
    private double getDelta(Segment s1, Segment s2, double minov, double maxov) {
        // Delta = min[(maxov-minov), minov, len(s1)/2, len(s2)/2]
        double d1 = maxov - minov;
        double d2 = minov;
        double d3 = (double) s1.length() / 2;
        double d4 = (double) s2.length() / 2;
        return Math.min(Math.min(d1, d2), Math.min(d3, d4));
    }

    //SOV = Segment Overlap measure
    //basierend auf Strukturblöcken --> Segment Objekten!
    //basierend auf:
    //wo beginnen und enden Segmente
    //wie stark überlappen sie sich
    //wie groß ist die Abweichung bei den Rändern (Abweichung ist besser als Unterbrechungen!)
    //Formel: (100/N) * totalSum
    //mit N (Nominierung) bestraft wenn Segment nicht gefunden oder zerstückelt
    //und totalSum += ((minov + delta) / maxov) * s1.length()
    public double calculateSOV(String predicition, String reference) {

        StringBuilder sbPred = new StringBuilder();
        StringBuilder sbRef  = new StringBuilder();
        for (int i = 0; i < predicition.length(); i++) {
            if (predicition.charAt(i) != '-') {
                sbPred.append(predicition.charAt(i));
                sbRef.append(reference.charAt(i));
            }
        }
        String cleanPred = sbPred.toString();
        String cleanRef  = sbRef.toString();
        List<Segment> preds = getSegments(cleanPred);
        List<Segment> refs  = getSegments(cleanRef);

        double totalSum = 0;
        int N = 0;

        for (char state : new char[]{'H', 'E', 'C'}) {
            List<Segment> s1List = filterByType(refs, state);
            List<Segment> s2List = filterByType(preds, state);

            for (Segment s1 : s1List) {
                List<Segment> overlapsOfS1 = new ArrayList<>();
                for (Segment s2 : s2List) {
                    if (isOverlapping(s1, s2)) {
                        overlapsOfS1.add(s2);
                    }
                }

                if (overlapsOfS1.isEmpty()) {
                    //wenn segment verpasst
                    N += s1.length();
                } else {
                    //pro Überlappung die referenzlänge
                    N += s1.length() * overlapsOfS1.size();
                }

                for (Segment s2 : overlapsOfS1) {
                    double minov = getMinov(s1, s2);
                    double maxov = getMaxov(s1, s2);
                    double delta = getDelta(s1, s2, minov, maxov);
                    totalSum += ((minov + delta) / maxov) * s1.length();
                }
            }
        }

        if (N == 0) return -1.0;
        return (100.0 / N) * totalSum;
    }

    public double calculateStateSOV(String predicition, String reference, char state) {

        StringBuilder sbPred = new StringBuilder();
        StringBuilder sbRef  = new StringBuilder();
        for (int i = 0; i < predicition.length(); i++) {
            if (predicition.charAt(i) != '-') {
                sbPred.append(predicition.charAt(i));
                sbRef.append(reference.charAt(i));
            }
        }
        String cleanPred = sbPred.toString();
        String cleanRef  = sbRef.toString();
        List<Segment> preds = filterByType(getSegments(cleanPred), state);
        List<Segment> refs  = filterByType(getSegments(cleanRef), state);

        double stateSum = 0;
        int Nstate = 0;

        for (Segment s1 : refs) {
            List<Segment> overlapsOfS1 = new ArrayList<>();
            for (Segment s2 : preds) {
                if (isOverlapping(s1, s2)) {
                    overlapsOfS1.add(s2);
                }
            }

            if (overlapsOfS1.isEmpty()) {
                Nstate += s1.length();
            } else {
                Nstate += s1.length() * overlapsOfS1.size();
            }

            for (Segment s2 : overlapsOfS1) {
                double minov = getMinov(s1, s2);
                double maxov = getMaxov(s1, s2);
                double delta = getDelta(s1, s2, minov, maxov);
                stateSum += ((minov + delta) / maxov) * s1.length();
            }
        }

        if (Nstate == 0) return -1.0;
        return (100.0 / Nstate) * stateSum;
    }
}

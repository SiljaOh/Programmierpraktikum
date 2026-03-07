package alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Main {

    private static final Pattern CATH_ID =
            Pattern.compile("[0-9][A-Za-z0-9]{3}[A-Za-z0-9][0-9]{2}");

    private static String normalizeId(String id) {
        if (id == null) return null;
        String t = id.trim();

        Matcher m = CATH_ID.matcher(t);
        if (m.find()) return m.group().toLowerCase(Locale.ROOT);

        String cleaned = t.replaceAll("[^A-Za-z0-9]", "");
        m = CATH_ID.matcher(cleaned);
        if (m.find()) return m.group().toLowerCase(Locale.ROOT);

        return cleaned.isEmpty() ? null : cleaned.toLowerCase(Locale.ROOT);
    }

    private static Map<String, String> readSeqLib(String file) throws IOException {
        Map<String, String> seq = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                String[] parts = line.split(":");
                seq.put(parts[0], parts[1]);

            }
        }
        return seq;
    }

    private static void printAli(String id1Raw, String id2Raw, String a1, String a2, double score) {
        System.out.println(">" + id1Raw + " " + id2Raw + " " + score);
        System.out.println(id1Raw + ": " + a1);
        System.out.println(id2Raw + ": " + a2);
    }

    public static void main(String[] args) throws Exception {

        String mode = "global";
        String format = "scores";
        String pairsFile = null;
        String seqlibFile = null;
        String matrixFile = null;

        boolean nw = false;

        long go = 0;
        long ge = 0;

        for (int i = 0; i < args.length; i++) {
            String a = args[i];

            if (a.equals("--mode") && i + 1 < args.length) {
                mode = args[++i].trim();
            } else if (a.equals("--format") && i + 1 < args.length) {
                format = args[++i].trim();
            } else if (a.equals("--pairs") && i + 1 < args.length) {
                pairsFile = args[++i];
            } else if (a.equals("--seqlib") && i + 1 < args.length) {
                seqlibFile = args[++i];
            } else if (a.equals("-m") && i + 1 < args.length) {
                matrixFile = args[++i];
            } else if (a.equals("--nw")) {
                nw = true;
            } else if (a.equals("--go") && i + 1 < args.length) {
                go = Long.parseLong(args[++i]);
            } else if (a.equals("--ge") && i + 1 < args.length) {
                ge = Long.parseLong(args[++i]);
            } else if (a.equals("--dpmatrices") && i + 1 < args.length) {
                i++;
            } else if (a.equals("--check")) {
                // akzeptieren
            } else {
                // unbekannte Optionen ignorieren
            }
        }

        if (!nw || pairsFile == null || seqlibFile == null || matrixFile == null) {
            System.out.println(
                    "Usage: java -jar alignment.jar --nw --pairs <pairs> --seqlib <seqlib> -m <matrix> " +
                            "--go <v> --ge <v> --mode <global|local|freeshift> --format <scores|ali|html>"
            );
            return;
        }

        Map<String, String> sequences = readSeqLib(seqlibFile);

        SubstMatrix mat = SubstMatrix.load(Path.of(matrixFile));
        int scale = mat.scale;
        go *= scale;
        ge *= scale;

        GlobalAlignment ga = new GlobalAlignment(mat, go, ge);
        LocalAlignment la = new LocalAlignment(mat, go, ge);
        FreeShiftAlignment fa = new FreeShiftAlignment(mat, go, ge);

        try (BufferedReader br = new BufferedReader(new FileReader(pairsFile))) {
            String line;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                String[] p = line.split("\\s+");
                if (p.length < 2) continue;

                String id1Raw = p[0];
                String id2Raw = p[1];

                //String s1 = sequences.get(normalizeId(id1Raw));
                //String s2 = sequences.get(normalizeId(id2Raw));
                String s1 = sequences.get(id1Raw);
                String s2 = sequences.get(id2Raw);
                if (s1 == null || s2 == null) continue;

                long scoreScaled;
                String[] ali = null;

                if (mode.equals("global")) {
                    scoreScaled = ga.computeScore(s1, s2);
                    if (format.equals("ali")) ali = ga.alignment(s1, s2);

                } else if (mode.equals("local")) {
                    scoreScaled = la.computeScore(s1, s2);
                    if (format.equals("ali")) ali = la.alignment(s1, s2);

                } else if (mode.equals("freeshift")) {
                    scoreScaled = fa.computeScore(s1, s2);
                    if (format.equals("ali")) ali = fa.alignment(s1, s2);

                } else {
                    continue;
                }

                double score = scoreScaled / (double) scale;

                if (format.equals("scores")) {
                    System.out.printf(Locale.US, "%s %s %.4f%n", id1Raw, id2Raw, score);

                } else if (format.equals("ali")) {
                    if (ali == null || ali.length != 2) {
                        throw new IllegalStateException("alignment() must return String[2].");
                    }
                    printAli(id1Raw, id2Raw, ali[0], ali[1], score);

                } else if (format.equals("html")) {
                    System.out.printf(Locale.US, "<p>%s %s %.4f</p>%n", id1Raw, id2Raw, score);
                }
            }
        }
    }
}
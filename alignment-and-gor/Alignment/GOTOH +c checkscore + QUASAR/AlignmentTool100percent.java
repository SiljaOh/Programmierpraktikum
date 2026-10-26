package AlignmentTool;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/*
commandline Alignment Tool für Pairwise Alignments von Proteinsequenzen
Sequenzen aus einer seqlib lädt, Paare aus pairs-Datei lädt, eine Substitutionsmatrix lädt
Gotoh berechnet pro Paar ein Alignment(affine Gaps)
Ausgabe: Score v Score + Alignment
 */

public class AlignmentTool100percent {

    enum Mode { GLOBAL, LOCAL, FREESHIFT } //erlaubte Werte für -- Mode
    enum Format { SCORES, ALI } // erlaubte Werte für -- format


    private static final long NEG_INF = Long.MIN_VALUE / 4; // Use long everywhere to avoid overflow. /4 gibt puffer für zustände die "unmöglich" sind

    public static void main(String[] args) throws Exception {
        //Argumente parsen
        Args a = Args.parse(args);
        if (a == null) {
            System.out.println(Args.help());
            return;
        }
        //Daten einlesen, baut map, SubstMatrix.load : liest SUbstitutionsmatirx und baut Scoring-Objekt
        SeqLib lib = SeqLib.load(a.seqlibPath);
        SubstMatrix mat = SubstMatrix.load(a.matrixPath);
        int scale = mat.getScoreScale();
        long goScaled = a.go * scale; //damit Dezimalwerte auf long WErte
        long geScaled = a.ge * scale;

        List<Pair> pairs = readPairs(a.pairsPath);

        //über alle Paare iterieren, Alignment rechnen, Output bauen
        StringBuilder out = new StringBuilder();
        for (Pair p : pairs) {
            String s1 = lib.get(p.id1);
            String s2 = lib.get(p.id2);
            if (s1 == null || s2 == null) {
                throw new IllegalArgumentException("Pair references unknown id: " + p.id1 + " or " + p.id2);
            }

            GotohResult res = gotohAlign(s1, s2, mat, goScaled, geScaled, a.mode);
            //Ausgabe-Format scores
            if (a.format == Format.SCORES) {
                out.append(p.id1).append(" ").append(p.id2).append(" ")
                        .append(String.format(Locale.US, "%.4f", (double) res.score / scale))
                        .append("\n");
            } else {
                //Ausgabe-Format ali
                String fullA = res.fullAlignedA;
                String fullB = res.fullAlignedB;
                out.append(">").append(p.id1).append(" ").append(p.id2).append(" ")
                        .append(String.format(Locale.US, "%.4f", (double) res.score / scale))
                        .append("\n");
                out.append(p.id1).append(": ").append(fullA).append("\n");
                out.append(p.id2).append(": ").append(fullB).append("\n");

                // Kein #checkscore in Ausgabe – Server erwartet nur Header + 2 Sequenzzeilen
            }
        }

        System.out.print(out);
    }


    static class GotohResult {
        //zwei Versionen alignedA/B: nur der Teil den Backtracking tatsächlich gebaut hat
        //fullAlignedA/B: Prefix/Suffix dranbasteln damit Output die komplette Sequenzen enthält
        final long score;
        final String alignedA;
        final String alignedB;

        final String fullAlignedA;
        final String fullAlignedB;
        GotohResult(long score, String alignedA, String alignedB, String fullAlignedA, String fullAlignedB) {
            this.score = score;
            this.alignedA = alignedA;
            this.alignedB = alignedB;
            this.fullAlignedA = fullAlignedA;
            this.fullAlignedB = fullAlignedB;
        }
    }

    // Trace codes for backtracking
    // M: came from diag of (M/Ix/Iy), or from ZERO (local)
    private static final byte TM_FROM_M  = 1;
    private static final byte TM_FROM_IX = 2;
    private static final byte TM_FROM_IY = 3;
    private static final byte TM_FROM_0  = 4; // local reset

    // Ix: came from M (open) or Ix (extend) or ZERO (local)
    private static final byte TX_FROM_M  = 1;
    private static final byte TX_FROM_IX = 2;
    private static final byte TX_FROM_0  = 3;

    // Iy: came from M (open) or Iy (extend) or ZERO (local)
    private static final byte TY_FROM_M  = 1;
    private static final byte TY_FROM_IY = 2;
    private static final byte TY_FROM_0  = 3;

    // Zustands-IDs für Backtracking
    private static final int ST_M  = 0;
    private static final int ST_IX = 1;
    private static final int ST_IY = 2;


    static GotohResult gotohAlign(String A, String B, SubstMatrix mat, long go, long ge, Mode mode) {
        int n = A.length();
        int m = B.length();

        // 1-based chars (to simplify i-1)
        char[] a = (" " + A).toCharArray();
        char[] b = (" " + B).toCharArray();

        //Matrizen anlegen
        long[][] M  = new long[n + 1][m + 1];
        long[][] Ix = new long[n + 1][m + 1];
        long[][] Iy = new long[n + 1][m + 1];

        byte[][] tM  = new byte[n + 1][m + 1];
        byte[][] tIx = new byte[n + 1][m + 1];
        byte[][] tIy = new byte[n + 1][m + 1];

        // Initialisierung
        initMatrices(M, Ix, Iy, tM, tIx, tIy, n, m, go, ge, mode);

        // Fill-loop und bestScore Tracking
        long bestScore = (mode == Mode.LOCAL) ? 0 : NEG_INF;
        int bestI = 0, bestJ = 0, bestState = ST_M;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {

                // Ix (gap in B): (i-1,j) -> (i,j)
                long fromM_open = safeAdd(M[i - 1][j], go + ge);
                long fromX_ext  = safeAdd(Ix[i - 1][j], ge);

                long ixVal;
                byte ixTrace;
                if (fromM_open >= fromX_ext) {
                    ixVal = fromM_open;
                    ixTrace = TX_FROM_M;
                } else {
                    ixVal = fromX_ext;
                    ixTrace = TX_FROM_IX;
                }

                // Iy (gap in A): (i,j-1) -> (i,j)
                long fromM_openY = safeAdd(M[i][j - 1], go + ge);
                long fromY_ext   = safeAdd(Iy[i][j - 1], ge);

                long iyVal;
                byte iyTrace;
                if (fromM_openY >= fromY_ext) {
                    iyVal = fromM_openY;
                    iyTrace = TY_FROM_M;
                } else {
                    iyVal = fromY_ext;
                    iyTrace = TY_FROM_IY;
                }

                // M (match/mismatch): (i-1,j-1) -> (i,j)
                long diagM  = M[i - 1][j - 1];
                long diagIx = Ix[i - 1][j - 1];
                long diagIy = Iy[i - 1][j - 1];

                long bestDiag = diagM;
                byte mTrace = TM_FROM_M;
                if (diagIx > bestDiag) { bestDiag = diagIx; mTrace = TM_FROM_IX; }
                if (diagIy > bestDiag) { bestDiag = diagIy; mTrace = TM_FROM_IY; }

                long sub = mat.score(a[i], b[j]);
                long mVal = safeAdd(bestDiag, sub);

                // Local-Mode Anpassung(Reset auf 0)(wie smith-waterman Prinzip negative Scores auf 0)
                if (mode == Mode.LOCAL) {
                    // allow reset to 0
                    if (ixVal < 0) { ixVal = 0; ixTrace = TX_FROM_0; }
                    if (iyVal < 0) { iyVal = 0; iyTrace = TY_FROM_0; }
                    if (mVal  < 0) { mVal  = 0; mTrace  = TM_FROM_0; }
                }
                //Werte speichern
                Ix[i][j] = ixVal; tIx[i][j] = ixTrace;
                Iy[i][j] = iyVal; tIy[i][j] = iyTrace;
                M[i][j]  = mVal;  tM[i][j]  = mTrace;

                // Track best for local, or for freeshift we handle after fill.
                if (mode == Mode.LOCAL) {
                    // best anywhere across all 3
                    if (M[i][j] > bestScore) { bestScore = M[i][j]; bestI = i; bestJ = j; bestState = ST_M; }
                    if (Ix[i][j] > bestScore){ bestScore = Ix[i][j]; bestI = i; bestJ = j; bestState = ST_IX; }
                    if (Iy[i][j] > bestScore){ bestScore = Iy[i][j]; bestI = i; bestJ = j; bestState = ST_IY; }
                }
            }
        }

        // GLOBAL nach dem fill wählen wo ich starte
        int startI, startJ, startState;
        long finalScore;

        if (mode == Mode.GLOBAL) {
            startI = n; startJ = m;
            long sM = M[n][m], sX = Ix[n][m], sY = Iy[n][m];
            if (sM >= sX && sM >= sY) { finalScore = sM; startState = ST_M; }
            else if (sX >= sY)        { finalScore = sX; startState = ST_IX; }
            else                      { finalScore = sY; startState = ST_IY; }
        } else if (mode == Mode.LOCAL) {
            startI = bestI; startJ = bestJ; startState = bestState; finalScore = bestScore;
        } else { // FREESHIFT (semi-global): best on last row or last column (across states)
            long best = NEG_INF;
            int bi = n, bj = m, bs = ST_M;

            // last row: i=n, j=0..m
            for (int j = 0; j <= m; j++) {
                if (M[n][j]  > best) { best = M[n][j];  bi = n; bj = j; bs = ST_M; }
                if (Ix[n][j] > best) { best = Ix[n][j]; bi = n; bj = j; bs = ST_IX; }
                if (Iy[n][j] > best) { best = Iy[n][j]; bi = n; bj = j; bs = ST_IY; }
            }
            // last col: j=m, i=0..n
            for (int i = 0; i <= n; i++) {
                if (M[i][m]  > best) { best = M[i][m];  bi = i; bj = m; bs = ST_M; }
                if (Ix[i][m] > best) { best = Ix[i][m]; bi = i; bj = m; bs = ST_IX; }
                if (Iy[i][m] > best) { best = Iy[i][m]; bi = i; bj = m; bs = ST_IY; }
            }

            startI = bi; startJ = bj; startState = bs; finalScore = best;
        }

        // Backtrack
        StringBuilder alA = new StringBuilder();
        StringBuilder alB = new StringBuilder();
        int i = startI, j = startJ, st = startState;

        while (i > 0 || j > 0) {
            if (mode == Mode.LOCAL) {
                // stop at reset
                if (st == ST_M && tM[i][j] == TM_FROM_0) break;
                if (st == ST_IX && tIx[i][j] == TX_FROM_0) break;
                if (st == ST_IY && tIy[i][j] == TY_FROM_0) break;
            }
            // a[0] und b[0] sind Leerzeichen (1-based Arrays) – nie anhängen
            if (st == ST_M && (i <= 0 || j <= 0)) break;
            if (st == ST_IX && i <= 0) break;
            if (st == ST_IY && j <= 0) break;

            if (st == ST_M) {
                byte tr = tM[i][j];
                // diag move
                alA.append(a[i]);
                alB.append(b[j]);
                i--; j--;
                if (tr == TM_FROM_M) st = ST_M;
                else if (tr == TM_FROM_IX) st = ST_IX;
                else if (tr == TM_FROM_IY) st = ST_IY;
                else break; // local reset
            } else if (st == ST_IX) {
                byte tr = tIx[i][j];
                // up move: A[i] vs '-'
                alA.append(a[i]);
                alB.append('-');
                i--;
                if (tr == TX_FROM_M) st = ST_M;
                else if (tr == TX_FROM_IX) st = ST_IX;
                else break; // local reset
            } else { // ST_IY
                byte tr = tIy[i][j];
                // left move: '-' vs B[j]
                alA.append('-');
                alB.append(b[j]);
                j--;
                if (tr == TY_FROM_M) st = ST_M;
                else if (tr == TY_FROM_IY) st = ST_IY;
                else break; // local reset
            }

            // Safety (shouldn't happen, but avoids infinite loop if traces are broken)
            if (i < 0 || j < 0) break;
        }

        String alignedA = alA.reverse().toString();
        String alignedB = alB.reverse().toString();

        // Full alignment bauen (Prefix/Suffix bei local/freeshift)
        String fullA = alignedA;
        String fullB = alignedB;
        if ((mode == Mode.LOCAL || mode == Mode.FREESHIFT) && (i > 0 || j > 0 || startI < n || startJ < m)) {
            // i,j can be 0 or negative when loop breaks; clamp to avoid String.repeat(-1) and substring issues
            int endI = Math.max(0, i);
            int endJ = Math.max(0, j);
            int si = Math.min(startI, n);
            int sj = Math.min(startJ, m);
            int suffixGapsA = Math.max(0, m - sj);
            int suffixGapsB = Math.max(0, n - si);
            String prefixA = A.substring(0, endI) + "-".repeat(endJ);
            String prefixB = "-".repeat(endI) + B.substring(0, endJ);
            String suffixA = A.substring(si, n) + "-".repeat(suffixGapsA);
            String suffixB = "-".repeat(suffixGapsB) + B.substring(sj, m);
            fullA = prefixA + alignedA + suffixA;
            fullB = prefixB + alignedB + suffixB;
        }

        return new GotohResult(finalScore, alignedA, alignedB, fullA, fullB);
    }

    private static void initMatrices(
            long[][] M, long[][] Ix, long[][] Iy,
            byte[][] tM, byte[][] tIx, byte[][] tIy,
            int n, int m, long go, long ge, Mode mode
    ) {
        // Fill with NEG_INF by default
        for (int i = 0; i <= n; i++) {
            Arrays.fill(M[i], NEG_INF);
            Arrays.fill(Ix[i], NEG_INF);
            Arrays.fill(Iy[i], NEG_INF);
        }

        if (mode == Mode.GLOBAL) {
            M[0][0] = 0;
            Ix[0][0] = Iy[0][0] = NEG_INF;

            // first column: align A[1..i] with gaps in B => Ix valid
            for (int i = 1; i <= n; i++) {
                Ix[i][0] = go + (long) i * ge;
                tIx[i][0] = (i == 1) ? TX_FROM_M : TX_FROM_IX; // not used much at borders
                // M[i][0], Iy[i][0] stay NEG_INF
            }
            // first row: align B[1..j] with gaps in A => Iy valid
            for (int j = 1; j <= m; j++) {
                Iy[0][j] = go + (long) j * ge;
                tIy[0][j] = (j == 1) ? TY_FROM_M : TY_FROM_IY;
                // M[0][j], Ix[0][j] stay NEG_INF
            }
        }
        else if (mode == Mode.LOCAL) {
            // local: allow start anywhere -> 0 baseline
            for (int i = 0; i <= n; i++) {
                M[i][0] = 0; Ix[i][0] = 0; Iy[i][0] = 0;
                tM[i][0] = TM_FROM_0; tIx[i][0] = TX_FROM_0; tIy[i][0] = TY_FROM_0;
            }
            for (int j = 0; j <= m; j++) {
                M[0][j] = 0; Ix[0][j] = 0; Iy[0][j] = 0;
                tM[0][j] = TM_FROM_0; tIx[0][j] = TX_FROM_0; tIy[0][j] = TY_FROM_0;
            }
        }
        else { // FREESHIFT (semi-global): free end gaps and free start shifts (common variant)
            // Make starts along borders free (0) so you can "shift" without penalty.
            for (int i = 0; i <= n; i++) {
                M[i][0] = 0;
                Ix[i][0] = 0;
                Iy[i][0] = NEG_INF;
                tM[i][0] = TM_FROM_0;
                tIx[i][0] = TX_FROM_0;
            }
            for (int j = 0; j <= m; j++) {
                M[0][j] = 0;
                Iy[0][j] = 0;
                Ix[0][j] = NEG_INF;
                tM[0][j] = TM_FROM_0;
                tIy[0][j] = TY_FROM_0;
            }
            // Note: This variant is often "good enough" for course freeshift.
            // If your reference scores differ, adjust init + where you take the max at the end.
        }
    }

        // NEG_INF + irgendwas = NEG_INF
    private static long safeAdd(long x, long y) {
        if (x <= NEG_INF / 2) return NEG_INF;
        // We don't expect overflow with realistic scores, but keep it safe-ish:
        long r = x + y;
        // If overflow happens, clamp (very unlikely here)
        if (((x ^ r) & (y ^ r)) < 0) {
            return (r < 0) ? Long.MIN_VALUE / 2 : Long.MAX_VALUE / 2;
        }
        return r;
    }


    //score eines gegebenen Alignments nachrechnen, unabhängig von DP prüfen: ist Score consistent
    static long checkscore(String alignedA, String alignedB, SubstMatrix mat, long go, long ge) {
        if (alignedA.length() != alignedB.length()) {
            throw new IllegalArgumentException("Alignment strings must have equal length.");
        }
        long score = 0;

        boolean inGapA = false; // currently in a run where A has '-' (gap in A)
        boolean inGapB = false; // currently in a run where B has '-' (gap in B)

        for (int k = 0; k < alignedA.length(); k++) {
            char ca = alignedA.charAt(k);
            char cb = alignedB.charAt(k);

            if (ca == '-' && cb == '-') {
                throw new IllegalArgumentException("Invalid alignment: both positions are gaps at k=" + k);
            }

            if (ca == '-') {
                // gap in A (align '-' with B char)
                if (!inGapA) score += (go + ge);
                else score += ge;
                inGapA = true;
                inGapB = false;
            } else if (cb == '-') {
                // gap in B (align A char with '-')
                if (!inGapB) score += (go + ge);
                else score += ge;
                inGapB = true;
                inGapA = false;
            } else {
                // match/mismatch
                score += mat.score(ca, cb);
                inGapA = false;
                inGapB = false;
            }
        }

        return score;
    }


    static class Pair {
        final String id1, id2;
        Pair(String id1, String id2) { this.id1 = id1; this.id2 = id2; }
    }

    static List<Pair> readPairs(Path p) throws IOException {
        List<String> lines = Files.readAllLines(p, StandardCharsets.UTF_8);
        List<Pair> res = new ArrayList<>();
        for (String line : lines) {
            line = line.trim();
            if (line.isEmpty() || line.startsWith("#")) continue;
            String[] tok = line.split("\\s+");
            if (tok.length < 2) continue;
            res.add(new Pair(tok[0], tok[1]));
        }
        return res;
    }

    static class SeqLib {
        private final Map<String, String> idToSeq = new HashMap<>();

        static SeqLib load(Path p) throws IOException {
            SeqLib lib = new SeqLib();
            List<String> lines = Files.readAllLines(p, StandardCharsets.UTF_8);
            for (String line : lines) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                int idx = line.indexOf(':');
                if (idx < 1) continue;
                String id = line.substring(0, idx).trim();
                String seq = line.substring(idx + 1).trim();
                lib.idToSeq.put(id, seq);
            }
            return lib;
        }

        String get(String id) { return idToSeq.get(id); }
    }

    static class SubstMatrix {
        private final Map<Character, Integer> index = new HashMap<>();
        private final long[][] score;
        private final int scoreScale;

        private SubstMatrix(long[][] score, Map<Character, Integer> index, int scoreScale) {
            this.score = score;
            this.index.putAll(index);
            this.scoreScale = scoreScale;
        }

        int getScoreScale() { return scoreScale; }

        long score(char a, char b) {
            Integer ia = index.get(Character.toUpperCase(a));
            Integer ib = index.get(Character.toUpperCase(b));
            if (ia == null || ib == null) {
                throw new IllegalArgumentException("Unknown residue in scoring matrix: '" + a + "' or '" + b + "'");
            }
            return score[ia][ib];
        }


        static SubstMatrix load(Path p) throws IOException {
            List<String> raw = Files.readAllLines(p, StandardCharsets.UTF_8);

            String rowIndex = null;
            String colIndex = null;

            // We will read the lower-triangular matrix rows.
            List<double[]> triRows = new ArrayList<>();

            for (String line : raw) {
                String t = line.trim();
                if (t.isEmpty()) continue;
                if (t.startsWith("#") || t.startsWith(";") || t.startsWith("//")) continue;

                if (t.startsWith("ROWINDEX")) {
                    // e.g. ROWINDEX ARNDCQEGHILKMFPSTWYV
                    String[] parts = t.split("\\s+");
                    if (parts.length >= 2) rowIndex = parts[1].trim();
                    continue;
                }
                if (t.startsWith("COLINDEX")) {
                    String[] parts = t.split("\\s+");
                    if (parts.length >= 2) colIndex = parts[1].trim();
                    continue;
                }

                if (t.startsWith("MATRIX")) {
                    // e.g. MATRIX     -2.      8.   (triangular numeric row)
                    String rest = t.substring("MATRIX".length()).trim();
                    if (rest.isEmpty()) continue;

                    String[] nums = rest.split("\\s+");
                    double[] row = new double[nums.length];
                    for (int i = 0; i < nums.length; i++) {
                        // values look like "6." or "-2." -> Double.parseDouble works
                        row[i] = Double.parseDouble(nums[i]);
                    }
                    triRows.add(row);
                }
            }

            if (rowIndex == null || colIndex == null) {
                throw new IllegalArgumentException("Matrix file missing ROWINDEX/COLINDEX: " + p);
            }
            if (rowIndex.length() != colIndex.length()) {
                throw new IllegalArgumentException("ROWINDEX and COLINDEX have different lengths: " + rowIndex.length() + " vs " + colIndex.length());
            }

            int n = rowIndex.length();
            if (triRows.size() < n) {
                throw new IllegalArgumentException("Matrix parse incomplete: expected at least " + n + " MATRIX rows, got " + triRows.size());
            }

            // Build index map
            Map<Character, Integer> idx = new HashMap<>();
            for (int i = 0; i < n; i++) {
                idx.put(Character.toUpperCase(rowIndex.charAt(i)), i);
            }

            // Dezimal-Matrizen (z. B. Dayhoff) in Zehntel-Einheiten speichern
            boolean hasDecimals = false;
            for (int r = 0; r < n; r++) {
                double[] row = triRows.get(r);
                for (int c = 0; c <= r && c < row.length; c++) {
                    if (Math.abs(row[c] - Math.round(row[c])) > 1e-9) {
                        hasDecimals = true;
                        break;
                    }
                }
                if (hasDecimals) break;
            }
            int scoreScale = hasDecimals ? 10 : 1;

            // Fill full symmetric matrix from lower triangle
            long[][] sc = new long[n][n];
            for (int r = 0; r < n; r++) {
                double[] row = triRows.get(r);
                if (row.length < r + 1) {
                    throw new IllegalArgumentException("MATRIX row " + r + " too short: expected " + (r + 1) + " numbers, got " + row.length);
                }
                for (int c = 0; c <= r; c++) {
                    long val = scoreScale == 10 ? Math.round(row[c] * 10) : Math.round(row[c]);
                    sc[r][c] = val;
                    sc[c][r] = val;
                }
            }

            return new SubstMatrix(sc, idx, scoreScale);
        }
    }



    static class Args {
        Path pairsPath;
        Path seqlibPath;
        Path matrixPath;

        long go = -12;
        long ge = -1;

        Mode mode = Mode.GLOBAL;
        Format format = Format.ALI;

        static Args parse(String[] args) {
            Args a = new Args();
            for (int i = 0; i < args.length; i++) {
                String s = args[i];

                switch (s) {
                    case "--pairs":
                        a.pairsPath = Paths.get(requireNext(args, ++i, "--pairs"));
                        break;
                    case "--seqlib":
                        a.seqlibPath = Paths.get(requireNext(args, ++i, "--seqlib"));
                        break;
                    case "-m":
                        a.matrixPath = Paths.get(requireNext(args, ++i, "-m"));
                        break;
                    case "--go":
                        a.go = Long.parseLong(requireNext(args, ++i, "--go"));
                        break;
                    case "--ge":
                        a.ge = Long.parseLong(requireNext(args, ++i, "--ge"));
                        break;
                    case "--mode":
                        String md = requireNext(args, ++i, "--mode").toLowerCase(Locale.ROOT);
                        if (md.equals("global")) a.mode = Mode.GLOBAL;
                        else if (md.equals("local")) a.mode = Mode.LOCAL;
                        else if (md.equals("freeshift")) a.mode = Mode.FREESHIFT;
                        else return null;
                        break;
                    case "--format":
                        String f = requireNext(args, ++i, "--format").toLowerCase(Locale.ROOT);
                        if (f.equals("scores")) a.format = Format.SCORES;
                        else if (f.equals("ali")) a.format = Format.ALI;
                        else return null;
                        break;
                    case "-h":
                    case "--help":
                        return null;
                    default:
                        // ignore unknown flags in this minimal tester OR fail:
                        // return null;
                        break;
                }
            }

            if (a.pairsPath == null || a.seqlibPath == null || a.matrixPath == null) return null;
            return a;
        }

        static String requireNext(String[] args, int i, String flag) {
            if (i >= args.length) throw new IllegalArgumentException("Missing value for " + flag);
            return args[i];
        }

        static String help() {
            return ""
                    + "Usage:\n"
                    + "  java AlignmentTool --pairs <pairs> --seqlib <seqlib> -m <matrix> [--go -12] [--ge -1]\n"
                    + "                   [--mode global|local|freeshift] [--format scores|ali]\n"
                    + "\n"
                    + "Input formats:\n"
                    + "  seqlib:  id:SEQUENCE   (one per line)\n"
                    + "  pairs:   id1 id2 [ignored...]\n"
                    + "\n"
                    + "Gap model:\n"
                    + "  Gap(k) = go + k*ge  (first gap char costs go+ge)\n";
        }
    }
}


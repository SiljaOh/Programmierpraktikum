package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.methods.*;
import de.lmu.ifi.bio.gor.model.Alignment;
import de.lmu.ifi.bio.gor.model.Gor3Matrix;
import de.lmu.ifi.bio.gor.model.Gor4Matrix;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.io.AlnParser;


import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.List;


public class Predict {


    public static void main(String[] args) throws IOException {

        String modelPath = null;
        String seqPath = null;
        //Default
        String format = "txt";
        boolean probabilities = false;
        boolean showScores = false;
        String alnPath = null;
        //für parameter tuning:
        double pseudoCount = 0.8;
        int minH = 4;
        int minE = 2;

        //jedes Argument durch gehen
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--model":
                    if (i + 1 < args.length) modelPath = args[++i];
                    break;
                case "--seq":
                    if (i + 1 < args.length) seqPath = args[++i];
                    break;
                case "--format":
                    if (i + 1 < args.length) format = args[++i];
                    break;
                case "--probabilities":
                    probabilities = true;
                    break;
                case "--show-scores":
                    showScores = true;
                    break;
                case "--maf":
                    if (i + 1 < args.length) alnPath = args[++i];
                    break;
                case "--minH":
                    if (i + 1 < args.length) minH = Integer.parseInt(args[++i]);
                    break;
                case "--minE":
                    if (i + 1 < args.length) minE = Integer.parseInt(args[++i]);
                    break;
                case "--pseudo":
                    if (i + 1 < args.length) pseudoCount = Double.parseDouble(args[++i]);
                    break;
                //wenn nichts angegeben
                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.err.println("Usage: java -jar predict.jar [--probabilities]\n" +
                            "--model <model file> --format <txt|html>\n" +
                            "{--seq <fasta file>|--maf <multiple-alignment-folder>}");
                    System.exit(1);
            }
        }
        //Pflichtangaben
        if (modelPath == null || (seqPath == null && alnPath == null)) {
            System.err.println("Usage: java -jar predict.jar [--probabilities] [--show-scores] \n" +
                    "--model <model file> --format <txt|html>\n" +
                    "{--seq <fasta file>|--maf <multiple-alignment-folder>}");
            System.exit(1);
        }

        String methodType = GorMethod.whichMethod(modelPath);
        GorMethod gor;

        switch (methodType) {
            case "GOR3":
                gor = new Gor3();
                break;
            case "GOR4":
                gor = new Gor4();
                break;
            default:
                gor = new Gor1();
        }

        gor.getMatrix().loadModel(modelPath);

        gor.setPseudo(pseudoCount);

        if (gor instanceof Gor3) {
            ((Gor3Matrix) gor.getMatrix()).calculateScores3();
        } else if (gor instanceof Gor4) {
            ((Gor4Matrix) gor.getMatrix()).calculateScores4();
        }

        //bonus parameter tuning
        gor.setMinH(minH);
        gor.setMinE(minE);
        //Bonus probs
        gor.setUseProbability(probabilities);
        gor.getMatrix().setUseProbability(probabilities);
        gor.setShowScores(showScores);


        //für html
        if ("html".equalsIgnoreCase(format)) {
            System.out.println("<div class='gor-results' style='font-family:monospace; line-height:1.2;'>");
        }

        //entweder --seq oder --maf
        if (seqPath != null) {
            List<Protein> proteins = GorMethod.readFasta(seqPath);
            for (Protein p : proteins) {
                if (p.getSequence().length() < 17) {
                    continue;
                }
                String result = gor.predict(p.getSequence());

                //BONUS postprocessing
                String processedresult = gor.postProcess(result);
                //JMOL etc auf website iengebaut!!
                printResult(p.getId(), p.getSequence(), processedresult, probabilities, showScores, gor, format);


            }
        } else {
            File mafDir = new File(alnPath);
            File[] alnFiles = mafDir.listFiles(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".aln");
                }
            });
            //jetzt alle files die da drin
            if (alnFiles != null) {
                for (File alnFile : alnFiles) {
                    //jede file inein pbjekt mit parser
                    Alignment alignment = AlnParser.parseAlignment(alnFile.getPath());
                    if (alignment == null) continue;
                    String result = gor.predictGor5(alignment);

                    //BONUS postprocessing
                    String processedresult = gor.postProcess(result);

                    printResult(alignment.getId(), alignment.getReferenceSeq(), processedresult, probabilities, showScores, gor, format);
                }
            }
        }
        if ("html".equalsIgnoreCase(format)) {
            System.out.println("</div>");
        }
    }

    //output funktion
    private static void printResult(String id, String as, String ps, boolean useProbabilities, boolean showScores, GorMethod gor, String format) {
        if ("html".equalsIgnoreCase(format)) {
            System.out.println("<div style='margin-bottom:20px; border-left: 4px solid #0d6b5c; padding-left:10px;'>");
            System.out.println("  <div>&gt; " + id + "</div>");
            System.out.println("  <div>AS " + as + "</div>");

            // PS mit Farben für die Web-Vorschau
            System.out.print("  <div>PS ");
            for (char c : ps.toCharArray()) {
                String color = (c == 'H') ? "#e63946" : (c == 'E' ? "#457b9d" : "#2a9d8f");
                System.out.print("<span style='color:" + color + "; font-weight:bold;'>" + c + "</span>");
            }
            System.out.println("</div>");

            if (useProbabilities || showScores) {
                printProbabilitiesHtml(as, gor, useProbabilities);
            }
            System.out.println("</div>");

        } else {
            // txt Output
            System.out.println("> " + id);
            System.out.println("AS " + as);
            System.out.println("PS " + ps);
            if (useProbabilities || showScores) {
                printExtra(as, gor, useProbabilities);
            }
            System.out.println();
        }
    }

    private static void printExtra(String sequence, GorMethod gor, boolean useProbability) {
        double[][] raw = gor.getScores(sequence);
        String[] labels;
        if (useProbability) {
            labels = new String[]{"PH", "PE", "PC"};
        } else {
            labels = new String[]{"IH", "IE", "IC"};
        }

        int m = 8;
        //jeder zustand
        for (int s = 0; s < 3; s++) {
            StringBuilder sb = new StringBuilder(labels[s] + " --------");
            //ohne ränder
            for (int i = m; i < raw.length - m; i++) {
                int digit;
                if (useProbability) {
                    //sigmoid
                    double p = 1.0 / (1.0 + Math.exp(-raw[i][s]));
                    //richitg runden
                    digit = (int) Math.round(p * 10);
                } else {
                    double[] exp = {Math.exp(raw[i][0]), Math.exp(raw[i][1]), Math.exp(raw[i][2])};
                    double sum = exp[0] + exp[1] + exp[2];
                    //softmax
                    digit = (int) Math.round(exp[s] / sum * 10);
                }
                //keine 10 ausgegeben falss wäre aufgerundet 10
                sb.append(Math.min(9, digit));
            }
            sb.append("--------");
            System.out.println(sb);
        }
    }

    private static void printProbabilitiesHtml(String sequence, GorMethod gor, boolean useProbability) {
        double[][] raw = gor.getScores(sequence);
        String[] labels;
        if (useProbability) {
            labels = new String[]{"PH", "PE", "PC"};
        } else {
            labels = new String[]{"IH", "IE", "IC"};
        }
        //H ist rot, E ist blau, C ist grün
        String[] colors = {"#e63946", "#457b9d", "#2a9d8f"};
        int m = 8;

        for (int s = 0; s < 3; s++) {
            System.out.print("  <div>" + labels[s] + " --------");
            for (int i = m; i < raw.length - m; i++) {
                int val;
                if (useProbability) {
                    // wie oben sigmoid
                    double p = 1.0 / (1.0 + Math.exp(-raw[i][s]));
                    val = Math.min(9, (int) Math.round(p * 10));
                } else {
                    double[] exp = {Math.exp(raw[i][0]), Math.exp(raw[i][1]), Math.exp(raw[i][2])};
                    double sum = exp[0] + exp[1] + exp[2];
                    val = Math.min(9, (int) Math.round(exp[s] / sum * 10));
                }
                double opacity = 0.15 + (val / 9.0) * 0.85;
                System.out.print(
                        "<span style='color:" + colors[s] + "; font-weight:bold; opacity:" +
                                String.format("%.2f", opacity) + ";'>" + val + "</span>"
                );
            }
            System.out.println("--------</div>");
        }
    }
}



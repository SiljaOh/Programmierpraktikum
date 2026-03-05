package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.methods.Gor1;
import de.lmu.ifi.bio.gor.methods.Gor3;
import de.lmu.ifi.bio.gor.methods.Gor4;
import de.lmu.ifi.bio.gor.methods.GorMethod;
import de.lmu.ifi.bio.gor.model.Alignment;
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
        String alnPath = null;

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
                case "--maf":
                    if (i + 1 < args.length) alnPath = args[++i];
                    break;
                //wenn nichts angegeben
                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.err.println("Usage: java -jar predict.jar --model <file> --seq <file> [--format <txt|html>]");
                    System.err.println("Or:    java -jar predict.jar --model <file> --maf <folder> [--format <txt|html>]");
                    System.exit(1);
            }
        }
        //Pflichtangaben
        if (modelPath == null || (seqPath == null && alnPath == null)) {
            System.err.println("Usage: java -jar predict.jar --model <file> --seq <file> [--format <txt|html>]");
            System.err.println("Or:    java -jar predict.jar --model <file> --maf <folder> [--format <txt|html>]");
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

        //entweder --seq oder --maf
        if (seqPath != null) {
            List<Protein> proteins = GorMethod.readFasta(seqPath);
            for (Protein p : proteins) {
                if (p.getSequence().length() < 17) continue;
                String result = gor.predict(p.getSequence());

                //BONUS postprocessing
                String processedresult = gor.postProcess(result);

                printResult(p.getId(), p.getSequence(), processedresult, probabilities, gor);
            }
        } else {
            File mafDir = new File(alnPath);
            File[] alnFiles = mafDir.listFiles(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".aln");
                }
            });
            //jetzt alle fies da drin
            if (alnFiles != null) {
                for (File alnFile : alnFiles) {
                    //jede file inein pbjekt mit parser
                    Alignment alignment = AlnParser.parseAlignment(alnFile.getPath());
                    if (alignment == null) continue;
                    String result = gor.predictGor5(alignment);

                    //BONUS postprocessing
                    String processedresult = gor.postProcess(result);

                    printResult(alignment.getId(), alignment.getReferenceSeq(), processedresult, probabilities, gor);
                }
            }
        }
    }

    //output funktion
    private static void printResult(String id, String as, String ps, boolean probabilities, GorMethod gor) {
        System.out.println("> " + id);
        System.out.println("AS " + as);
        System.out.println("PS " + ps);
        if (probabilities) {
            printProbabilities(as, gor);
            System.out.println();
        }
    }

    private static void printProbabilities(String sequence, GorMethod gor) {
        double[][] raw = gor.getScores(sequence);
        String[] labels = {"PH", "PE", "PC"};
        int m = 8;
        //jeder zustand
        for (int s = 0; s < 3; s++) {
            StringBuilder sb = new StringBuilder(labels[s] + " --------");
            //ohne ränder
            for (int i = m; i < raw.length - m; i++) {
                //softmax auf alle 3 scores -> 0-1 wslk aus log odds wert, alle pos
                double[] exp = {Math.exp(raw[i][0]), Math.exp(raw[i][1]), Math.exp(raw[i][2])};
                //alle werte
                double sum = exp[0] + exp[1] + exp[2];
                // P(S|a) = exp(score_s) / summe, mal 10 und auf ganze zahl runden → 0-9
                sb.append((int)(exp[s] / sum * 10));
            }
            sb.append("--------");
            System.out.println(sb);
        }
    }
}



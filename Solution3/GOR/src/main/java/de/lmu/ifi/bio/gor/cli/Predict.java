package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.methods.Gor1;
import de.lmu.ifi.bio.gor.methods.Gor3;
import de.lmu.ifi.bio.gor.methods.Gor4;
import de.lmu.ifi.bio.gor.methods.GorMethod;
import de.lmu.ifi.bio.gor.model.Alignment;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.io.AlnParser;

import java.io.File;
import java.io.IOException;
import java.util.List;

import static de.lmu.ifi.bio.gor.io.AlnParser.parseAlignment;

public class Predict {

    public static void main(String[] args) throws IOException {

        //java -jar predict.jar --model <file> --seq <file> [--format <txt|html>] [--probabilities]

        String modelPath = null;
        String seqPath = null;
        //Default
        String format = "txt";
        boolean probabilities = false;
        String alnPath = null;

        for (int i = 0; i < args.length; i++) {
            //jedes Argument durch gehen
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
                    //Vorhanden oder nicht vorhanden
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

        // ... nach der Argumente-Verarbeitung ...

        String methodType = GorMethod.whichMethod(modelPath);
        GorMethod gor;

        if (methodType.equals("GOR1")) {
            gor = new Gor1();
            gor.getMatrix().loadModel(modelPath);
        } else if (methodType.equals("GOR3")) {
            gor = new Gor3();
            ((Gor3) gor).getMatrix3().loadModel3(modelPath);
        } else {
            gor = new Gor4();
            ((Gor4) gor).getMatrix4().loadModel4(modelPath);
        }

        //entweder --seq oder --maf
        if (seqPath != null) {
            // Standard GOR Vorhersage
            List<Protein> proteins = GorMethod.readFasta(seqPath);
            for (Protein p : proteins) {
                if (p.getSequence().length() < 17) continue;
                String result = gor.predict(p.getSequence());
                printResult(p.getId(), p.getSequence(), result);
            }
        } else if (alnPath != null) {
            // GOR V
            File mafDir = new File(alnPath);
            File[] alnFiles = mafDir.listFiles((dir, name) -> name.endsWith(".aln"));

            if (alnFiles != null) {
                for (File alnFile : alnFiles) {
                    Alignment alignment = AlnParser.parseAlignment(alnFile.getPath());
                    if (alignment == null) continue;

                    String result = gor.predictGor5(alignment);
                    printResult(alignment.getId(), alignment.getReferenceSeq(), result);
                }
            }
        }
    }

    //output funktion
    private static void printResult(String id, String as, String ps) {
        System.out.println("> " + id);
        System.out.println("AS " + as);
        System.out.println("PS " + ps + "\n");
    }

}



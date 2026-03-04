package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.methods.Gor1;
import de.lmu.ifi.bio.gor.methods.Gor3;
import de.lmu.ifi.bio.gor.methods.Gor4;
import de.lmu.ifi.bio.gor.methods.GorMethod;
import de.lmu.ifi.bio.gor.model.Protein;

import java.io.IOException;
import java.util.List;

public class Predict {

    public static void main(String[] args) throws IOException {

        //java -jar predict.jar --model <file> --seq <file> [--format <txt|html>] [--probabilities]

        String modelPath = null;
        String seqPath = null;
        //Default
        String format = "txt";
        boolean probabilities = false;

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
                //nicht sangegeben
                default:
                    System.exit(1);
            }
        }
        // Pflichtargumente prüfen
        if (modelPath == null || seqPath == null) {
            System.exit(1);
        }



        String GOR = GorMethod.whichMethod(modelPath);

        if (GOR.equals("GOR1")) {

            Gor1 gor1 = new Gor1();
            gor1.getMatrix().loadModel(modelPath);

            List<Protein> proteins = GorMethod.readFasta(seqPath);

            for (Protein protein : proteins) {
                //um auf selbe Anzahl Protein zu kommen
                if (protein.getSequence().length() < 17) continue;
                String result = gor1.predict(protein.getSequence());
                System.out.println("> " + protein.getId());
                System.out.println("AS " + protein.getSequence());
                System.out.println("PS " + result + "\n");
            }

        } else if (GOR.equals("GOR3")) {
            Gor3 gor3 = new Gor3();
            gor3.getMatrix3().loadModel3(modelPath);
            List<Protein> proteins = GorMethod.readFasta(seqPath);
            for (Protein protein : proteins) {
                if (protein.getSequence().length() < 17) continue;
                String result = gor3.predict(protein.getSequence());
                System.out.println("> " + protein.getId());
                System.out.println("AS " + protein.getSequence());
                System.out.println("PS " + result + "\n");
            }

        }else if (GOR.equals("GOR4")) {
            Gor4 gor4 = new Gor4();
            gor4.getMatrix4().loadModel4(modelPath);
            List<Protein> proteins = GorMethod.readFasta(seqPath);
            for (Protein protein : proteins) {
                String result = gor4.predict(protein.getSequence());
                if (protein.getSequence().length() < 17) continue;
                System.out.println("> " + protein.getId());
                System.out.println("AS " + protein.getSequence());
                System.out.println("PS " + result + "\n");
            }
        }

    }

}



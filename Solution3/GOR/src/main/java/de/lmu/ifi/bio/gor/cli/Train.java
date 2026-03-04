package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.io.SeclibParser;
import de.lmu.ifi.bio.gor.methods.Gor1;
import de.lmu.ifi.bio.gor.methods.Gor3;
import de.lmu.ifi.bio.gor.methods.Gor4;
import de.lmu.ifi.bio.gor.methods.GorMethod;
import de.lmu.ifi.bio.gor.model.Protein;

import java.io.IOException;
import java.util.List;


public class Train {

    public static void main(String[] args) throws IOException {

        String trainPath = null; // path to seclib training file
        String methodSelection = "1"; // <gor1|gor3|gor4>
        String modelOutputPath = null; //model file output


        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--db":
                    if (i + 1 < args.length) trainPath = args[++i];
                    break;
                case "--method":
                    if (i + 1 < args.length) methodSelection = args[++i];
                    break;
                case "--model":
                    if (i + 1 < args.length) modelOutputPath = args[++i];
                    break;

                default:
                    System.exit(1);
            }
        }

        if (trainPath == null || modelOutputPath == null) {
            System.out.println("Usage: java -jar train.jar --db <path> --method <1> --model <output>");
            return;
        }

        //richtige gor methode
        GorMethod gor = null;
        switch (methodSelection.toLowerCase()) {
            case "gor3":
                gor = new Gor3();
                break;
            case "gor4":
                gor = new Gor4();
                break;
            default:
                gor = new Gor1();
                break;
        }

        //seclib
        SeclibParser parser = new SeclibParser();
        List<Protein> data = parser.parseSeclib(trainPath); // Pfad aus den CLI-Args

        //train, get counts, no log scores
        gor.train(data);

        // 3. Speichern (schreibt die Datei, die du mir per cat gezeigt hast)
        gor.getMatrix().save(modelOutputPath);

        //Tausch ab GOR4


    }
}



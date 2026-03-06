package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.io.SeclibParser;
import de.lmu.ifi.bio.gor.methods.Evaluator;
import de.lmu.ifi.bio.gor.methods.Gor1;
import de.lmu.ifi.bio.gor.methods.Gor3;
import de.lmu.ifi.bio.gor.methods.Gor4;
import de.lmu.ifi.bio.gor.methods.GorMethod;
import de.lmu.ifi.bio.gor.model.Protein;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class CrossValidation {

    public static void main(String[] args) throws IOException {

        String dbPath = null;
        String methodSelection = "gor1";

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--db":
                    if (i + 1 < args.length) dbPath = args[++i];
                    break;
                case "--method":
                    if (i + 1 < args.length) methodSelection = args[++i];
                    break;
                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.err.println("Usage: java CrossValidation --db <file> --method <gor1|gor3|gor4>");
                    System.exit(1);
            }
        }

        if (dbPath == null) {
            System.err.println("Usage: java CrossValidation --db <file> --method <gor1|gor3|gor4>");
            System.exit(1);
        }

        // alle Proteine laden
        List<Protein> all = SeclibParser.parseSeclib(dbPath);
        System.out.println("Loaded " + all.size() + " proteins");
        System.out.println("Method: " + methodSelection.toUpperCase());
        System.out.println();

        // shuffle mit fixem seed
        Random rng = new Random(42);
        for (int i = all.size() - 1; i > 0; i--) {
            int j = rng.nextInt(i + 1);
            Protein temp = all.get(i);
            all.set(i, all.get(j));
            all.set(j, temp);
        }

        Evaluator eval = new Evaluator();
        double[] q3Scores  = new double[5];
        double[] sovScores = new double[5];
        double[] qHScores  = new double[5];
        double[] qEScores  = new double[5];
        double[] qCScores  = new double[5];

        int foldSize = all.size() / 5;

        for (int fold = 0; fold < 5; fold++) {

            // test = fold i, train = alle anderen
            int start = fold * foldSize;
            int end   = (fold == 4) ? all.size() : start + foldSize;

            List<Protein> test  = all.subList(start, end);
            List<Protein> train = new ArrayList<>();
            for (int j = 0; j < all.size(); j++) {
                if (j < start || j >= end) train.add(all.get(j));
            }

            // richtige Methode
            GorMethod gor;
            switch (methodSelection.toLowerCase()) {
                case "gor3": gor = new Gor3(); break;
                case "gor4": gor = new Gor4(); break;
                default:     gor = new Gor1(); break;
            }

            // trainieren
            gor.train(train);

            // predicten + evaluieren
            double sumQ3 = 0, sumSOV = 0, sumQH = 0, sumQE = 0, sumQC = 0;
            int n = 0;

            for (Protein p : test) {
                if (p.getSequence().length() < 17) continue;

                String predicted = gor.postProcess(gor.predict(p.getSequence()));
                String reference = p.getStructure();

                sumQ3  += eval.calculateQ3(predicted, reference);
                sumSOV += eval.calculateSOV(predicted, reference);
                sumQH  += eval.calculateQobs(predicted, reference, 'H');
                sumQE  += eval.calculateQobs(predicted, reference, 'E');
                sumQC  += eval.calculateQobs(predicted, reference, 'C');
                n++;
            }

            q3Scores[fold]  = sumQ3  / n;
            sovScores[fold] = sumSOV / n;
            qHScores[fold]  = sumQH  / n;
            qEScores[fold]  = sumQE  / n;
            qCScores[fold]  = sumQC  / n;

            System.out.printf("Fold %d | Q3: %.2f | SOV: %.2f | qH: %.2f | qE: %.2f | qC: %.2f%n",
                    fold + 1, q3Scores[fold], sovScores[fold],
                    qHScores[fold], qEScores[fold], qCScores[fold]);
        }

        // Gesamtergebnis
        System.out.println();
        System.out.println("=== 5-CV Result ===");
        System.out.printf("Q3  | mean: %.2f | std: %.2f%n", mean(q3Scores),  std(q3Scores));
        System.out.printf("SOV | mean: %.2f | std: %.2f%n", mean(sovScores), std(sovScores));
        System.out.printf("qH  | mean: %.2f | std: %.2f%n", mean(qHScores),  std(qHScores));
        System.out.printf("qE  | mean: %.2f | std: %.2f%n", mean(qEScores),  std(qEScores));
        System.out.printf("qC  | mean: %.2f | std: %.2f%n", mean(qCScores),  std(qCScores));
    }

    private static double mean(double[] arr) {
        double sum = 0;
        for (double v : arr) sum += v;
        return sum / arr.length;
    }

    private static double std(double[] arr) {
        double m = mean(arr);
        double sum = 0;
        for (double v : arr) sum += (v - m) * (v - m);
        return Math.sqrt(sum / arr.length);
    }
}
package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.io.SeclibParser;
import de.lmu.ifi.bio.gor.methods.Evaluator;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.model.ValidationResult;
import java.io.*;
import java.util.*;

public class validateGor {

    public static void main(String[] args) {

        String predictionFile = null;
        String seclibFile = null;
        String format = "txt";
        String summary= null;
        String detailed = null;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-p":
                    if (i + 1 < args.length) predictionFile = args[++i];
                    break;
                case "-r":
                    if (i + 1 < args.length) seclibFile = args[++i];
                    break;
                case "-f":
                    if (i + 1 < args.length) format = args[++i];
                    break;
                case "-s":
                    if (i + 1 < args.length) summary = args[++i];
                    break;
                case "-d":
                    if (i + 1 < args.length) detailed = args[++i];
                    break;
                //wenn nichts angegeben
                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.err.println("Usage: java -jar .jar --model <file> --seq <file> [--format <txt|html>]");
                    System.err.println("Or:    java -jar predict.jar --model <file> --maf <folder> [--format <txt|html>]");
                    System.exit(1);
            }

        }
        if (predictionFile == null || seclibFile == null) {
            System.err.println("Prediction file and Reference file are required");
            System.exit(1);
        }

        try {

            List<Protein> refList = SeclibParser.parseSeclib(seclibFile);
            Map<String, Protein> references = new HashMap<>();
            for (Protein p : refList) {
                references.put(p.getId(), p);
            }
            List<Protein> predictions = SeclibParser.parseSeclib(predictionFile);

            Evaluator eval = new Evaluator();
            List<ValidationResult> results = new ArrayList<>();

            //berechnung von q3,sov,qh,qe,qc, sov_h, sov_e, sov_c für jedes protein
            for (Protein pred : predictions) {

                String id = pred.getId();

                if (references.containsKey(id)) {

                    Protein ref = references.get(id);

                    String as = pred.getSequence();
                    String ps = pred.getStructure();
                    String ss = ref.getStructure();

                    double q3 = eval.calculateQ3(ps, ss);
                    double sov = eval.calculateSOV(ps, ss);

                    double qH = eval.calculateQobs(ps, ss, 'H');
                    double qE = eval.calculateQobs(ps, ss, 'E');
                    double qC = eval.calculateQobs(ps, ss, 'C');

                    double sovH = eval.calculateStateSOV(ps, ss, 'H');
                    double sovE = eval.calculateStateSOV(ps, ss, 'E');
                    double sovC = eval.calculateStateSOV(ps, ss, 'C');

                    results.add(new ValidationResult(id, q3, sov, qH, qE, qC, sovH, sovE, sovC, as, ps ,ss));
                }
            }

            //output
            writeDetailedOutput(detailed, results);
            writeSummaryOutput(summary, results);

        } catch (IOException e) {
            System.err.println("Error reading files: " + e.getMessage());
        }
    }

    //SUMMARYOUTPUT:
    //Statistic for protein validation
    //
    //Number of Proteins:           513
    //Mean Protein Length:          164.0
    //Sum of Protein Length:        84119
    //Sum of Predicted Positions:   75911
    //
    //Tabelle mit allen metrics: q3,qObs_H, qObs_E, qObs_C, SOV, SOV_H, SOV_E, SOV_C
    //und für jede Metric: Mean, Dev, Min, Max, Median, Quantil_25, Quantil_75, Quantil_5, Quantil_95
    private static void writeSummaryOutput(String path, List<ValidationResult> Nresults) throws IOException {
        if (path == null) return;

        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {

            // NEU: nur Proteine mit validen Q3-Werten
            List<ValidationResult> results = Nresults.stream()
                    .filter(r -> r.getQ3() >= 0)
                    .toList();

            pw.println("Statistic for protein validation");
            pw.println();
            pw.println("Number of Proteins:           " + results.size());

            //mean , sum protein length
            List<Double> lengths = results.stream().map(r -> (double) r.getReference().length()).toList();
            double sumLength = lengths.stream().mapToDouble(Double::doubleValue).sum();
            pw.printf(Locale.US, "Mean Protein Length:          %.1f%n", ValidationResult.getMean(lengths));
            pw.printf(Locale.US, "Sum of Protein Length:        %.0f%n", sumLength);

            //sum of predictet pos
            long predictedPositions = results.stream()
                    .flatMapToLong(r -> r.getPrediction().chars().filter(c -> c != '-').mapToLong(c -> 1))
                    .sum();
            pw.println("Sum of Predicted Positions:   " + predictedPositions);
            pw.println();

            //QScores und SOV scores
            printMetricLine(pw, "q3", results.stream().map(ValidationResult::getQ3).toList());
            printMetricLine(pw, "qObs_H", results.stream().map(ValidationResult::getQH).toList());
            printMetricLine(pw, "qObs_E", results.stream().map(ValidationResult::getQE).toList());
            printMetricLine(pw, "qObs_C", results.stream().map(ValidationResult::getQC).toList());
            pw.println();
            printMetricLine(pw, "SOV", results.stream().map(ValidationResult::getSov).toList());
            printMetricLine(pw, "SOV_H", results.stream().map(ValidationResult::getSovH).toList());
            printMetricLine(pw, "SOV_E", results.stream().map(ValidationResult::getSovE).toList());
            printMetricLine(pw, "SOV_C", results.stream().map(ValidationResult::getSovC).toList());
        }
    }

    private static void printMetricLine(PrintWriter pw, String name, List<Double> values) {

        double mean = ValidationResult.getMean(values);
        double dev = ValidationResult.getStdDev(values, mean);

        // Min/Max nur von validen Werten (>= 0)
        List<Double> valid = values.stream().filter(v -> v >= 0).toList();
        double min = valid.stream().mapToDouble(v -> v).min().orElse(0.0);
        double max = valid.stream().mapToDouble(v -> v).max().orElse(0.0);

        pw.printf(Locale.US, "%-8s: Mean: %6.1f   Dev: %6.1f   Min: %6.1f   Max: %6.1f   Median: %6.1f   " +
                        "Quantil_25: %6.1f   Quantil_75: %6.1f   Quantil_5: %6.1f   Quantil_95: %6.1f%n",
                name, mean, dev, min, max,
                ValidationResult.getMedian(values),
                ValidationResult.getQuantile(values, 0.25),
                ValidationResult.getQuantile(values, 0.75),
                ValidationResult.getQuantile(values, 0.05),
                ValidationResult.getQuantile(values, 0.95));
    }

    //DETAILEDOUTPUT:
    //
    //>id Q3  SOV  QH  QE  QC  SOV_H  SOV_E  SOV_C
    //AS GLPDVQTKEAAAAAKKKGGLPRRTKLAAAAAAAA
    //PS --------HHHHHHHHHEEEEECCCCCC--------
    //SS CHHHHHCHHHHHHHHHEEEEECCEEEECCCHHHH

    private static void writeDetailedOutput(String path, List<ValidationResult> nresults) throws IOException {
        if (path == null) return;

        try  (PrintWriter pw = new PrintWriter(new FileWriter(path))) {

            List<ValidationResult> results = nresults.stream()
                    .filter(r -> r.getQ3() >= 0)
                    .toList();

            for (ValidationResult r : results) {
                //head >id Q3 SOV QH QE QC SOV_H SOV_E SOV_C
                pw.printf(Locale.US, "> %-10s %s %s %s %s %s %s %s %s%n",
                        r.getId(),
                        formatScore(r.getQ3()),
                        formatScore(r.getSov()),
                        formatScore(r.getQH()),
                        formatScore(r.getQE()),
                        formatScore(r.getQC()),
                        formatScore(r.getSovH()),
                        formatScore(r.getSovE()),
                        formatScore(r.getSovC())
                );

                pw.println("AS " + r.getSequence());
                pw.println("PS " + r.getPrediction());
                pw.println("SS " + r.getReference());

                pw.println();
            }
        }
    }

    //gibt bei -1.0 den Strich - aus
    private static String formatScore(double value) {
        if (value < 0) {
            return " - ";
        }
        return String.format(Locale.US, "%5.1f", value);
    }
}

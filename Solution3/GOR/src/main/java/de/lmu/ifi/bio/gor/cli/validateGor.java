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
        String summary = null;
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
                default:
                    System.err.println("Unknown argument: " + args[i]);
                    System.exit(1);
            }
        }

        if (predictionFile == null || seclibFile == null) {
            System.err.println("Prediction file and Reference file are required");
            System.exit(1);
        }

        try {
            List<Protein> refList = SeclibParser.parseSeclib(seclibFile);
            Map<String, Protein> predMap = parsePrd(predictionFile);

            Evaluator eval = new Evaluator();
            List<ValidationResult> results = new ArrayList<>();

            for (Protein ref : refList) {
                String id = ref.getId();
                Protein pred = predMap.get(id);

                String ss = ref.getStructure();
                String as = ref.getSequence();

                if (pred == null) {
                    String emptyPs = "-".repeat(as.length());
                    results.add(new ValidationResult(id, 0.0, 0.0, Double.NaN, Double.NaN, Double.NaN,
                            Double.NaN, Double.NaN, Double.NaN, as, emptyPs, ss));
                    continue;
                }

                String ps = pred.getStructure();

                double q3 = eval.calculateQ3(ps, ss);
                double sov = eval.calculateSOV(ps, ss);
                double qH = eval.calculateQobs(ps, ss, 'H');
                double qE = eval.calculateQobs(ps, ss, 'E');
                double qC = eval.calculateQobs(ps, ss, 'C');
                double sovH = eval.calculateStateSOV(ps, ss, 'H');
                double sovE = eval.calculateStateSOV(ps, ss, 'E');
                double sovC = eval.calculateStateSOV(ps, ss, 'C');

                results.add(new ValidationResult(id, q3, sov, qH, qE, qC, sovH, sovE, sovC, as, ps, ss));
            }

            writeDetailedOutput(detailed, results);
            writeSummaryOutput(summary, results);

        } catch (IOException e) {
            System.err.println("Error reading files: " + e.getMessage());
        }
    }

    private static void writeSummaryOutput(String path, List<ValidationResult> results) throws IOException {
        if (path == null) return;

        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {
            pw.println("Statistic for protein validation");
            pw.println();
            pw.println("Number of Proteins:           " + results.size());

            List<Double> lengths = results.stream().map(r -> (double) r.getReference().length()).toList();
            double sumLength = lengths.stream().mapToDouble(Double::doubleValue).sum();
            pw.printf(Locale.US, "Mean Protein Length:          %.1f%n", ValidationResult.getMean(lengths));
            pw.printf(Locale.US, "Sum of Protein Length:        %.0f%n", sumLength);

            long predictedPositions = results.stream()
                    .flatMapToLong(r -> r.getPrediction().chars().filter(c -> c != '-').mapToLong(c -> 1))
                    .sum();
            pw.println("Sum of Predicted Positions:   " + predictedPositions);
            pw.println();

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

    private static void writeDetailedOutput(String path, List<ValidationResult> results) throws IOException {
        if (path == null) return;

        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {
            for (ValidationResult r : results) {
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

    //Abfangen der -1.0 Fehler von SOV und Qobs --> -
    //Umformung in 1 Nachkommastelle
    private static String formatScore(double value) {
        if (value < 0) return " - ";
        return String.format(Locale.US, "%5.1f", value);
    }

    //gibt aus: Map<ID, AS+PS>
    //liesst ein: > 1234
    //            AS GPK....
    //            SS CCC....
    private static Map<String, Protein> parsePrd(String path) throws IOException {
        Map<String, Protein> map = new LinkedHashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line, id = null, as = "", ps = "";
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("> ")) {
                    if (id != null) map.put(id, new Protein(id, as, ps));
                    id = line.substring(2).trim().split("\\s+")[0];
                    as = "";
                    ps = "";
                } else if (line.startsWith("AS ")) {
                    as = line.substring(3).trim();
                } else if (line.startsWith("PS ")) {
                    ps = line.substring(3).trim();
                }
            }
            if (id != null) map.put(id, new Protein(id, as, ps));
        }
        return map;
    }
}
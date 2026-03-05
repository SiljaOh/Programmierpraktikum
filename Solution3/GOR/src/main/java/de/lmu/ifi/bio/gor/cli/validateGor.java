package de.lmu.ifi.bio.gor.cli;

import de.lmu.ifi.bio.gor.methods.Evaluator;
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
            //parse fasta
            Map<String, String[]> predictions = parseValidationFasta(predictionFile);
            Map<String, String[]> references = parseValidationFasta(seclibFile);

            Evaluator eval = new Evaluator();
            List<ValidationResult> results = new ArrayList<>();

            //berechnung von q3,sov,qh,qe,qc für jedes protein
            for (String id : predictions.keySet()) {
                if (references.containsKey(id)) {

                    String[] predData = predictions.get(id);
                    String[] refData = references.get(id);

                    String as = predData[0];
                    String ps = predData[1];
                    String ss = refData[1];

                    double q3 = eval.calculateQ3(ps, ss);
                    double sov = eval.calculateSOV(ps, ss);
                    double qH = eval.calculateQobs(ps, ss, 'H');
                    double qE = eval.calculateQobs(ps, ss, 'E');
                    double qC = eval.calculateQobs(ps, ss, 'C');

                    results.add(new ValidationResult(id, q3, sov, qH, qE, qC, as, ps ,ss));
                }
            }

            //output
            writeDetailedOutput(detailed, results);
            writeSummaryOutput(summary, results);

            System.out.println("Validation completed");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void writeSummaryOutput(String path, List<ValidationResult> results) throws IOException {
        if (path == null) return;

        List<Double> q3s = results.stream().map(r -> r.getQ3()).toList();
        List<Double> sovs = results.stream().map(r -> r.getSov()).toList();

        double meanQ3 = ValidationResult.getMean(q3s);
        double stdQ3 = ValidationResult.getStdDev(q3s, meanQ3);

        try (PrintWriter pw = new PrintWriter(new FileWriter(path))) {
            pw.println("Summary Statistics:");
            pw.printf("Q3:  Mean = %.2f, StdDev = %.2f%n", meanQ3, stdQ3);
            pw.printf("SOV: Mean = %.2f%n", ValidationResult.getMean(sovs));
            printMetricLine(pw, "qH", results.stream().map(ValidationResult::getQH).toList());
            printMetricLine(pw, "qE", results.stream().map(ValidationResult::getQE).toList());
            printMetricLine(pw, "qC", results.stream().map(ValidationResult::getQC).toList());
        }
    }

    private static void writeDetailedOutput(String path, List<ValidationResult> results) throws IOException {
        //TODO
    }

    //TODO
    private static Map<String, String[]> parseValidationFasta(String file) throws IOException {
        Map<String, String[]> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {

                    String id = line.substring(1).split("\\s+")[0];

                    String as = br.readLine(); // Aminosäuresequenz
                    String str = br.readLine(); // Struktur PS oder SS

                    if (id != null && as != null && str != null) {
                        //array: [0]=AS, [1]=Struktur
                        map.put(id, new String[]{as.trim(), str.trim()});
                    }
                }
            }
        }
        return map;
    }
    private static void printMetricLine(PrintWriter pw, String name, List<Double> values) {
        double mean = ValidationResult.getMean(values);
        double std = ValidationResult.getStdDev(values, mean);
        pw.printf(Locale.US, "%-4s: Mean = %6.2f, StdDev = %6.2f%n", name, mean, std);
    }

}

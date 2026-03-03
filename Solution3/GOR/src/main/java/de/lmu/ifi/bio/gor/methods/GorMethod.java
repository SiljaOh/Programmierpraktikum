package de.lmu.ifi.bio.gor.methods;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.model.GorMatrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public abstract class GorMethod {
    //Vererbungsklasse für Gor1,3,4,5
    //alle Methoden die alle Gor benötigen
    //3 x 17 x 20 matrix
    protected GorMatrix matrix;

    public abstract void train(List<Protein> data);

    public abstract String predict(String sequence);

    public abstract void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result);

    public GorMatrix getMatrix() {
        return this.matrix; // Stellt sicher, dass Predict auf die Matrix zugreifen kann
    }


    public static List<Protein> readFasta(String filePath) throws IOException {
        List<Protein> proteins = new ArrayList<Protein>();
        String id = null;
        StringBuilder sequence = new StringBuilder();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {

                    if (id != null && sequence.length() > 0) {
                        proteins.add(new Protein(id, sequence.toString(), " ".repeat(sequence.length())));
                        sequence.setLength(0); // StringBuilder leeren für das nächste Protein
                    }

                    id = line.substring(1).trim();
                } else {
                    sequence.append(line.replaceAll("\\s+", ""));
                }
            }
            if (id != null && sequence.length() > 0) {
                proteins.add(new Protein(id, sequence.toString(), " ".repeat(sequence.length())));
            }
        }

        if (proteins.isEmpty()) {
            throw new IOException("Datei '" + filePath + "' enthält keinen gültigen FASTA-Header (startet nicht mit '>').");
        }
        return proteins;
    }

    //sliding window anpassen
    public String slidingWindow(Protein protein, boolean isTraining) {
        String sequence = protein.getSequence();
        String structure = protein.getStructure();
        StringBuilder result = new StringBuilder();

        int windowSize = 17;
        int m = windowSize / 2;

        for (int i = m; i < sequence.length() - m; i++) {

            // Fenster ausschneiden (17 Zeichen lang)
            // substring(start, end) -> end ist exklusiv, daher + m + 1
            String window = sequence.substring(i - m, i + m + 1);

            char target = ' ';
            if (isTraining) {
                target = structure.charAt(i);
            }

            handleWindow(window, i, target, isTraining, result);
        }

        // predict: ignorierte Ränder mit Platzhaltern auffüllen
        // damit die Länge des Resultats stimmt.
        if (!isTraining) {
            String paddedResult = "-".repeat(m) + result.toString() + "-".repeat(m);
            return paddedResult;
        }

        return result.toString();
    }

}

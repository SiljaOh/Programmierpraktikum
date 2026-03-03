package de.lmu.ifi.bio.gor.methods;
import de.lmu.ifi.bio.gor.model.Protein;
import de.lmu.ifi.bio.gor.model.GorMatrix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

public abstract class GorMethod {
    //Vererbungsklasse für Gor1,3,4,5
    //alle Methoden die alle Gor benötigen
    //3 x 17 x 20 matrix
    protected GorMatrix matrix = new GorMatrix();

    public abstract void train(List<Protein> data);

    public abstract String predict(String sequence);

    public abstract void handleWindow(String window, int position, char target, boolean isTraining, StringBuilder result);

    public GorMatrix getMatrix() {
        return this.matrix; // Stellt sicher, dass Predict auf die Matrix zugreifen kann
    }


    protected String applyPadding(String sequence) {
        String padding = "XXXXXXXX"; //17er Fenster; 8 links, Zentrum, 8 rechts
        return padding + sequence + padding;
    }

    public static Protein readFasta(String filePath) throws IOException {
        String id = null;
        StringBuilder sequence = new StringBuilder();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    // Falls schon eine ID da ist, haben wir es mit einer Multi-Fasta Datei zu tun.
                    // Hier brechen wir ab oder verarbeiten nur das erste Protein.
                    if (id != null) break;

                    // ID extrahieren (nimmt alles nach dem '>')
                    id = line.substring(1).trim();
                } else {
                    // Sequenz-Zeilen bereinigen (Leerzeichen entfernen, falls vorhanden)
                    sequence.append(line.replaceAll("\\s+", ""));
                }
            }
        }

        // Validierung: Fail-Fast
        if (id == null) {
            throw new IOException("Datei '" + filePath + "' enthält keinen gültigen FASTA-Header (startet nicht mit '>').");
        }
        if (sequence.length() == 0) {
            throw new IOException("Datei '" + filePath + "' enthält keine Sequenzdaten für ID: " + id);
        }

        // Protein mit ID, Sequenz und (wichtig!) der Platzhalter-Struktur zurückgeben
        return new Protein(id, sequence.toString(), " ".repeat(sequence.length()));
    }


    //sliding window anpassen
    public String slidingWindow(Protein protein, boolean isTraining) {

        String sequence = protein.getSequence();
        String structure = protein.getStructure();
        String paddedSequence = applyPadding(protein.getSequence());

        StringBuilder result = new StringBuilder();

        for(int i = 0; i < sequence.length(); i++) {
            //windows ausschneiden
            String window = paddedSequence.substring(i, i + 17); //Window size = 17

            //struktur im zentrum
            char target = ' ';

            if (isTraining) {
                target = structure.charAt(i);
            }

            //window individuell weiterverarbeiten
            handleWindow(window, i, target, isTraining, result);
        }
        return result.toString();
    }
}

package de.lmu.ifi.bio.gor.io;

import de.lmu.ifi.bio.gor.model.Protein;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class SeclibParser {
    //Verbindungsstelle zu get_seclib.py
    //Seclib: > ID
    //          AS
    //          SS

    public List<Protein> parseSeclib(String filePath) throws IOException {

        List<Protein> proteinList = new ArrayList<>();
        File file = new File(filePath);

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            String id = null;
            String as = null;
            String ss = null;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {

                    if (id != null && as != null && ss != null) {
                        proteinList.add(new Protein(id, as, ss));
                    }
                    id = line.substring(1).trim(); //ID ohne >
                    as = null; // Reset
                    ss = null;

                } else if (line.startsWith("AS ")) {
                    as = line.substring(3).trim();
                } else if (line.startsWith("SS ")) {
                    ss = line.substring(3).trim();
                }

                //letztes protein
                if (id != null && as != null && ss != null) {
                    proteinList.add(new Protein(id, as, ss));
                }
            }

            return proteinList;
        }
    }
}






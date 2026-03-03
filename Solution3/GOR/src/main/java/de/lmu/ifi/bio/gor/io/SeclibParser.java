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

                String trimmedLine = line.trim();
                if (trimmedLine.isEmpty()) continue;

                if (trimmedLine.startsWith(">")) {

                    if (id != null && as != null && ss != null) {
                        proteinList.add(new Protein(id, as, ss));
                    }
                    // Reset für das neue Protein
                    id = trimmedLine.substring(1).trim();
                    as = null;
                    ss = null;
                } else if (trimmedLine.startsWith("AS ")) {

                    as = trimmedLine.substring(3).trim().replaceAll("[^A-Z]", "");
                } else if (trimmedLine.startsWith("SS ")) {

                    ss = trimmedLine.substring(3).trim().replaceAll("[^A-Z]","");
                }
            }


            if (id != null && as != null && ss != null) {
                proteinList.add(new Protein(id, as, ss));
            }
        }
        return proteinList;

    }
}






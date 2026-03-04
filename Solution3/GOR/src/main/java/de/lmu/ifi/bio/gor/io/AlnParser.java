package de.lmu.ifi.bio.gor.io;

import de.lmu.ifi.bio.gor.model.Alignment;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class AlnParser {

    //format: > 9apia (Zielprotein)
    //        AS AKHSEAFTVNFGDTEE... (AS sequence)
    //        SS CCCHHHHHHHHHHCCC... (sekundärstruktur
    //         1 AK-S-AFTVN-----E... (nummerierte AS sequenzen mit Gaps)
    //         .......

    public static Alignment parseAlignment(String filePath) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                if (line.startsWith(">")) {
                    String id = line.substring(1).trim();
                    String asSeq = null;
                    String ssSeq = null;


                    String nextLine;
                    while ((nextLine = br.readLine()) != null) {
                        nextLine = nextLine.trim();
                        if (nextLine.startsWith("AS")) {
                            asSeq = nextLine.substring(3).trim();
                        } else if (nextLine.startsWith("SS")) {
                            ssSeq = nextLine.substring(3).trim();
                            break;
                        }
                    }

                    Alignment aln = new Alignment(id, asSeq, ssSeq);

                    //homologues
                    while ((line = br.readLine()) != null) {
                        line = line.trim();
                        if (line.isEmpty()) continue;

                        String[] parts = line.split("\\s+", 2);
                        if (parts.length == 2) {
                            aln.addHomologues(parts[1]);
                        }
                    }

                    return aln;
                }
            }
        }
        return null;
    }
}



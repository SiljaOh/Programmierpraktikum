package de.lmu.ifi.bio.gor.io;

import de.lmu.ifi.bio.gor.model.Alignment;

import java.io.IOException;

public class AlnParser {

    //format: > 9apia (Zielprotein)
    //        AS AKHSEAFTVNFGDTEE... (AS sequence)
    //        SS CCCHHHHHHHHHHCCC... (sekundärstruktur
    //         1 AK-S-AFTVN-----E... (nummerierte AS sequenzen mit Gaps)
    //         .......

    public static Alignment parseAlignment(String alignment) throws IOException {

    }

}

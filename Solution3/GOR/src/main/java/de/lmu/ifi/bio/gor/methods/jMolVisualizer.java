package de.lmu.ifi.bio.gor.methods;

import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class jMolVisualizer {


    public static void writeHtml(String id, String sequence, String realSS, String predictedSS, String outputPath) throws IOException {

        //findet die falsch predicteten
        StringBuilder wrong = new StringBuilder();
        for (int i = 0; i < realSS.length(); i++) {
            if (realSS.charAt(i) != predictedSS.charAt(i)) {
                if (wrong.length() > 0) {
                    //in JMOL: select resno=5 or resno=12 or resno=47
                    wrong.append(" or resno=");
                }
                //welche stelle, 0 anfnag
                wrong.append(i + 1);
            }
        }


        PrintWriter pw = new PrintWriter(new FileWriter(outputPath));
        pw.println("<!DOCTYPE html><html><head><title>" + id + "</title>");
        pw.println("<script src='https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js'></script>");
        pw.println("</head><body style='background:#111;color:#eee;font-family:monospace;padding:20px'>");
        //ID als titel
        pw.println("<h2>" + id + "</h2>");
        //2 nebeneinander zum vergleich
        pw.println("<div style='display:flex;gap:20px'>");
        pw.println("<div><b>Echte Struktur</b><div id='r'></div></div>");
        pw.println("<div><b>Vorhergesagt</b><div id='p'></div></div>");
        pw.println("</div>");
        //html text ausgabe
        pw.println("<pre style='margin-top:20px'>");
        pw.println("AS   " + sequence);
        pw.println("REAL " + realSS);
        pw.print  ("PRED ");
        for (int i = 0; i < predictedSS.length(); i++) {
            boolean isWrong = i < realSS.length() && realSS.charAt(i) != predictedSS.charAt(i);
            if (isWrong) {
                //flach in pink
                pw.print("<span style='color:pink'>" + predictedSS.charAt(i) + "</span>");
            }
            else {
                pw.print(predictedSS.charAt(i));
            }
        }
        pw.println("\n</pre>");

        // JMol Script
        pw.println("<script>");
        pw.println("var Info = {width:400,height:400,use:'HTML5',color:'#111',");
        pw.println("  j2sPath:'https://chemapps.stolaf.edu/jmol/jsmol/j2s',script:''};");
        //id loaden
        String colorScript = "load " + id + "; background black; cartoons on; wireframe off;" + "select helix; color blue; select sheet; color yellow; select coil; color grey;";

        // Viewer echte Struktur
        pw.println("Info.script = \"" + colorScript + "\";");
        pw.println("var r = Jmol.getApplet('r', Info);");
        pw.println("document.getElementById('r').innerHTML = Jmol.getAppletHtml(r);");

        // Viewer vorhergesagt, falsche positionen pink
        String predScript = colorScript;
        if (wrong.length() > 0) {
            //zusammenführen mit starr uhd farbe
            predScript += "select resno=" + wrong + "; color pink;";
        }
        pw.println("Info.script = \"" + predScript + "\";");
        pw.println("var p = Jmol.getApplet('p', Info);");
        pw.println("document.getElementById('p').innerHTML = Jmol.getAppletHtml(p);");
        pw.println("</script></body></html>");

        pw.close();

        //öffnen
        // nach pw.close():
        Desktop.getDesktop().browse(new File(outputPath).toURI());


    }
}

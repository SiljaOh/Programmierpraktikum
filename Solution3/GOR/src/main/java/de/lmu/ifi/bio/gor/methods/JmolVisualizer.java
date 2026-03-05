package de.lmu.ifi.bio.gor.methods;

import java.awt.Desktop;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class JmolVisualizer {

    public static void show(String id, String sequence, String realSS, String predictedSS) throws IOException {

        File temp = File.createTempFile("gor_" + id, ".html");
        temp.deleteOnExit();

        // falsche Positionen sammeln (1-basiert für JMol)
        List<Integer> wrongList = new ArrayList<>();
        for (int i = 0; i < Math.min(realSS.length(), predictedSS.length()); i++) {
            if (realSS.charAt(i) != predictedSS.charAt(i)) {
                wrongList.add(i + 1);
            }
        }
        //Jmol falsche posis
        String wrongSelect = "";
        if (!wrongList.isEmpty()) {
            StringBuilder sb = new StringBuilder(" select ");
            for (int i = 0; i < wrongList.size(); i++) {
                if (i > 0) {
                    sb.append(" or ");
                }
                sb.append("resno=").append(wrongList.get(i));
            }
            sb.append("; color orange;");
            wrongSelect = sb.toString();
        }

        String colorScript = "load =" + id + "; background black; cartoons on; wireframe off;" + "select helix; color blue; select sheet; color yellow; select coil; color grey;";

        PrintWriter pw = new PrintWriter(new FileWriter(temp));
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
        //html text ausgabe Sequenzen
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

        // JMol
        pw.println("<script>");
        pw.println("var Info = {width:400,height:400,use:'HTML5',color:'#111',");
        pw.println("  j2sPath:'https://chemapps.stolaf.edu/jmol/jsmol/j2s',script:''};");

        // Viewer echte Struktur
        pw.println("Info.script = \"" + colorScript + "\";");
        pw.println("var r = Jmol.getApplet('r', Info);");
        pw.println("document.getElementById('r').innerHTML = Jmol.getAppletHtml(r);");

        pw.println("Info.script = \"" + colorScript + wrongSelect + "\";");
        pw.println("var p = Jmol.getApplet('p', Info);");
        pw.println("document.getElementById('p').innerHTML = Jmol.getAppletHtml(p);");
        pw.println("</script></body></html>");

        pw.close();
        Desktop.getDesktop().browse(temp.toURI());
    }
}

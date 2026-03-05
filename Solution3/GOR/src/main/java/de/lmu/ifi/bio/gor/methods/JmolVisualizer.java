package de.lmu.ifi.bio.gor.methods;

import java.awt.Desktop;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class JmolVisualizer {

    // alle Proteine speichern
    private static final List<String[]> proteins = new ArrayList<>();

    public static void addProtein(String id, String sequence, String realSS, String predictedSS) {
        proteins.add(new String[]{id, sequence, realSS, predictedSS});
    }

    public static void showAll() throws IOException {
        File temp = File.createTempFile("gor_all", ".html");
        temp.deleteOnExit();

        PrintWriter pw = new PrintWriter(new FileWriter(temp));
        pw.println("<!DOCTYPE html><html><head><title>GOR Visualisierung</title>");
        pw.println("<script src='https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js'></script>");
        pw.println("</head><body style='background:#111;color:#eee;font-family:monospace;padding:20px'>");

        // Liste aller Proteine als Buttons
        pw.println("<h2>GOR Visualisierung</h2>");
        pw.println("<div style='margin-bottom:20px;display:flex;flex-wrap:wrap;gap:8px'>");
        for (String[] p : proteins) {
            pw.println("<button onclick='show(\"" + p[0] + "\")' style='background:#333;color:#eee;border:none;padding:6px 12px;cursor:pointer'>" + p[0] + "</button>");
        }
        pw.println("</div>");

        pw.println("<div id='content'></div>");

        // alle protein daten
        pw.println("<script>");
        pw.println("var data = {};");
        for (String[] p : proteins) {
            String id = p[0], seq = p[1], real = p[2], pred = p[3];

            // falsche positionen
            List<Integer> wrongList = new ArrayList<>();
            for (int i = 0; i < Math.min(real.length(), pred.length()); i++) {
                if (real.charAt(i) != pred.charAt(i)) wrongList.add(i + 1);
            }
            String wrongSelect = "";
            if (!wrongList.isEmpty()) {
                StringBuilder sb = new StringBuilder(" select ");
                for (int i = 0; i < wrongList.size(); i++) {
                    if (i > 0) sb.append(" or ");
                    sb.append("resno=").append(wrongList.get(i));
                }
                sb.append("; color pink;");
                wrongSelect = sb.toString();
            }

            // pred sequenz mit pink
            StringBuilder predHtml = new StringBuilder();
            for (int i = 0; i < pred.length(); i++) {
                boolean isWrong = i < real.length() && real.charAt(i) != pred.charAt(i);
                if (isWrong) predHtml.append("<span style='color:pink'>").append(pred.charAt(i)).append("</span>");
                else         predHtml.append(pred.charAt(i));
            }

            pw.println("data['" + id + "'] = {");
            pw.println("  seq: '" + seq + "',");
            pw.println("  real: '" + real + "',");
            pw.println("  predHtml: '" + predHtml + "',");
            pw.println("  wrongSelect: '" + wrongSelect + "'");
            pw.println("};");
        }

        // show funktion
        pw.println("function show(id) {");
        pw.println("  var d = data[id];");
        pw.println("  var colorScript = 'load =' + id + '; background black; cartoons on; wireframe off;" +
                "select helix; color red; select sheet; color yellow; select coil; color grey;';");
        //überschriften und so
        pw.println("  document.getElementById('content').innerHTML = " +
                "'<h3>' + id + '</h3>" +
                "<div style=\"display:flex;gap:20px\">" +
                "<div><b>Echte Struktur</b><div id=\\'r\\'></div></div>" +
                "<div><b>Vorhergesagt</b><div id=\\'p\\'></div></div>" +
                "</div>" +
                "<pre style=\"margin-top:20px\">AS   ' + d.seq + '\\nREAL ' + d.real + '\\nPRED ' + d.predHtml + '</pre>';");
        pw.println("  var Info = {width:400,height:400,use:'HTML5',color:'#111'," +
                "j2sPath:'https://chemapps.stolaf.edu/jmol/jsmol/j2s',script:''};");
        pw.println("  Info.script = colorScript;");
        pw.println("  var r = Jmol.getApplet('r', Info);");
        pw.println("  document.getElementById('r').innerHTML = Jmol.getAppletHtml(r);");
        pw.println("  Info.script = colorScript + d.wrongSelect;");
        pw.println("  var p = Jmol.getApplet('p', Info);");
        pw.println("  document.getElementById('p').innerHTML = Jmol.getAppletHtml(p);");
        pw.println("}");

        // erstes protein direkt anzeigen
        if (!proteins.isEmpty()) {
            pw.println("show('" + proteins.get(0)[0] + "');");
        }

        pw.println("</script></body></html>");
        pw.close();

        Desktop.getDesktop().browse(temp.toURI());
    }
}
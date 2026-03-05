package de.lmu.ifi.bio.gor.cli;

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
    }

}

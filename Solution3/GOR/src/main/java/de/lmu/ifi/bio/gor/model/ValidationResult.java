package de.lmu.ifi.bio.gor.model;

import java.util.List;

public class ValidationResult {

   private String id;
   private String sequence;
   private String prediction;
   private String reference;
   private double q3;
   private double sov;
   private double qH;
   private double qE;
   private double qC;

   public ValidationResult(String id, double q3, double sov, double qH, double qE, double qC,
                           String sequence, String prediction, String reference) {
      this.id = id;
      this.q3 = q3;
      this.sov = sov;
      this.qH = qH;
      this.qE = qE;
      this.qC = qC;
      this.sequence = sequence;
      this.prediction = prediction;
      this.reference = reference;
   }
   //TODO SOV per state
   public double getQ3() { return q3; }
   public double getSov() { return sov; }
   public double getQH() { return qH; }
   public double getQE() { return qE; }
   public double getQC() { return qC; }
   public String getId() { return id; }
   public String getSequence() { return sequence; }
   public String getPrediction() { return prediction; }
   public String getReference() { return reference; }

   //Mean
   public static double getMean(List<Double> values) {
      return values.stream().mapToDouble(d -> d).average().orElse(0.0);
   }

   //Standard Deviation
   public static double getStdDev(List<Double> values, double mean) {
      double temp = 0;
      for (double d : values) temp += (d - mean) * (d - mean);
      return Math.sqrt(temp / values.size());
   }
}

package de.lmu.ifi.bio.gor.model;

import java.util.ArrayList;
import java.util.Collections;
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
   private double sovH;
   private double sovE;
   private double sovC;

   public ValidationResult(String id, double q3, double sov, double qH, double qE, double qC,
                           double sovH, double sovE, double sovC,
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
      this.sovH = sovH;
      this.sovE = sovE;
      this.sovC = sovC;
   }

   public double getQ3() { return q3; }
   public double getSov() { return sov; }
   public double getQH() { return qH; }
   public double getQE() { return qE; }
   public double getQC() { return qC; }
   public String getId() { return id; }
   public String getSequence() { return sequence; }
   public String getPrediction() { return prediction; }
   public String getReference() { return reference; }
   public double getSovH() { return sovH; }
   public double getSovE() { return sovE; }
   public double getSovC() { return sovC; }

   //Mean
   public static double getMean(List<Double> values) {
      return values.stream().filter(d -> d >= 0).mapToDouble(d -> d).average().orElse(0.0);
   }

   //Standard Deviation
   public static double getStdDev(List<Double> values, double mean) {
      List<Double> valid = values.stream().filter(d -> d >= 0).toList();
      if (valid.isEmpty()) return 0.0;
      double temp = 0;
      for (double d : valid) temp += (d - mean) * (d - mean);
      return Math.sqrt(temp / valid.size());
   }

   //median = 50% Quantile
   public static double getMedian(List<Double> values) {
      return getQuantile(values, 0.5);
   }

   //Quantile
   public static double getQuantile(List<Double> values, double quantile) {
      List<Double> validValues = new ArrayList<>(values.stream().filter(d -> d >= 0).toList());
      if (validValues.isEmpty()) return 0.0;

      Collections.sort(validValues);

      int index = (int) Math.ceil(quantile * validValues.size()) - 1;
      return validValues.get(Math.max(0, index));
   }
}

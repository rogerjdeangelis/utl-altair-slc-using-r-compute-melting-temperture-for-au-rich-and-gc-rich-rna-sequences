%let pgm=utl-altair-slc-using-r-compute-melting-temperture-for-au-rich-and-gc-rich-rna-sequences;

%stop-submission;

Altair slc using r compute melting temperture for au rich and gc rich rna sequences

Too long to post on a list, see github
https://github.com/rogerjdeangelis/utl-altair-slc-using-r-compute-melting-temperture-for-au-rich-and-gc-rich-rna-sequences

[SIEMENS FORUM](https://support.industry.siemens.com/cs/document/109780236)
https://support.industry.siemens.com/cs/document/109780236/how-to-configure-rna?dti=0&lc=en-US

WHAT WE WANT TO COMPUTE MELTING TEMPERATURE FOR RNA SEQUENCES

   Given RNA Sequences compare melting temperatures (unfolding the helix?)

   AU-rich (20% GC)  = AAAAUUUUAAAAUUUUAAAAUUUU
   GC-rich (80% GC)  = GGGGCCCCGGGGCCCCGGGGCCCC


   1. AU-rich RNA melts at: 57.5°C
      GC-rich RNA melts at: 96.9°C
      Difference: 39.4°C

   2. The GC-rich sequence requires MUCH higher temperature to melt
      because each G-C pair has 3 hydrogen bonds (vs 2 in A-U pairs)

   Melting Temperature (Tm) = The temperature where the RNA structure
   is unfolded,

   RNA folding is the process where a linear chain of RNA nucleotides curls up into a
   specific 3D shape by forming base pairs with itself.

   RNA unfolding in bacteria can KILL them or trigger their DEATH -
   but not directly from the RNA unfolding itself. Instead, RNA unfolding
   disrupts essential bacterial processes, leading to cell death

CONTENTS

  1  Preparation
  2  R RNA process
  3  R RNA output
  4  R log

Python does support viennaRNA, but it seems to requite coda python and
seems directed toward linux?

/*                                        _   _
/ |  _ __  _ __ ___ _ __   __ _ _ __ __ _| |_(_) ___  _ __
| | | `_ \| `__/ _ \ `_ \ / _` | `__/ _` | __| |/ _ \| `_ \
| | | |_) | | |  __/ |_) | (_| | | | (_| | |_| | (_) | | | |
|_| | .__/|_|  \___| .__/ \__,_|_|  \__,_|\__|_|\___/|_| |_|
    |_|            |_|
*/

 Step-by-Step Installation Guide for viennaRNA

 1. Download the Binary Package

  First, navigate to the
  official ViennaRNA website and download the correct installer for Windows. The official source is
  the University of Vienna's Theoretical Biochemistry group.

  Official Website:
  https://www.tbi.univie.ac.at/RNA/

  What to look for: On the website, locate the "Download"
  section and choose the binary installer for Windows (e.g., a .exe file). Ensure you select the
  version matching your system architecture (32-bit or 64-bit). As of late 2024, the latest stable
  release is version 2.7.0,


 2. Run the Installer

  Once the download is complete:

  Run
  the .exe file. You may need to grant administrator permissions if a User Account Control (UAC)
  prompt appears.

  Follow the on-screen instructions. The installer will guide you through the
  setup process. The default installation directory is typically something like C:\Program
  Files\ViennaRNA.


 3. Add to System PATH (Crucial Step)

  Find the installation path: Locate the bin folder inside your
  ViennaRNA installation directory (e.g., C:\Program Files\ViennaRNA\bin). Inside, you should see
  files like RNAfold.exe.

  Open System Properties: Press Win + R, type sysdm.cpl, and press
  Enter. Go to the Advanced tab and click Environment Variables.

  Edit the PATH: Under "System
  variables", scroll down and select the Path variable, then click Edit.

  Add the path: Click
  New and paste the full path to the bin folder (e.g., C:\Program Files (x86)\ViennaRNA Package).

  Save:
  Click OK on all open windows.


 4. Verify the Installation

  To confirm everything is working:

  Press Win + R, type cmd, and press Enter to open the Command Prompt.

  Type RNAfold

   Input string (upper or lower case); @ to quit
   ....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8

   Type @ to exit

   FYI: If you open task manager and restart 'Windows Explorer' the path to viennaRNA
   eill be efftuated, Yow will need to close and reopen any app you are using,

/*___
|___ \   _ __   _ __ _ __   __ _   _ __  _ __ ___   ___ ___  ___ ___
  __) | | `__| | `__| `_ \ / _` | | `_ \| `__/ _ \ / __/ _ \/ __/ __|
 / __/  | |    | |  | | | | (_| | | |_) | | | (_) | (_|  __/\__ \__ \
|_____| |_|    |_|  |_| |_|\__,_| | .__/|_|  \___/ \___\___||___/___/
                                  |_|
 _ __  _ __ ___   ___ ___  ___ ___
| `_ \| `__/ _ \ / __/ _ \/ __/ __|
| |_) | | | (_) | (_|  __/\__ \__ \
| .__/|_|  \___/ \___\___||___/___/
|_|
*/

/*--- you need to install wiennaRNA


options validvarname=v7;
options set=RHOME "C:\Progra~1\R\R-4.5.2\bin\r";
proc r;
submit;
# R program to estimate melting temperatures of RNA sequences
# Clear workspace
rm(list = ls())

# RNA folding function using ViennaRNA
fold_rna <- function(sequence, temp_celsius, rnafold_path = "C:/Program Files (x86)/ViennaRNA Package/RNAfold.exe") {
  input_file <- tempfile(fileext = ".seq")
  writeLines(sequence, input_file)
  cmd <- sprintf('"%s" -T %.1f --noPS -i "%s"', rnafold_path, temp_celsius, input_file)
  output <- system(cmd, intern = TRUE)
  unlink(input_file)

  # Extract structure and count base pairs
  for (line in output) {
    if (grepl("^[\\.\\(\\)]+", trimws(line))) {
      structure <- strsplit(trimws(line), "\\s+")[[1]][1]
      bp_count <- nchar(gsub("[^()]", "", structure))
      return(list(structure = structure, bp_count = bp_count))
    }
  }
  return(list(structure = "", bp_count = 0))
}

# Function to estimate melting temperature
estimate_Tm <- function(sequence, sequence_name) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("Estimating Melting Temperature (Tm) for: %s\n", sequence_name))
  cat(sprintf("Sequence: %s\n", sequence))
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  # Calculate GC content
  seq_upper <- toupper(sequence)
  gc_count <- nchar(gsub("[^GC]", "", seq_upper))
  au_count <- nchar(gsub("[^AU]", "", seq_upper))
  total <- nchar(sequence)
  gc_percent <- (gc_count / total) * 100
  au_percent <- (au_count / total) * 100

  cat(sprintf("Sequence Length: %d nucleotides\n", total))
  cat(sprintf("GC Content: %.1f%% (%d G/C pairs)\n", gc_percent, gc_count))
  cat(sprintf("AU Content: %.1f%% (%d A/U pairs)\n", au_percent, au_count))
  cat("\n")

  # Method 1: Experimental temperature scan (most accurate)
  cat("METHOD 1: Experimental Temperature Scan\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  # Determine temperature range based on GC content
  if (gc_percent > 70) {
    temp_range <- seq(20, 100, by = 5)
  } else if (gc_percent < 30) {
    temp_range <- seq(0, 60, by = 5)
  } else {
    temp_range <- seq(10, 90, by = 5)
  }

  results <- data.frame(
    Temperature = temp_range,
    Structure = character(length(temp_range)),
    BP_Count = integer(length(temp_range)),
    stringsAsFactors = FALSE
  )

  max_bp <- 0
  for (i in seq_along(temp_range)) {
    temp <- temp_range[i]
    result <- fold_rna(sequence, temp)
    results$Structure[i] <- result$structure
    results$BP_Count[i] <- result$bp_count
    if (result$bp_count > max_bp) max_bp <- result$bp_count
  }

  # Calculate fraction folded (normalized to max base pairs)
  results$Fraction_Folded <- results$BP_Count / max_bp

  # Find Tm (temperature where fraction folded = 0.5)
  tm_estimated <- NA
  for (i in 1:(nrow(results)-1)) {
    if (results$Fraction_Folded[i] >= 0.5 && results$Fraction_Folded[i+1] <= 0.5) {
      # Linear interpolation
      t1 <- results$Temperature[i]
      t2 <- results$Temperature[i+1]
      f1 <- results$Fraction_Folded[i]
      f2 <- results$Fraction_Folded[i+1]
      tm_estimated <- t1 + (0.5 - f1) * (t2 - t1) / (f2 - f1)
      break
    }
  }

  # Display results table
  cat(sprintf("%-10s | %-30s | %10s | %12s\n",
              "Temp (°C)", "Structure", "Base Pairs", "Fraction Folded"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in seq_along(temp_range)) {
    cat(sprintf("%3d°C       | %-30s | %10d | %11.2f\n",
                results$Temperature[i],
                results$Structure[i],
                results$BP_Count[i],
                results$Fraction_Folded[i]))
  }

  cat("\n")

  # Method 2: Thermodynamic formula (nearest-neighbor)
  cat("METHOD 2: Thermodynamic Formula (Nearest-Neighbor Model)\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  # Calculate melting temperature using empirical formula
  # Based on: Tm = 81.5 + 41*(GC%) - 675/Length (for DNA/RNA in 1M NaCl)
  tm_formula <- 81.5 + (41 * (gc_percent/100)) - (675 / total)

  # Adjust for RNA (slightly more stable than DNA)
  tm_rna_adjusted <- tm_formula + 2.5

  cat(sprintf("Formula-based Tm (DNA): %.1f°C\n", tm_formula))
  cat(sprintf("Formula-based Tm (RNA): %.1f°C (adjusted for RNA stability)\n", tm_rna_adjusted))
  cat("\n")

  # Method 3: GC content rule of thumb
  cat("METHOD 3: GC Content Rule of Thumb\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  # Each GC pair adds ~4°C to melting temperature
  # Each AU pair adds ~2°C
  tm_gc_ruler <- (gc_count * 4) + (au_count * 2)
  cat(sprintf("GC rule-based Tm: %.1f°C\n", tm_gc_ruler))
  cat("  (Each GC pair = 4°C, each AU pair = 2°C)\n")
  cat("\n")

  # Method 4: Wallace rule (simplest, for short oligos)
  cat("METHOD 4: Wallace Rule (for sequences < 25nt)\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  # Tm = 2°C × (A+T) + 4°C × (G+C)
  if (total <= 25) {
    tm_wallace <- (2 * au_count) + (4 * gc_count)
    cat(sprintf("Wallace rule Tm: %.1f°C\n", tm_wallace))
    cat("  (Tm = 2°C x (A+T) + 4°C x (G+C))\n")
  } else {
    tm_wallace <- NA
    cat("  (Not applicable for sequences > 25nt)\n")
  }
  cat("\n")

  # Determine final Tm estimate
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FINAL MELTING TEMPERATURE ESTIMATE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  # Use experimental scan if successful, otherwise use formulas
  if (!is.na(tm_estimated)) {
    final_tm <- tm_estimated
    cat(sprintf("\nExperimental Tm (50%% melted): %.1f°C\n", final_tm))
  } else {
    # If no clear 50% point, estimate from formula
    final_tm <- tm_rna_adjusted
    cat(sprintf("\nPredicted Tm (from formula): %.1f°C\n", final_tm))
  }

  # Additional interpretation
  cat("\n")
  cat("INTERPRETATION:\n")
  if (final_tm < 40) {
    cat("  * This RNA melts at low temperatures (body temperature may melt it)\n")
    cat("  * Likely to be unfolded at physiological temperature (37°C)\n")
  } else if (final_tm < 60) {
    cat("  * This RNA has moderate thermal stability\n")
    cat("  * Partially stable at body temperature\n")
  } else {
    cat("  * This RNA is thermally stable (requires high heat to melt)\n")
    cat("  * Remains structured at physiological temperature\n")
  }

  cat("\n")

  # Return results
  return(list(
    sequence = sequence,
    name = sequence_name,
    gc_percent = gc_percent,
    tm_experimental = tm_estimated,
    tm_formula = tm_rna_adjusted,
    tm_gc_ruler = tm_gc_ruler,
    tm_final = final_tm,
    melting_curve = results
  ))
}

# Define the two sequences
sequences <- list(
  "AU-rich (20% GC)" = "AAAAUUUUAAAAUUUUAAAAUUUU",
  "GC-rich (80% GC)" = "GGGGCCCCGGGGCCCCGGGGCCCC"
)

# Store all results
all_tm_results <- list()

# Analyze each sequence
for (seq_name in names(sequences)) {
  result <- estimate_Tm(sequences[[seq_name]], seq_name)
  all_tm_results[[seq_name]] <- result
}

# Create comparative summary
cat("\n\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("COMPARATIVE SUMMARY: Melting Temperatures\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

comparison <- data.frame(
  Sequence = names(sequences),
  GC_Content = sapply(all_tm_results, function(x) sprintf("%.1f%%", x$gc_percent)),
  Tm_Experimental = sapply(all_tm_results, function(x)
    ifelse(is.na(x$tm_experimental), "N/A", sprintf("%.1f°C", x$tm_experimental))),
  Tm_Formula = sapply(all_tm_results, function(x) sprintf("%.1f°C", x$tm_formula)),
  Final_Tm = sapply(all_tm_results, function(x) sprintf("%.1f°C", x$tm_final)),
  stringsAsFactors = FALSE
)

print(comparison)

cat("\n")
cat("KEY FINDINGS:\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

au_tm <- all_tm_results[["AU-rich (20% GC)"]]$tm_final
gc_tm <- all_tm_results[["GC-rich (80% GC)"]]$tm_final

cat(sprintf("\n1. AU-rich RNA melts at: %.1f°C\n", au_tm))
cat(sprintf("   GC-rich RNA melts at: %.1f°C\n", gc_tm))
cat(sprintf("   Difference: %.1f°C\n", gc_tm - au_tm))
cat("\n2. The GC-rich sequence requires MUCH higher temperature to melt\n")
cat("   because each G-C pair has 3 hydrogen bonds (vs 2 in A-U pairs)\n")
cat("\n3. Physiological temperature (37°C) will:\n")
cat(sprintf("   * %s the AU-rich RNA (melted)\n", ifelse(au_tm < 37, "MELT", "NOT melt")))
cat(sprintf("   * %s the GC-rich RNA (still structured)\n", ifelse(gc_tm < 37, "melt", "NOT melt")))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("BIOLOGICAL SIGNIFICANCE:\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n")
cat("* Thermophilic bacteria use GC-rich RNA/DNA to maintain structure at high temperatures\n")
cat("* RNA thermometers (AU-rich) melt at body temperature to control gene expression\n")
cat("* PCR primers are designed with higher GC content for thermal stability\n")
cat("* RNA structure predictions must account for temperature effects\n")
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# Create a melting curve plot
plot_melting_curves <- function(all_results) {
  # Find common temperature range
  all_temps <- c()
  all_fractions <- list()

  for (seq_name in names(all_results)) {
    if (!is.null(all_results[[seq_name]]$melting_curve)) {
      curve <- all_results[[seq_name]]$melting_curve
      all_temps <- curve$Temperature
      all_fractions[[seq_name]] <- curve$Fraction_Folded
    }
  }

  if (length(all_fractions) > 0) {
    # Create plot
    plot(all_temps, all_fractions[[1]],
         type = "l",
         col = "red",
         lwd = 3,
         xlim = range(all_temps),
         ylim = c(0, 1),
         xlab = "Temperature (°C)",
         ylab = "Fraction Folded (Structured)",
         main = "RNA Melting Curves Comparison",
         sub = "Higher GC content = Higher melting temperature")

    # Add other sequences
    colors <- c("red", "blue")
    for (i in 2:length(all_fractions)) {
      lines(all_temps, all_fractions[[i]], col = colors[i], lwd = 3)
    }

    # Add horizontal line at 50% folded
    abline(h = 0.5, lty = 2, col = "gray", lwd = 2)

    # Add legend
    legend("topright",
           legend = names(all_fractions),
           col = colors[1:length(all_fractions)],
           lwd = 3,
           bg = "white")

    # Add text annotation
    text(mean(range(all_temps)), 0.55, "Melting Temperature (Tm)", cex = 0.9)
    text(mean(range(all_temps)), 0.45, "50% of structure remains", cex = 0.8)
  }
}

# Create the plot
plot_melting_curves(all_tm_results)
endsubmit;
run;

/*____                      _               _
|___ /   _ __    ___  _   _| |_ _ __  _   _| |_
  |_ \  | `__|  / _ \| | | | __| `_ \| | | | __|
 ___) | | |    | (_) | |_| | |_| |_) | |_| | |_
|____/  |_|     \___/ \__,_|\__| .__/ \__,_|\__|
                               |_|
*/

/**************************************************************************************************************************/
/*   Altair SLC                                                                                                           */
/*                                                                                                                        */
/*   ======================================================================                                               */
/*  Estimating Melting Temperature (Tm) for: AU-rich (20% GC)                                                             */
/*  Sequence: AAAAUUUUAAAAUUUUAAAAUUUU                                                                                    */
/*  ======================================================================                                                */
/*  Sequence Length: 24 nucleotides                                                                                       */
/*  GC Content: 0.0% (0 G/C pairs)                                                                                        */
/*  AU Content: 100.0% (24 A/U pairs)                                                                                     */
/*  METHOD 1: Experimental Temperature Scan                                                                               */
/*  ----------------------------------------------------------------------                                                */
/*  Temp (°C) | Structure                      | Base Pairs | Fraction Folded                                             */
/*  ----------------------------------------------------------------------                                                */
/*    0°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*    5°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   10°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   15°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   20°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   25°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   30°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   35°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   40°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   45°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   50°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   55°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   60°C       | ........................       |          0 |        0.00                                               */
/*  METHOD 2: Thermodynamic Formula (Nearest-Neighbor Model)                                                              */
/*  ----------------------------------------------------------------------                                                */
/*  Formula-based Tm (DNA): 53.4°C                                                                                        */
/*  Formula-based Tm (RNA): 55.9°C (adjusted for RNA stability)                                                           */
/*  METHOD 3: GC Content Rule of Thumb                                                                                    */
/*  ----------------------------------------------------------------------                                                */
/*  GC rule-based Tm: 48.0°C                                                                                              */
/*    (Each GC pair = 4°C, each AU pair = 2°C)                                                                            */
/*  METHOD 4: Wallace Rule (for sequences < 25nt)                                                                         */
/*  ----------------------------------------------------------------------                                                */
/*  Wallace rule Tm: 48.0°C                                                                                               */
/*    (Tm = 2°C x (A+T) + 4°C x (G+C))                                                                                    */
/*  ======================================================================                                                */
/*  FINAL MELTING TEMPERATURE ESTIMATE                                                                                    */
/*  ======================================================================                                                */
/*  Experimental Tm (50% melted): 57.5°C                                                                                  */
/*  INTERPRETATION:                                                                                                       */
/*    * This RNA has moderate thermal stability                                                                           */
/*    * Partially stable at body temperature                                                                              */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*   ======================================================================                                               */
/*  Estimating Melting Temperature (Tm) for: GC-rich (80% GC)                                                             */
/*  Sequence: GGGGCCCCGGGGCCCCGGGGCCCC                                                                                    */
/*  ======================================================================                                                */
/*  Sequence Length: 24 nucleotides                                                                                       */
/*  GC Content: 100.0% (24 G/C pairs)                                                                                     */
/*  AU Content: 0.0% (0 A/U pairs)                                                                                        */
/*  METHOD 1: Experimental Temperature Scan                                                                               */
/*  ----------------------------------------------------------------------                                                */
/*  Temp (°C) | Structure                      | Base Pairs | Fraction Folded                                             */
/*  ----------------------------------------------------------------------                                                */
/*   20°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   25°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   30°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   35°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   40°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   45°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   50°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   55°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   60°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   65°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   70°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   75°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   80°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   85°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   90°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*   95°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*  100°C       | ((((((((((....))))))))))       |         20 |        1.00                                               */
/*  METHOD 2: Thermodynamic Formula (Nearest-Neighbor Model)                                                              */
/*  ----------------------------------------------------------------------                                                */
/*  Formula-based Tm (DNA): 94.4°C                                                                                        */
/*  Formula-based Tm (RNA): 96.9°C (adjusted for RNA stability)                                                           */
/*  METHOD 3: GC Content Rule of Thumb                                                                                    */
/*  ----------------------------------------------------------------------                                                */
/*  GC rule-based Tm: 96.0°C                                                                                              */
/*    (Each GC pair = 4°C, each AU pair = 2°C)                                                                            */
/*  METHOD 4: Wallace Rule (for sequences < 25nt)                                                                         */
/*  ----------------------------------------------------------------------                                                */
/*  Wallace rule Tm: 96.0°C                                                                                               */
/*    (Tm = 2°C x (A+T) + 4°C x (G+C))                                                                                    */
/*  ======================================================================                                                */
/*  FINAL MELTING TEMPERATURE ESTIMATE                                                                                    */
/*  ======================================================================                                                */
/*  Predicted Tm (from formula): 96.9°C                                                                                   */
/*  INTERPRETATION:                                                                                                       */
/*    * This RNA is thermally stable (requires high heat to melt)                                                         */
/*    * Remains structured at physiological temperature                                                                   */
/*  ================================================================================                                      */
/*  COMPARATIVE SUMMARY: Melting Temperatures                                                                             */
/*  ================================================================================                                      */
/*                           Sequence GC_Content Tm_Experimental Tm_Formula Final_Tm                                      */
/*  AU-rich (20% GC) AU-rich (20% GC)       0.0%         57.5Ã‚°C    55.9Ã‚°C  57.5Ã‚°C                                   */
/*  GC-rich (80% GC) GC-rich (80% GC)     100.0%             N/A    96.9Ã‚°C  96.9Ã‚°C                                    */
/*  KEY FINDINGS:                                                                                                         */
/*  --------------------------------------------------------------------------------                                      */
/*  1. AU-rich RNA melts at: 57.5°C                                                                                       */
/*     GC-rich RNA melts at: 96.9°C                                                                                       */
/*     Difference: 39.4°C                                                                                                 */
/*  2. The GC-rich sequence requires MUCH higher temperature to melt                                                      */
/*     because each G-C pair has 3 hydrogen bonds (vs 2 in A-U pairs)                                                     */
/*  3. Physiological temperature (37°C) will:                                                                             */
/*     * NOT melt the AU-rich RNA (melted)                                                                                */
/*     * NOT melt the GC-rich RNA (still structured)                                                                      */
/*  ================================================================================                                      */
/*  BIOLOGICAL SIGNIFICANCE:                                                                                              */
/*  ================================================================================                                      */
/*  * Thermophilic bacteria use GC-rich RNA/DNA to maintain structure at high temperatures                                */
/*  * RNA thermometers (AU-rich) melt at body temperature to control gene expression                                      */
/*  * PCR primers are designed with higher GC content for thermal stability                                               */
/*  * RNA structure predictions must account for temperature effects                                                      */
/*  ================================================================================                                      */
/**************************************************************************************************************************/

/*  _     _
| || |   | | ___   __ _
| || |_  | |/ _ \ / _` |
|__   _| | | (_) | (_| |
   |_|   |_|\___/ \__, |
                  |___/
*/

1                                          Altair SLC         12:39 Tuesday, April 28, 2026

NOTE: Copyright 2002-2025 World Programming, an Altair Company
NOTE: Altair SLC 2026 (05.26.01.00.000758)
      Licensed to Roger DeAngelis
NOTE: This session is executing on the X64_WIN11PRO platform and is running in 64 bit mode

NOTE: AUTOEXEC processing beginning; file is C:\wpsoto\autoexec.sas
NOTE: AUTOEXEC source line
1       +  ï»¿ods _all_ close;
           ^
NOTE: AUTOEXEC processing completed

1
2         options validvarname=v7;
3         options set=RHOME "C:\Progra~1\R\R-4.5.2\bin\r";
4         proc r;
5         submit;
6         # R program to estimate melting temperatures of RNA sequences
7         # Clear workspace
8         rm(list = ls())
9
10        # RNA folding function using ViennaRNA
11        fold_rna <- function(sequence, temp_celsius, rnafold_path = "C:/Program Files (x86)/ViennaRNA Package/RNAfold.exe") {
12          input_file <- tempfile(fileext = ".seq")
13          writeLines(sequence, input_file)
14          cmd <- sprintf('"%s" -T %.1f --noPS -i "%s"', rnafold_path, temp_celsius, input_file)
15          output <- system(cmd, intern = TRUE)
16          unlink(input_file)
17
18          # Extract structure and count base pairs
19          for (line in output) {
20            if (grepl("^[\\.\\(\\)]+", trimws(line))) {
21              structure <- strsplit(trimws(line), "\\s+")[[1]][1]
22              bp_count <- nchar(gsub("[^()]", "", structure))
23              return(list(structure = structure, bp_count = bp_count))
24            }
25          }
26          return(list(structure = "", bp_count = 0))
27        }
28
29        # Function to estimate melting temperature
30        estimate_Tm <- function(sequence, sequence_name) {
31          cat("\n", paste(rep("=", 70), collapse = ""), "\n")
32          cat(sprintf("Estimating Melting Temperature (Tm) for: %s\n", sequence_name))
33          cat(sprintf("Sequence: %s\n", sequence))
34          cat(paste(rep("=", 70), collapse = ""), "\n\n")
35
36          # Calculate GC content
37          seq_upper <- toupper(sequence)
38          gc_count <- nchar(gsub("[^GC]", "", seq_upper))
39          au_count <- nchar(gsub("[^AU]", "", seq_upper))
40          total <- nchar(sequence)
41          gc_percent <- (gc_count / total) * 100
42          au_percent <- (au_count / total) * 100
43
44          cat(sprintf("Sequence Length: %d nucleotides\n", total))
45          cat(sprintf("GC Content: %.1f%% (%d G/C pairs)\n", gc_percent, gc_count))
46          cat(sprintf("AU Content: %.1f%% (%d A/U pairs)\n", au_percent, au_count))
47          cat("\n")
48
49          # Method 1: Experimental temperature scan (most accurate)
50          cat("METHOD 1: Experimental Temperature Scan\n")
51          cat(paste(rep("-", 70), collapse = ""), "\n")
52
53          # Determine temperature range based on GC content
54          if (gc_percent > 70) {
55            temp_range <- seq(20, 100, by = 5)
56          } else if (gc_percent < 30) {
57            temp_range <- seq(0, 60, by = 5)
58          } else {
59            temp_range <- seq(10, 90, by = 5)
60          }
61
62          results <- data.frame(
63            Temperature = temp_range,
64            Structure = character(length(temp_range)),
65            BP_Count = integer(length(temp_range)),
66            stringsAsFactors = FALSE
67          )
68
69          max_bp <- 0
70          for (i in seq_along(temp_range)) {
71            temp <- temp_range[i]
72            result <- fold_rna(sequence, temp)
73            results$Structure[i] <- result$structure
74            results$BP_Count[i] <- result$bp_count
75            if (result$bp_count > max_bp) max_bp <- result$bp_count
76          }
77
78          # Calculate fraction folded (normalized to max base pairs)
79          results$Fraction_Folded <- results$BP_Count / max_bp
80
81          # Find Tm (temperature where fraction folded = 0.5)
82          tm_estimated <- NA
83          for (i in 1:(nrow(results)-1)) {
84            if (results$Fraction_Folded[i] >= 0.5 && results$Fraction_Folded[i+1] <= 0.5) {
85              # Linear interpolation
86              t1 <- results$Temperature[i]
87              t2 <- results$Temperature[i+1]
88              f1 <- results$Fraction_Folded[i]
89              f2 <- results$Fraction_Folded[i+1]
90              tm_estimated <- t1 + (0.5 - f1) * (t2 - t1) / (f2 - f1)
91              break
92            }
93          }
94
95          # Display results table
96          cat(sprintf("%-10s | %-30s | %10s | %12s\n",
97                      "Temp (Â°C)", "Structure", "Base Pairs", "Fraction Folded"))
98          cat(paste(rep("-", 70), collapse = ""), "\n")
99
100         for (i in seq_along(temp_range)) {
101           cat(sprintf("%3dÂ°C       | %-30s | %10d | %11.2f\n",
102                       results$Temperature[i],
103                       results$Structure[i],
104                       results$BP_Count[i],
105                       results$Fraction_Folded[i]))
106         }
107
108         cat("\n")
109
110         # Method 2: Thermodynamic formula (nearest-neighbor)
111         cat("METHOD 2: Thermodynamic Formula (Nearest-Neighbor Model)\n")
112         cat(paste(rep("-", 70), collapse = ""), "\n")
113
114         # Calculate melting temperature using empirical formula
115         # Based on: Tm = 81.5 + 41*(GC%) - 675/Length (for DNA/RNA in 1M NaCl)
116         tm_formula <- 81.5 + (41 * (gc_percent/100)) - (675 / total)
117
118         # Adjust for RNA (slightly more stable than DNA)
119         tm_rna_adjusted <- tm_formula + 2.5
120
121         cat(sprintf("Formula-based Tm (DNA): %.1fÂ°C\n", tm_formula))
122         cat(sprintf("Formula-based Tm (RNA): %.1fÂ°C (adjusted for RNA stability)\n", tm_rna_adjusted))
123         cat("\n")
124
125         # Method 3: GC content rule of thumb
126         cat("METHOD 3: GC Content Rule of Thumb\n")
127         cat(paste(rep("-", 70), collapse = ""), "\n")
128
129         # Each GC pair adds ~4Â°C to melting temperature
130         # Each AU pair adds ~2Â°C
131         tm_gc_ruler <- (gc_count * 4) + (au_count * 2)
132         cat(sprintf("GC rule-based Tm: %.1fÂ°C\n", tm_gc_ruler))
133         cat("  (Each GC pair = 4Â°C, each AU pair = 2Â°C)\n")
134         cat("\n")
135
136         # Method 4: Wallace rule (simplest, for short oligos)
137         cat("METHOD 4: Wallace Rule (for sequences < 25nt)\n")
138         cat(paste(rep("-", 70), collapse = ""), "\n")
139
140         # Tm = 2Â°C Ã— (A+T) + 4Â°C Ã— (G+C)
141         if (total <= 25) {
142           tm_wallace <- (2 * au_count) + (4 * gc_count)
143           cat(sprintf("Wallace rule Tm: %.1fÂ°C\n", tm_wallace))
144           cat("  (Tm = 2Â°C x (A+T) + 4Â°C x (G+C))\n")
145         } else {
146           tm_wallace <- NA
147           cat("  (Not applicable for sequences > 25nt)\n")
148         }
149         cat("\n")
150
151         # Determine final Tm estimate
152         cat(paste(rep("=", 70), collapse = ""), "\n")
153         cat("FINAL MELTING TEMPERATURE ESTIMATE\n")
154         cat(paste(rep("=", 70), collapse = ""), "\n")
155
156         # Use experimental scan if successful, otherwise use formulas
157         if (!is.na(tm_estimated)) {
158           final_tm <- tm_estimated
159           cat(sprintf("\nExperimental Tm (50%% melted): %.1fÂ°C\n", final_tm))
160         } else {
161           # If no clear 50% point, estimate from formula
162           final_tm <- tm_rna_adjusted
163           cat(sprintf("\nPredicted Tm (from formula): %.1fÂ°C\n", final_tm))
164         }
165
166         # Additional interpretation
167         cat("\n")
168         cat("INTERPRETATION:\n")
169         if (final_tm < 40) {
170           cat("  * This RNA melts at low temperatures (body temperature may melt it)\n")
171           cat("  * Likely to be unfolded at physiological temperature (37Â°C)\n")
172         } else if (final_tm < 60) {
173           cat("  * This RNA has moderate thermal stability\n")
174           cat("  * Partially stable at body temperature\n")
175         } else {
176           cat("  * This RNA is thermally stable (requires high heat to melt)\n")
177           cat("  * Remains structured at physiological temperature\n")
178         }
179
180         cat("\n")
181
182         # Return results
183         return(list(
184           sequence = sequence,
185           name = sequence_name,
186           gc_percent = gc_percent,
187           tm_experimental = tm_estimated,
188           tm_formula = tm_rna_adjusted,
189           tm_gc_ruler = tm_gc_ruler,
190           tm_final = final_tm,
191           melting_curve = results
192         ))
193       }
194
195       # Define the two sequences
196       sequences <- list(
197         "AU-rich (20% GC)" = "AAAAUUUUAAAAUUUUAAAAUUUU",
198         "GC-rich (80% GC)" = "GGGGCCCCGGGGCCCCGGGGCCCC"
199       )
200
201       # Store all results
202       all_tm_results <- list()
203
204       # Analyze each sequence
205       for (seq_name in names(sequences)) {
206         result <- estimate_Tm(sequences[[seq_name]], seq_name)
207         all_tm_results[[seq_name]] <- result
208       }
209
210       # Create comparative summary
211       cat("\n\n")
212       cat(paste(rep("=", 80), collapse = ""), "\n")
213       cat("COMPARATIVE SUMMARY: Melting Temperatures\n")
214       cat(paste(rep("=", 80), collapse = ""), "\n\n")
215
216       comparison <- data.frame(
217         Sequence = names(sequences),
218         GC_Content = sapply(all_tm_results, function(x) sprintf("%.1f%%", x$gc_percent)),
219         Tm_Experimental = sapply(all_tm_results, function(x)
220           ifelse(is.na(x$tm_experimental), "N/A", sprintf("%.1fÂ°C", x$tm_experimental))),
221         Tm_Formula = sapply(all_tm_results, function(x) sprintf("%.1fÂ°C", x$tm_formula)),
222         Final_Tm = sapply(all_tm_results, function(x) sprintf("%.1fÂ°C", x$tm_final)),
223         stringsAsFactors = FALSE
224       )
225
226       print(comparison)
227
228       cat("\n")
229       cat("KEY FINDINGS:\n")
230       cat(paste(rep("-", 80), collapse = ""), "\n")
231
232       au_tm <- all_tm_results[["AU-rich (20% GC)"]]$tm_final
233       gc_tm <- all_tm_results[["GC-rich (80% GC)"]]$tm_final
234
235       cat(sprintf("\n1. AU-rich RNA melts at: %.1fÂ°C\n", au_tm))
236       cat(sprintf("   GC-rich RNA melts at: %.1fÂ°C\n", gc_tm))
237       cat(sprintf("   Difference: %.1fÂ°C\n", gc_tm - au_tm))
238       cat("\n2. The GC-rich sequence requires MUCH higher temperature to melt\n")
239       cat("   because each G-C pair has 3 hydrogen bonds (vs 2 in A-U pairs)\n")
240       cat("\n3. Physiological temperature (37Â°C) will:\n")
241       cat(sprintf("   * %s the AU-rich RNA (melted)\n", ifelse(au_tm < 37, "MELT", "NOT melt")))
242       cat(sprintf("   * %s the GC-rich RNA (still structured)\n", ifelse(gc_tm < 37, "melt", "NOT melt")))
243
244       cat("\n")
245       cat(paste(rep("=", 80), collapse = ""), "\n")
246       cat("BIOLOGICAL SIGNIFICANCE:\n")
247       cat(paste(rep("=", 80), collapse = ""), "\n")
248       cat("\n")
249       cat("* Thermophilic bacteria use GC-rich RNA/DNA to maintain structure at high temperatures\n")
250       cat("* RNA thermometers (AU-rich) melt at body temperature to control gene expression\n")
251       cat("* PCR primers are designed with higher GC content for thermal stability\n")
252       cat("* RNA structure predictions must account for temperature effects\n")
253       cat("\n")
254       cat(paste(rep("=", 80), collapse = ""), "\n")
255
256       # Create a melting curve plot
257       plot_melting_curves <- function(all_results) {
258         # Find common temperature range
259         all_temps <- c()
260         all_fractions <- list()
261
262         for (seq_name in names(all_results)) {
263           if (!is.null(all_results[[seq_name]]$melting_curve)) {
264             curve <- all_results[[seq_name]]$melting_curve
265             all_temps <- curve$Temperature
266             all_fractions[[seq_name]] <- curve$Fraction_Folded
267           }
268         }
269
270         if (length(all_fractions) > 0) {
271           # Create plot
272           plot(all_temps, all_fractions[[1]],
273                type = "l",
274                col = "red",
275                lwd = 3,
276                xlim = range(all_temps),
277                ylim = c(0, 1),
278                xlab = "Temperature (Â°C)",
279                ylab = "Fraction Folded (Structured)",
280                main = "RNA Melting Curves Comparison",
281                sub = "Higher GC content = Higher melting temperature")
282
283           # Add other sequences
284           colors <- c("red", "blue")
285           for (i in 2:length(all_fractions)) {
286             lines(all_temps, all_fractions[[i]], col = colors[i], lwd = 3)
287           }
288
289           # Add horizontal line at 50% folded
290           abline(h = 0.5, lty = 2, col = "gray", lwd = 2)
291
292           # Add legend
293           legend("topright",
294                  legend = names(all_fractions),
295                  col = colors[1:length(all_fractions)],
296                  lwd = 3,
297                  bg = "white")
298
299           # Add text annotation
300           text(mean(range(all_temps)), 0.55, "Melting Temperature (Tm)", cex = 0.9)
301           text(mean(range(all_temps)), 0.45, "50% of structure remains", cex = 0.8)
302         }
303       }
304
305       # Create the plot
306       plot_melting_curves(all_tm_results)
307       endsubmit;
NOTE: Using R version 4.5.2 (2025-10-31 ucrt) from C:\Program Files\R\R-4.5.2

NOTE: Submitting statements to R:

> # R program to estimate melting temperatures of RNA sequences
> # Clear workspace
> rm(list = ls())
>
> # RNA folding function using ViennaRNA
> fold_rna <- function(sequence, temp_celsius, rnafold_path = "C:/Program Files (x86)/ViennaRNA Package/RNAfold.exe") {
+   input_file <- tempfile(fileext = ".seq")
+   writeLines(sequence, input_file)
+   cmd <- sprintf('"%s" -T %.1f --noPS -i "%s"', rnafold_path, temp_celsius, input_file)
+   output <- system(cmd, intern = TRUE)
+   unlink(input_file)
+
+   # Extract structure and count base pairs
+   for (line in output) {
+     if (grepl("^[\\.\\(\\)]+", trimws(line))) {
+       structure <- strsplit(trimws(line), "\\s+")[[1]][1]
+       bp_count <- nchar(gsub("[^()]", "", structure))
+       return(list(structure = structure, bp_count = bp_count))
+     }
+   }
+   return(list(structure = "", bp_count = 0))
+ }
>
> # Function to estimate melting temperature
> estimate_Tm <- function(sequence, sequence_name) {
+   cat("\n", paste(rep("=", 70), collapse = ""), "\n")
+   cat(sprintf("Estimating Melting Temperature (Tm) for: %s\n", sequence_name))
+   cat(sprintf("Sequence: %s\n", sequence))
+   cat(paste(rep("=", 70), collapse = ""), "\n\n")
+
+   # Calculate GC content
+   seq_upper <- toupper(sequence)
+   gc_count <- nchar(gsub("[^GC]", "", seq_upper))
+   au_count <- nchar(gsub("[^AU]", "", seq_upper))
+   total <- nchar(sequence)
+   gc_percent <- (gc_count / total) * 100
+   au_percent <- (au_count / total) * 100
+
+   cat(sprintf("Sequence Length: %d nucleotides\n", total))
+   cat(sprintf("GC Content: %.1f%% (%d G/C pairs)\n", gc_percent, gc_count))
+   cat(sprintf("AU Content: %.1f%% (%d A/U pairs)\n", au_percent, au_count))
+   cat("\n")
+
+   # Method 1: Experimental temperature scan (most accurate)
+   cat("METHOD 1: Experimental Temperature Scan\n")
+   cat(paste(rep("-", 70), collapse = ""), "\n")
+
+   # Determine temperature range based on GC content
+   if (gc_percent > 70) {
+     temp_range <- seq(20, 100, by = 5)
+   } else if (gc_percent < 30) {
+     temp_range <- seq(0, 60, by = 5)
+   } else {
+     temp_range <- seq(10, 90, by = 5)
+   }
+
+   results <- data.frame(
+     Temperature = temp_range,
+     Structure = character(length(temp_range)),
+     BP_Count = integer(length(temp_range)),
+     stringsAsFactors = FALSE
+   )
+
+   max_bp <- 0
+   for (i in seq_along(temp_range)) {
+     temp <- temp_range[i]
+     result <- fold_rna(sequence, temp)
+     results$Structure[i] <- result$structure
+     results$BP_Count[i] <- result$bp_count
+     if (result$bp_count > max_bp) max_bp <- result$bp_count
+   }
+
+   # Calculate fraction folded (normalized to max base pairs)
+   results$Fraction_Folded <- results$BP_Count / max_bp
+
+   # Find Tm (temperature where fraction folded = 0.5)
+   tm_estimated <- NA
+   for (i in 1:(nrow(results)-1)) {
+     if (results$Fraction_Folded[i] >= 0.5 && results$Fraction_Folded[i+1] <= 0.5) {
+       # Linear interpolation
+       t1 <- results$Temperature[i]
+       t2 <- results$Temperature[i+1]
+       f1 <- results$Fraction_Folded[i]
+       f2 <- results$Fraction_Folded[i+1]
+       tm_estimated <- t1 + (0.5 - f1) * (t2 - t1) / (f2 - f1)
+       break
+     }
+   }
+
+   # Display results table
+   cat(sprintf("%-10s | %-30s | %10s | %12s\n",
+               "Temp (Â°C)", "Structure", "Base Pairs", "Fraction Folded"))
+   cat(paste(rep("-", 70), collapse = ""), "\n")
+
+   for (i in seq_along(temp_range)) {
+     cat(sprintf("%3dÂ°C       | %-30s | %10d | %11.2f\n",
+                 results$Temperature[i],
+                 results$Structure[i],
+                 results$BP_Count[i],
+                 results$Fraction_Folded[i]))
+   }
+
+   cat("\n")
+
+   # Method 2: Thermodynamic formula (nearest-neighbor)
+   cat("METHOD 2: Thermodynamic Formula (Nearest-Neighbor Model)\n")
+   cat(paste(rep("-", 70), collapse = ""), "\n")
+
+   # Calculate melting temperature using empirical formula
+   # Based on: Tm = 81.5 + 41*(GC%) - 675/Length (for DNA/RNA in 1M NaCl)
+   tm_formula <- 81.5 + (41 * (gc_percent/100)) - (675 / total)
+
+   # Adjust for RNA (slightly more stable than DNA)
+   tm_rna_adjusted <- tm_formula + 2.5
+
+   cat(sprintf("Formula-based Tm (DNA): %.1fÂ°C\n", tm_formula))
+   cat(sprintf("Formula-based Tm (RNA): %.1fÂ°C (adjusted for RNA stability)\n", tm_rna_adjusted))
+   cat("\n")
+
+   # Method 3: GC content rule of thumb
+   cat("METHOD 3: GC Content Rule of Thumb\n")
+   cat(paste(rep("-", 70), collapse = ""), "\n")
+
+   # Each GC pair adds ~4Â°C to melting temperature
+   # Each AU pair adds ~2Â°C
+   tm_gc_ruler <- (gc_count * 4) + (au_count * 2)
+   cat(sprintf("GC rule-based Tm: %.1fÂ°C\n", tm_gc_ruler))
+   cat("  (Each GC pair = 4Â°C, each AU pair = 2Â°C)\n")
+   cat("\n")
+
+   # Method 4: Wallace rule (simplest, for short oligos)
+   cat("METHOD 4: Wallace Rule (for sequences < 25nt)\n")
+   cat(paste(rep("-", 70), collapse = ""), "\n")
+
+   # Tm = 2Â°C Ã— (A+T) + 4Â°C Ã— (G+C)
+   if (total <= 25) {
+     tm_wallace <- (2 * au_count) + (4 * gc_count)
+     cat(sprintf("Wallace rule Tm: %.1fÂ°C\n", tm_wallace))
+     cat("  (Tm = 2Â°C x (A+T) + 4Â°C x (G+C))\n")
+   } else {
+     tm_wallace <- NA
+     cat("  (Not applicable for sequences > 25nt)\n")
+   }
+   cat("\n")
+
+   # Determine final Tm estimate
+   cat(paste(rep("=", 70), collapse = ""), "\n")

2                                                                                                                         Altair SLC

+   cat("FINAL MELTING TEMPERATURE ESTIMATE\n")
+   cat(paste(rep("=", 70), collapse = ""), "\n")
+
+   # Use experimental scan if successful, otherwise use formulas
+   if (!is.na(tm_estimated)) {
+     final_tm <- tm_estimated
+     cat(sprintf("\nExperimental Tm (50%% melted): %.1fÂ°C\n", final_tm))
+   } else {
+     # If no clear 50% point, estimate from formula
+     final_tm <- tm_rna_adjusted
+     cat(sprintf("\nPredicted Tm (from formula): %.1fÂ°C\n", final_tm))
+   }
+
+   # Additional interpretation
+   cat("\n")
+   cat("INTERPRETATION:\n")
+   if (final_tm < 40) {
+     cat("  * This RNA melts at low temperatures (body temperature may melt it)\n")
+     cat("  * Likely to be unfolded at physiological temperature (37Â°C)\n")
+   } else if (final_tm < 60) {
+     cat("  * This RNA has moderate thermal stability\n")
+     cat("  * Partially stable at body temperature\n")
+   } else {
+     cat("  * This RNA is thermally stable (requires high heat to melt)\n")
+     cat("  * Remains structured at physiological temperature\n")
+   }
+
+   cat("\n")
+
+   # Return results
+   return(list(
+     sequence = sequence,
+     name = sequence_name,
+     gc_percent = gc_percent,
+     tm_experimental = tm_estimated,
+     tm_formula = tm_rna_adjusted,
+     tm_gc_ruler = tm_gc_ruler,
+     tm_final = final_tm,
+     melting_curve = results
+   ))
+ }
>
> # Define the two sequences
> sequences <- list(
+   "AU-rich (20% GC)" = "AAAAUUUUAAAAUUUUAAAAUUUU",
+   "GC-rich (80% GC)" = "GGGGCCCCGGGGCCCCGGGGCCCC"
+ )
>
> # Store all results
> all_tm_results <- list()
>
> # Analyze each sequence
> for (seq_name in names(sequences)) {
+   result <- estimate_Tm(sequences[[seq_name]], seq_name)
+   all_tm_results[[seq_name]] <- result
+ }
>
> # Create comparative summary
> cat("\n\n")
> cat(paste(rep("=", 80), collapse = ""), "\n")
> cat("COMPARATIVE SUMMARY: Melting Temperatures\n")
> cat(paste(rep("=", 80), collapse = ""), "\n\n")
>
> comparison <- data.frame(
+   Sequence = names(sequences),
+   GC_Content = sapply(all_tm_results, function(x) sprintf("%.1f%%", x$gc_percent)),
+   Tm_Experimental = sapply(all_tm_results, function(x)
+     ifelse(is.na(x$tm_experimental), "N/A", sprintf("%.1fÂ°C", x$tm_experimental))),
+   Tm_Formula = sapply(all_tm_results, function(x) sprintf("%.1fÂ°C", x$tm_formula)),
+   Final_Tm = sapply(all_tm_results, function(x) sprintf("%.1fÂ°C", x$tm_final)),
+   stringsAsFactors = FALSE
+ )
>
> print(comparison)
>
> cat("\n")
> cat("KEY FINDINGS:\n")
> cat(paste(rep("-", 80), collapse = ""), "\n")
>
> au_tm <- all_tm_results[["AU-rich (20% GC)"]]$tm_final
> gc_tm <- all_tm_results[["GC-rich (80% GC)"]]$tm_final
>
> cat(sprintf("\n1. AU-rich RNA melts at: %.1fÂ°C\n", au_tm))
> cat(sprintf("   GC-rich RNA melts at: %.1fÂ°C\n", gc_tm))
> cat(sprintf("   Difference: %.1fÂ°C\n", gc_tm - au_tm))
> cat("\n2. The GC-rich sequence requires MUCH higher temperature to melt\n")
> cat("   because each G-C pair has 3 hydrogen bonds (vs 2 in A-U pairs)\n")
> cat("\n3. Physiological temperature (37Â°C) will:\n")
> cat(sprintf("   * %s the AU-rich RNA (melted)\n", ifelse(au_tm < 37, "MELT", "NOT melt")))
> cat(sprintf("   * %s the GC-rich RNA (still structured)\n", ifelse(gc_tm < 37, "melt", "NOT melt")))
>
> cat("\n")
> cat(paste(rep("=", 80), collapse = ""), "\n")
> cat("BIOLOGICAL SIGNIFICANCE:\n")
> cat(paste(rep("=", 80), collapse = ""), "\n")
> cat("\n")
> cat("* Thermophilic bacteria use GC-rich RNA/DNA to maintain structure at high temperatures\n")
> cat("* RNA thermometers (AU-rich) melt at body temperature to control gene expression\n")
> cat("* PCR primers are designed with higher GC content for thermal stability\n")
> cat("* RNA structure predictions must account for temperature effects\n")
> cat("\n")
> cat(paste(rep("=", 80), collapse = ""), "\n")
>
> # Create a melting curve plot
> plot_melting_curves <- function(all_results) {
+   # Find common temperature range
+   all_temps <- c()
+   all_fractions <- list()
+
+   for (seq_name in names(all_results)) {
+     if (!is.null(all_results[[seq_name]]$melting_curve)) {
+       curve <- all_results[[seq_name]]$melting_curve
+       all_temps <- curve$Temperature
+       all_fractions[[seq_name]] <- curve$Fraction_Folded
+     }
+   }
+
+   if (length(all_fractions) > 0) {
+     # Create plot
+     plot(all_temps, all_fractions[[1]],
+          type = "l",
+          col = "red",
+          lwd = 3,
+          xlim = range(all_temps),
+          ylim = c(0, 1),
+          xlab = "Temperature (Â°C)",
+          ylab = "Fraction Folded (Structured)",
+          main = "RNA Melting Curves Comparison",
+          sub = "Higher GC content = Higher melting temperature")
+
+     # Add other sequences
+     colors <- c("red", "blue")
+     for (i in 2:length(all_fractions)) {
+       lines(all_temps, all_fractions[[i]], col = colors[i], lwd = 3)
+     }
+
+     # Add horizontal line at 50% folded
+     abline(h = 0.5, lty = 2, col = "gray", lwd = 2)
+
+     # Add legend
+     legend("topright",
+            legend = names(all_fractions),
+            col = colors[1:length(all_fractions)],
+            lwd = 3,
+            bg = "white")
+
+     # Add text annotation
+     text(mean(range(all_temps)), 0.55, "Melting Temperature (Tm)", cex = 0.9)
+     text(mean(range(all_temps)), 0.45, "50% of structure remains", cex = 0.8)
+   }
+ }
>
> # Create the plot
> plot_melting_curves(all_tm_results)
Error in xy.coords(x, y, xlabel, ylabel, log) :
  'x' and 'y' lengths differ
Calls: plot_melting_curves -> plot -> plot.default -> xy.coords

NOTE: Processing of R statements complete

308       run;
NOTE: Procedure r step took :
      real time : 7.117
      cpu time  : 0.015


309
ERROR: Error printed on page 1

NOTE: Submitted statements took :
      real time : 7.188
      cpu time  : 0.078

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/

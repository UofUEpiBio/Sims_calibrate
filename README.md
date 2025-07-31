

---

# Sims\_calibrate

This repository contains scripts for parameter generation, model loading, calibration, prediction, and analysis for epidemic simulation studies.

## File Descriptions

### `00-params.R`

Generates the simulation parameters used throughout the analysis.

### `01a-bilstm.R`

Loads the pre-trained BiLSTM model for epidemic curve prediction.

### `01b-abc.R`

Applies Approximate Bayesian Computation (ABC) for calibration.

### `02-abc-bilstm-prediction.R`

Combines ABC calibration with BiLSTM predictions and saves the predicted results.

### `03-epicurves-stats.R`

Uses the predicted parameters to generate epidemic curves for visualization and analysis.

### `04-parameter-stats.R`

Analyzes parameter predictions to compute bias and create boxplots for parameter distributions.

## Project Notes

* **R version:** Ensure you are running a compatible R version with necessary packages installed before running the scripts.
* **Order of execution:** Scripts are intended to be run in the order listed above for consistent results.
* **Output:** Prediction results, epidemic curves, and statistical plots are saved to the designated output folders as defined in the scripts.

## Getting Started

1. Generate parameters:

   ```bash
   Rscript 00-params.R
   ```
2. Load model and run calibration:

   ```bash
   Rscript 01a-bilstm.R
   Rscript 01b-abc.R
   ```
3. Combine calibration and predictions:

   ```bash
   Rscript 02-abc-bilstm-prediction.R
   ```
4. Generate epicurves:

   ```bash
   Rscript 03-epicurves-stats.R
   ```
5. Analyze parameter statistics:

   ```bash
   Rscript 04-parameter-stats.R
   ```

---

If you want, I can now give you a **final copy-paste–ready README.md file** so it renders perfectly on GitHub without you having to adjust anything.
Do you want me to do that now?


# 📊 Sims\_calibrate

This repository contains R scripts for **parameter generation**, **model loading**, **calibration**, **prediction**, and **analysis** for epidemic simulation studies.
Whether you want to run simulations 🧮, calibrate predictions 🔍, or visualize epidemic curves 📈 — this project has you covered.

---

## 📂 File Descriptions

### `00-params.R`

📝 Generates the **simulation parameters** used throughout the analysis.

### `01a-bilstm.R`

🤖 Loads the **pre-trained BiLSTM model** for epidemic curve prediction.

### `01b-abc.R`

🎯 Applies **Approximate Bayesian Computation (ABC)** for calibration.

### `02-abc-bilstm-prediction.R`

🔗 Combines **ABC calibration** with **BiLSTM predictions** and saves the results.

### `03-epicurves-stats.R`

📈 Generates **epidemic curves** from the predicted parameters.

### `04-parameter-stats.R`

📊 Analyzes predicted parameters to find **bias** and create **boxplots**.

---

## 🛠️ Project Notes

* **R version:** Make sure you’re using a compatible R version 📦 with all required packages installed.
* **Execution order:** Run scripts in the order listed above for consistent results.
* **Outputs:** Predictions, epidemic curves, and statistical plots are saved in the designated output folders.

---

## 🚀 Getting Started

1️⃣ **Generate parameters**

```bash
Rscript 00-params.R
```

2️⃣ **Load model & run calibration**

```bash
Rscript 01a-bilstm.R
Rscript 01b-abc.R
```

3️⃣ **Combine calibration & predictions**

```bash
Rscript 02-abc-bilstm-prediction.R
```

4️⃣ **Generate epicurves**

```bash
Rscript 03-epicurves-stats.R
```

5️⃣ **Analyze parameter statistics**

```bash
Rscript 04-parameter-stats.R
```


# Entire R script: load Bi-LSTM + scalers (with fallback for unfitted incidence) via reticulate

library(reticulate)

# 1) Write Python helper (model_utils.py) ---------------------------------------
python_code <- c(
  "# model_utils.py",
  "import torch",
  "import torch.nn as nn",
  "import joblib",
  "import numpy as np",
  "import traceback",
  "",
  "_model = None",
  "_scaler_add = None",
  "_scaler_tgt = None",
  "_scaler_inc = None",
  "_device = torch.device('cpu')",
  "",
  "class BiLSTMModel(nn.Module):",
  "    def __init__(self, input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout):",
  "        super().__init__()",
  "        self.bilstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True, dropout=dropout, bidirectional=True)",
  "        self.fc1 = nn.Linear(2 * hidden_dim + additional_dim, 64)",
  "        self.fc2 = nn.Linear(64, output_dim)",
  "        self.sigmoid = nn.Sigmoid()",
  "        self.softplus = nn.Softplus()",
  "",
  "    def forward(self, x, additional_inputs):",
  "        _, (h_n, _) = self.bilstm(x)",
  "        hid = torch.cat((h_n[-2], h_n[-1]), dim=1)",
  "        combined = torch.cat((hid, additional_inputs), dim=1)",
  "        x = self.fc1(combined)",
  "        out = self.fc2(x)",
  "        return torch.stack([",
  "            self.sigmoid(out[:,0]),",
  "            self.softplus(out[:,1]),",
  "            self.softplus(out[:,2])",
  "        ], dim=1)",
  "",
  "def load_model(model_path, scaler_additional_path, scaler_targets_path, scaler_incidence_path):",
  "    global _model, _scaler_add, _scaler_tgt, _scaler_inc",
  "    try:",
  "        _scaler_add = joblib.load(scaler_additional_path)",
  "        _scaler_tgt = joblib.load(scaler_targets_path)",
  "        _scaler_inc = joblib.load(scaler_incidence_path)",
  "        # instantiate with training hyperparams",
  "        _model = BiLSTMModel(input_dim=1, hidden_dim=160, num_layers=3, additional_dim=2, output_dim=3, dropout=0.5)",
  "        state = torch.load(model_path, map_location=_device)",
  "        _model.load_state_dict(state)",
  "        _model.to(_device).eval()",
  "    except Exception as e:",
  "        print('Error loading model/scalers:', e)",
  "        traceback.print_exc()",
  "        raise",
  "",
  "def predict(seq, additional_pair):",
  "    try:",
  "        x = np.asarray(seq, dtype=np.float32)",
  "        flat = x.reshape(x.shape[0], x.shape[1])",
  "        # fallback: if incidence scaler not fitted, use raw",
  "        try:",
  "            flat_scaled = _scaler_inc.transform(flat)",
  "        except Exception:",
  "            flat_scaled = flat",
  "        x_scaled = flat_scaled.reshape(x.shape[0], x.shape[1], 1)",
  "        add_np = np.array([additional_pair], dtype=np.float32)",
  "        add_scaled = _scaler_add.transform(add_np)",
  "        x_t = torch.tensor(x_scaled, dtype=torch.float32, device=_device)",
  "        add_t = torch.tensor(add_scaled, dtype=torch.float32, device=_device)",
  "        with torch.no_grad():",
  "            out = _model(x_t, add_t).cpu().numpy()",
  "        return _scaler_tgt.inverse_transform(out)[0].tolist()",
  "    except Exception as e:",
  "        print('Prediction error:', e)",
  "        traceback.print_exc()",
  "        raise"
)
writeLines(python_code, "model_utils.py")

# 2) Source helper --------------------------------------------------------------
source_python("model_utils.py")  # provides load_model() and predict()

# 3) Paths & load ---------------------------------------------------------------
base_dir <- "~/Desktop/Sims_calibrate"
model_path            <- normalizePath(file.path(base_dir, "model4_bilstm (2).pt"))
scaler_incidence_path <- normalizePath(file.path(base_dir, "scaler_incidence 1 (1).pkl"))
scaler_additional_path<- normalizePath(file.path(base_dir, "scaler_additional 1 (1).pkl"))
scaler_targets_path   <- normalizePath(file.path(base_dir, "scaler_targets 1 (1).pkl"))

load_model(model_path,
           scaler_additional_path,
           scaler_targets_path,
           scaler_incidence_path)

# 4) R wrapper ------------------------------------------------------------------
predict_with_bilstm <- function(time_series, n, recov) {
  stopifnot(is.numeric(time_series), length(time_series) == 60,
            is.numeric(n), is.numeric(recov))
  arr <- array(time_series, dim = c(1, 60, 1))
  out <- predict(arr, list(n, recov))
  names(out) <- c("ptran","crate","R0")
  out
}

# 5) Test -----------------------------------------------------------------------
set.seed(123)
example_ts <- runif(60, min = 0, max = 100)
print(predict_with_bilstm(example_ts, n = 5000, recov = 0.1))





# Load the epiworldR package
library(epiworldR)
R0=3
# --- Original model
model_sir <- ModelSIRCONN(
  name = "COVID-19",
  prevalence = 0.01,
  n = 5000,
  contact_rate = 10,
  transmission_rate = 0.03,
  recovery_rate = 0.1
)
run(model_sir,ndays = 60)
incidence <- plot_incidence(model_sir,plot=FALSE)
infected_original=incidence[,1]
?plot_incidence

saver <- make_saver("total_hist")
run_multiple(model_sir, ndays = 60, nsims = 1, saver = saver, nthread = 2)
ans <- run_multiple_get_results(model_sir)
df_original <- ans$total_hist
# Extract infected counts
infected_original <- df_original[df_original$state == "Infected", 5]
infected_original <- infected_original[-1]  # remove day 0
#length(infected_original)
#infected_df <- df_original %>%
#  filter(state == "Infected") %>%
#  arrange(sim_num, date)

#crate_true=(r0*recov)/ptran
#3*0.1/0.4
#1.5*0.1/0.03
# Saver and run

recov <- 0.1
n <- 5000
result <- predict_with_bilstm(infected_original, n, recov)
# --- Calibrate
#reticulate::py_install("torch")
#reticulate::py_install("joblib")

recov <- 0.1
n <- 5000
ptran <- result[1]
crate <- result[2]
r0 <- result[3]
crate_true=(r0*recov)/ptran
# --- Calibrated model
calibrated_model <- ModelSIRCONN(
  name = "Calibrated COVID-19",
  prevalence = 0.01,
  n = n,
  contact_rate = crate_true,
  transmission_rate = ptran,
  recovery_rate = recov
)

run_multiple(calibrated_model, ndays = 50, nsims = 1, saver = saver, nthread = 8)
results <- run_multiple_get_results(calibrated_model)
df_calibrated <- results$total_hist

infected_calibrated <- df_calibrated[df_calibrated$state == "Infected", 5]
infected_calibrated <- infected_calibrated[-1]

# --- Plot both curves
library(ggplot2)
library(dplyr)

# Align lengths first
min_len <- min(length(infected_original), length(infected_calibrated))
infected_original <- infected_original[1:min_len]
infected_calibrated <- infected_calibrated[1:min_len]
days <- 1:min_len

# Combine into a single data frame for ggplot
df_plot <- data.frame(
  Day = rep(days, 2),
  Infected = c(infected_original, infected_calibrated),
  Model = rep(c("Original", "Calibrated"), each = min_len)
)

# Plot using ggplot
ggplot(df_plot, aes(x = Day, y = Infected, color = Model, linetype = Model)) +
  geom_line(size = 1.2) +
  labs(
    title = "Infected Over Time: Original vs Calibrated Model",
    x = "Day",
    y = "Number of Infected Individuals"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Original" = "blue", "Calibrated" = "red")) +
  theme(legend.title = element_blank())



####
library(epiworldR)
library(dplyr)
library(ggplot2)

# --- your inputs -------------------------------------
recov <- 0.1
n     <- 5000
ptran <- result[1]
crate <- result[2]
r0    <- result[3]

# compute the “calibrated” contact rate so that R0 = r0_pred
crate_true <- (r0 * recov) / ptran

# --- define two models --------------------------------
model_uncalib <- ModelSIRCONN(
  name              = "Uncalibrated",
  prevalence        = 0.01,
  n                 = 5000,
  contact_rate      = 10,
  transmission_rate = 0.03,
  recovery_rate     = 0.1
)

model_calib <- ModelSIRCONN(
  name              = "Calibrated",
  prevalence        = 0.01,
  n                 = 5000,
  contact_rate      = crate_true,
  transmission_rate = ptran,
  recovery_rate     = recov
)

# --- set up a saver for the active (I) compartment ----
saver <- make_saver("total_hist")

# --- run 20 sims × 60 days in parallel -----------------
run_multiple(model_uncalib, ndays = 60, nsims = 20,
             saver = saver, nthread = 10)
run_multiple(model_calib,  ndays = 60, nsims = 20,
             saver = saver, nthread = 10)

# --- extract & label results --------------------------
res_uncalib <- run_multiple_get_results(model_uncalib)$total_hist %>%
  filter(state == "Infected") %>%
  rename(time = date, rep = sim_num, I = counts) %>%
  mutate(model = "Uncalibrated")

res_calib <- run_multiple_get_results(model_calib)$total_hist %>%
  filter(state == "Infected") %>%
  rename(time = date, rep = sim_num, I = counts) %>%
  mutate(model = "Calibrated")

all_results <- bind_rows(res_uncalib, res_calib)

# --- summarize mean + 95% CI --------------------------
summary_df <- all_results %>%
  group_by(model, time) %>%
  summarize(
    mean_I = mean(I),
    lower  = quantile(I, 0.025),
    upper  = quantile(I, 0.975),
    .groups = "drop"
  )

# --- plot ----------------------------------------------
p <- ggplot(summary_df, aes(x = time, y = mean_I, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    x     = "Day",
    y     = "Active Infected (I)",
    title = "True vs. Calibrated SIRCONN Simulations",
    color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14)

print(p)
View(summary_df)

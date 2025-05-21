library(reticulate)

# 1) Write Python helper (model_utils.py) with fixed incidence scaling ---------------------------------------
python_code <- c(
  "# model_utils.py",
  "import torch",
  "import torch.nn as nn",
  "import joblib",
  "import numpy as np",
  "import traceback",
  "from sklearn.preprocessing import MinMaxScaler",
  "",
  "_model = None",
  "_scaler_add = None",
  "_scaler_tgt = None",
  "_scaler_inc = None",
  "_device = torch.device('cpu')",
  "# Fixed incidence scaling parameters",
  "INCIDENCE_MIN = 0",
  "INCIDENCE_MAX = 5000",
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
  "def create_fixed_incidence_scaler(shape):",
  "    \"\"\"Create a MinMaxScaler with fixed parameters for incidence data\"\"\"",
  "    scaler = MinMaxScaler(feature_range=(0, 1))",
  "    # Set fixed parameters",
  "    scaler.data_min_ = np.zeros(shape)",
  "    scaler.data_max_ = np.ones(shape) * INCIDENCE_MAX",
  "    scaler.data_range_ = scaler.data_max_ - scaler.data_min_",
  "    scaler.scale_ = 1.0 / scaler.data_range_",
  "    scaler.min_ = 0 - scaler.data_min_ * scaler.scale_",
  "    return scaler",
  "",
  "def load_model(model_path, scaler_additional_path, scaler_targets_path, scaler_incidence_path=None):",
  "    global _model, _scaler_add, _scaler_tgt, _scaler_inc",
  "    try:",
  "        _scaler_add = joblib.load(scaler_additional_path)",
  "        _scaler_tgt = joblib.load(scaler_targets_path)",
  "        ",
  "        # Try to load the incidence scaler if provided",
  "        if scaler_incidence_path:",
  "            try:",
  "                _scaler_inc = joblib.load(scaler_incidence_path)",
  "                print('Loaded existing incidence scaler')",
  "            except Exception as e:",
  "                print(f'Could not load incidence scaler: {e}')",
  "                print('Will create fixed scaler during prediction')",
  "                _scaler_inc = None",
  "        else:",
  "            _scaler_inc = None",
  "            print('No incidence scaler provided, will create fixed scaler during prediction')",
  "        ",
  "        # instantiate with training hyperparams",
  "        _model = BiLSTMModel(input_dim=1, hidden_dim=160, num_layers=3, additional_dim=2, output_dim=3, dropout=0.5)",
  "        state = torch.load(model_path, map_location=_device)",
  "        _model.load_state_dict(state)",
  "        _model.to(_device).eval()",
  "        print('Model loaded successfully')",
  "    except Exception as e:",
  "        print('Error loading model/scalers:', e)",
  "        traceback.print_exc()",
  "        raise",
  "",
  "def predict(seq, additional_pair):",
  "    try:",
  "        global _scaler_inc",
  "        x = np.asarray(seq, dtype=np.float32)",
  "        if len(x.shape) == 1:",
  "            # If 1D array is passed, reshape to match expected format",
  "            x = x.reshape(1, -1, 1)",
  "        elif len(x.shape) == 2:",
  "            # If 2D array, reshape to 3D",
  "            x = x.reshape(x.shape[0], x.shape[1], 1)",
  "        ",
  "        # Check if we have a valid incidence scaler",
  "        if _scaler_inc is None:",
  "            print('Creating fixed incidence scaler with range [0, 5000]')",
  "            _scaler_inc = create_fixed_incidence_scaler(x.shape[1] if len(x.shape) > 1 else 1)",
  "        ",
  "        # Apply fixed scaling - with proper shape handling",
  "        if len(x.shape) == 3:",
  "            # Flatten to 2D for transformation",
  "            orig_shape = x.shape",
  "            flat = x.reshape(x.shape[0], x.shape[1])",
  "            try:",
  "                # Try to use the scaler",
  "                flat_scaled = _scaler_inc.transform(flat)",
  "            except Exception as e:",
  "                print(f'Scaling error: {e}')",
  "                print('Using manual scaling with range [0, 5000]')",
  "                # Manual scaling as fallback",
  "                flat_scaled = flat / INCIDENCE_MAX",
  "            # Reshape back to 3D",
  "            x_scaled = flat_scaled.reshape(orig_shape)",
  "        else:",
  "            # Directly apply scaling for 1D/2D",
  "            try:",
  "                x_scaled = _scaler_inc.transform(x)",
  "            except Exception:",
  "                x_scaled = x / INCIDENCE_MAX",
  "            x_scaled = x_scaled.reshape(-1, x_scaled.shape[-1], 1)",
  "        ",
  "        # Scale additional inputs",
  "        add_np = np.array([additional_pair], dtype=np.float32)",
  "        add_scaled = _scaler_add.transform(add_np)",
  "        ",
  "        # Convert to PyTorch tensors",
  "        x_t = torch.tensor(x_scaled, dtype=torch.float32, device=_device)",
  "        add_t = torch.tensor(add_scaled, dtype=torch.float32, device=_device)",
  "        ",
  "        # Run prediction",
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
base_dir <- "~/Desktop/Sims_calibrate/LSTM_model"
model_path            <- normalizePath(file.path(base_dir, "model4_bilstm.pt"))
scaler_additional_path<- normalizePath(file.path(base_dir, "scaler_additional.pkl"))
scaler_targets_path   <- normalizePath(file.path(base_dir, "scaler_targets.pkl"))
# Optional - can be NULL or missing if you want to use fixed scaling
scaler_incidence_path <- normalizePath(file.path(base_dir, "scaler_incidence.pkl"))

# Load model - will create fixed scaler if incidence scaler fails to load
load_model(model_path,
           scaler_additional_path,
           scaler_targets_path,
           scaler_incidence_path)

# 4) R wrapper ------------------------------------------------------------------
predict_with_bilstm <- function(time_series, n, recov) {
  stopifnot(is.numeric(time_series), is.numeric(n), is.numeric(recov))
  
  # Handle different shapes of input time series
  if (is.vector(time_series)) {
    # If it's a 1D vector
    arr <- matrix(time_series, nrow=1)
  } else if (is.matrix(time_series)) {
    # If it's already a 2D matrix
    arr <- time_series
  } else {
    stop("time_series must be a vector or matrix")
  }
  
  # Call Python predict function
  out <- predict(arr, list(n, recov))
  names(out) <- c("ptran","crate","R0")
  out
}

# 5) Test -----------------------------------------------------------------------
# Example with a random time series
set.seed(123)
example_ts <- runif(60, min = 0, max = 1000)
result <- predict_with_bilstm(example_ts, n = 5000, recov = 0.1)
print(result)

# 6) Real data application ------------------------------------------------------
# Using the model with epiworldR simulated data
library(epiworldR)
R0=3
# Run the original model to get incidence data
model_sir <- ModelSIRCONN(
  name = "COVID-19",
  prevalence = 0.01,
  n = 5000,
  contact_rate = 10,
  transmission_rate = 0.03,
  recovery_rate = 0.1
)
run(model_sir, ndays = 60)
incidence <- plot_incidence(model_sir, plot=FALSE)
infected_original <- incidence[,1]
infected_original <- infected_original[-1]  # remove day 0

# Input parameters
recov <- 0.1
n <- 5000

# Run the prediction with our model
result <- predict_with_bilstm(infected_original, n, recov)
print(result)

# Extract predicted parameters
ptran <- result[1]
crate <- result[2]
r0 <- result[3]

# Calculate the correct contact rate to maintain the desired R0
crate_true <- (r0 * recov) / ptran

# Create and run the calibrated model
calibrated_model <- ModelSIRCONN(
  name = "Calibrated COVID-19",
  prevalence = 0.01,
  n = n,
  contact_rate = crate_true,
  transmission_rate = ptran,
  recovery_rate = recov
)

run(calibrated_model, ndays = 60)
incidence_calib <- plot_incidence(calibrated_model, plot=FALSE)
infected_calibrated <- incidence_calib[,1]
infected_calibrated <- infected_calibrated[-1]  # remove day 0

# Plot comparison
library(ggplot2)
library(dplyr)

# Align lengths
min_len <- min(length(infected_original), length(infected_calibrated))
days <- 1:min_len

# Create plot data frame
df_plot <- data.frame(
  Day = rep(days, 2),
  Infected = c(infected_original[1:min_len], infected_calibrated[1:min_len]),
  Model = rep(c("Original", "Calibrated"), each = min_len)
)

# Plot
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

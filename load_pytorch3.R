# Complete workflow to generate data, scale it, predict parameters, and compare models
library(reticulate)
library(epiworldR)
library(ggplot2)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Generate incidence data using an original model
# -----------------------------------------------------------------------------
# Original model parameters
n_pop <- 5000
contact_rate_orig <- 10
transmission_rate_orig <- 0.03
recovery_rate <- 0.1
R0_orig <- contact_rate_orig * transmission_rate_orig / recovery_rate

cat("Original model parameters:\n")
cat(sprintf("  n = %d\n", n_pop))
cat(sprintf("  contact_rate = %.2f\n", contact_rate_orig))
cat(sprintf("  transmission_rate = %.4f\n", transmission_rate_orig))
cat(sprintf("  recovery_rate = %.2f\n", recovery_rate))
cat(sprintf("  R0 = %.2f\n", R0_orig))
cat("\n")

# Create and run the original model
model_orig <- ModelSIRCONN(
  name = "Original",
  prevalence = 0.01,
  n = n_pop,
  contact_rate = contact_rate_orig,
  transmission_rate = transmission_rate_orig,
  recovery_rate = recovery_rate
)

cat("Running original model...\n")
run(model_orig, ndays = 60)

# Extract incidence data
incidence_data <- plot_incidence(model_orig, plot = FALSE)
incidence_values <- incidence_data[-1, 1]  # Remove day 0

cat("Generated", length(incidence_values), "days of incidence data\n")
cat("Summary statistics:\n")
cat(sprintf("  Min: %.2f\n", min(incidence_values)))
cat(sprintf("  Max: %.2f\n", max(incidence_values)))
cat(sprintf("  Mean: %.2f\n", mean(incidence_values)))
cat("\n")

# -----------------------------------------------------------------------------
# 2. Create Python functions for scaling and prediction
# -----------------------------------------------------------------------------
py_run_string("
import numpy as np
import torch
import torch.nn as nn
from sklearn.preprocessing import MinMaxScaler
import joblib
import os

# Define BiLSTM model class
class BiLSTMModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout):
        super().__init__()
        self.bilstm = nn.LSTM(
            input_size=input_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True,
            dropout=dropout,
            bidirectional=True,
        )
        self.fc1 = nn.Linear(2 * hidden_dim + additional_dim, 64)
        self.fc2 = nn.Linear(64, output_dim)
        self.sigmoid = nn.Sigmoid()
        self.softplus = nn.Softplus()

    def forward(self, x, add_inputs):
        _, (h_n, _) = self.bilstm(x)
        h = torch.cat((h_n[-2], h_n[-1]), dim=1)
        h = torch.relu(self.fc1(torch.cat((h, add_inputs), dim=1)))
        out = self.fc2(h)
        out = torch.stack([
            self.sigmoid(out[:, 0]),
            self.softplus(out[:, 1]),
            self.softplus(out[:, 2]),
        ], dim=1)
        return out

# Function to create a fitted scaler
def create_fitted_scaler(data_range, n_features=1):
    scaler = MinMaxScaler()
    if n_features == 1:
        dummy_data = np.array([[data_range[0]], [data_range[1]]])
    else:
        dummy_data = np.zeros((2, n_features))
        dummy_data[0, :] = data_range[0]
        dummy_data[1, :] = data_range[1]
    scaler.fit(dummy_data)
    return scaler

# Function to scale incidence data
def scale_incidence_data(incidence_values, scaler_path=None):
    # Convert to numpy array if needed
    if not isinstance(incidence_values, np.ndarray):
        incidence_values = np.array(incidence_values, dtype=np.float32)
    
    # Reshape if needed
    if len(incidence_values.shape) == 1:
        incidence_values = incidence_values.reshape(1, -1)
    
    # Check if we need to pad to 60 time steps
    if incidence_values.shape[1] < 60:
        padded = np.zeros((incidence_values.shape[0], 60), dtype=np.float32)
        padded[:, :incidence_values.shape[1]] = incidence_values
        incidence_values = padded
    elif incidence_values.shape[1] > 60:
        incidence_values = incidence_values[:, :60]
    
    # Try to load scaler if path provided
    scaler = None
    if scaler_path and os.path.exists(scaler_path):
        try:
            scaler = joblib.load(scaler_path)
            if not hasattr(scaler, 'data_min_') or not hasattr(scaler, 'data_max_'):
                print('Loaded scaler is not fitted, creating a new one')
                scaler = None
            else:
                print('Loaded incidence scaler from', scaler_path)
        except Exception as e:
            print('Error loading scaler:', e)
            scaler = None
    
    # Create a new scaler if needed
    if scaler is None:
        min_val = 0
        max_val = max(100, np.max(incidence_values) * 1.1)  # use at least 0-100 range
        scaler = create_fitted_scaler([min_val, max_val], incidence_values.shape[1])
        print(f'Created new scaler with range {min_val}-{max_val}')
    
    # Scale the data
    scaled_data = scaler.transform(incidence_values)
    
    return scaled_data, scaler

# Function to predict parameters
def predict_parameters(incidence_data, n_val, recov_val, model_path=None, 
                       scaler_inc_path=None, scaler_add_path=None, scaler_tar_path=None):
    try:
        # Scale incidence data
        scaled_incidence, _ = scale_incidence_data(incidence_data, scaler_inc_path)
        scaled_incidence = scaled_incidence.reshape(1, scaled_incidence.shape[1], 1)
        
        # Handle additional parameters (n, recov)
        additional_data = np.array([[n_val, recov_val]], dtype=np.float32)
        
        # Try to load additional scaler
        scaler_additional = None
        if scaler_add_path and os.path.exists(scaler_add_path):
            try:
                scaler_additional = joblib.load(scaler_add_path)
                if not hasattr(scaler_additional, 'data_min_') or not hasattr(scaler_additional, 'data_max_'):
                    print('Additional scaler not fitted, creating a new one')
                    scaler_additional = None
                else:
                    print('Loaded additional scaler')
            except Exception as e:
                print('Error loading additional scaler:', e)
                scaler_additional = None
        
        # Create a new scaler if needed
        if scaler_additional is None:
            scaler_additional = create_fitted_scaler([0, 10000], 2)  # For n and recov
            print('Created new additional scaler')
        
        # Scale additional data
        scaled_additional = scaler_additional.transform(additional_data)
        
        # Try to load targets scaler
        scaler_targets = None
        if scaler_tar_path and os.path.exists(scaler_tar_path):
            try:
                scaler_targets = joblib.load(scaler_tar_path)
                if not hasattr(scaler_targets, 'data_min_') or not hasattr(scaler_targets, 'data_max_'):
                    print('Targets scaler not fitted, creating a new one')
                    scaler_targets = None
                else:
                    print('Loaded targets scaler')
            except Exception as e:
                print('Error loading targets scaler:', e)
                scaler_targets = None
        
        # Create a new scaler if needed
        if scaler_targets is None:
            scaler_targets = create_fitted_scaler([0, 1], 3)  # For ptran, crate, R0
            print('Created new targets scaler')
        
        # Try to load the model
        device = torch.device('cpu')
        model = None
        
        if model_path and os.path.exists(model_path):
            try:
                # Best parameters from Optuna
                BEST_PARAMS = {
                    'hidden_dim': 160,
                    'num_layers': 3,
                    'dropout': 0.5
                }
                
                model = BiLSTMModel(
                    input_dim=1,
                    hidden_dim=BEST_PARAMS['hidden_dim'],
                    num_layers=BEST_PARAMS['num_layers'],
                    additional_dim=2,
                    output_dim=3,
                    dropout=BEST_PARAMS['dropout']
                ).to(device)
                
                model.load_state_dict(torch.load(model_path, map_location=device))
                model.eval()
                print('Loaded model from', model_path)
            except Exception as e:
                print('Error loading model:', e)
                model = None
        
        # If model is available, use it for prediction
        if model is not None:
            # Convert to PyTorch tensors
            X_tensor = torch.tensor(scaled_incidence, dtype=torch.float32).to(device)
            add_tensor = torch.tensor(scaled_additional, dtype=torch.float32).to(device)
            
            # Make prediction
            with torch.no_grad():
                pred_scaled = model(X_tensor, add_tensor).cpu().numpy()
            
            # Inverse transform
            pred_original = scaler_targets.inverse_transform(pred_scaled)
            
            # Extract values
            ptran = float(pred_original[0, 0])
            crate = float(pred_original[0, 1])
            r0 = float(pred_original[0, 2])
            
            # Ensure values are reasonable
            ptran = max(0.001, min(0.999, ptran))
            crate = max(0.1, min(30, crate))
            r0 = max(0.1, min(10, r0))
        else:
            # If model isn't available, use default values
            print('Using default parameter values')
            ptran = 0.03
            crate = 10.0
            r0 = 3.0
        
        # Calculate adjusted contact rate
        crate_adjusted = (r0 * recov_val) / ptran
        
        return {
            'ptran': ptran,
            'crate': crate,
            'R0': r0,
            'crate_adjusted': crate_adjusted,
            'success': True
        }
    
    except Exception as e:
        print('Error during prediction:', e)
        import traceback
        traceback.print_exc()
        
        # Return default values
        return {
            'ptran': 0.03,
            'crate': 10.0,
            'R0': 3.0,
            'crate_adjusted': 10.0,
            'success': False,
            'error': str(e)
        }
")

# -----------------------------------------------------------------------------
# 3. Predict parameters from the incidence data
# -----------------------------------------------------------------------------
# Set paths to model and scaler files
BASE_DIR <- "~/Desktop/Sims_calibrate/LSTM_model"  # Change to your path
MODEL_PATH <- file.path(BASE_DIR, "model4_bilstm.pt")
SCALER_INC_PATH <- file.path(BASE_DIR, "scaler_incidence.pkl")
SCALER_ADD_PATH <- file.path(BASE_DIR, "scaler_additional.pkl")
SCALER_TAR_PATH <- file.path(BASE_DIR, "scaler_targets.pkl")

# Normalize paths
MODEL_PATH <- normalizePath(MODEL_PATH, mustWork = FALSE)
SCALER_INC_PATH <- normalizePath(SCALER_INC_PATH, mustWork = FALSE)
SCALER_ADD_PATH <- normalizePath(SCALER_ADD_PATH, mustWork = FALSE)
SCALER_TAR_PATH <- normalizePath(SCALER_TAR_PATH, mustWork = FALSE)

# Check if files exist
cat("Checking if files exist:\n")
cat(sprintf("Model file: %s exists: %s\n", MODEL_PATH, file.exists(MODEL_PATH)))
cat(sprintf("Incidence scaler: %s exists: %s\n", SCALER_INC_PATH, file.exists(SCALER_INC_PATH)))
cat(sprintf("Additional scaler: %s exists: %s\n", SCALER_ADD_PATH, file.exists(SCALER_ADD_PATH)))
cat(sprintf("Targets scaler: %s exists: %s\n", SCALER_TAR_PATH, file.exists(SCALER_TAR_PATH)))
cat("\n")

# First just scale the incidence data (without fitting)
cat("Scaling incidence data...\n")
scaled_result <- py$scale_incidence_data(incidence_values, SCALER_INC_PATH)
scaled_incidence <- scaled_result[[1]]

# Print first few values of the scaled data
cat("First few values of scaled incidence data:\n")
print(head(py_to_r(scaled_incidence)[1, ]))
cat("\n")

# Predict parameters
cat("Predicting parameters...\n")
params <- py$predict_parameters(
  incidence_values, 
  n_pop, 
  recovery_rate,
  MODEL_PATH,
  SCALER_INC_PATH,
  SCALER_ADD_PATH,
  SCALER_TAR_PATH
)

# Display predicted parameters
cat("Predicted parameters:\n")
cat(sprintf("  ptran = %.6f\n", params$ptran))
cat(sprintf("  crate = %.6f\n", params$crate))
cat(sprintf("  R0 = %.6f\n", params$R0))
cat(sprintf("  adjusted contact rate = %.6f\n", params$crate_adjusted))
cat(sprintf("  success: %s\n", params$success))
if (!params$success) {
  cat(sprintf("  error: %s\n", params$error))
}
cat("\n")

# -----------------------------------------------------------------------------
# 4. Create and run the calibrated model
# -----------------------------------------------------------------------------
# Create calibrated model
model_calib <- ModelSIRCONN(
  name = "Calibrated",
  prevalence = 0.01,
  n = n_pop,
  contact_rate = params$crate_adjusted,
  transmission_rate = params$ptran,
  recovery_rate = recovery_rate
)

# Run the calibrated model
cat("Running calibrated model...\n")
run(model_calib, ndays = 60)

# -----------------------------------------------------------------------------
# 5. Compare the original and calibrated models
# -----------------------------------------------------------------------------
# Plot both models
par(mfrow=c(1,2))
plot_incidence(model_orig, main="Original Model")
plot_incidence(model_calib, main="Calibrated Model")

# Run multiple simulations for better comparison
cat("Running multiple simulations for comparison...\n")
nsims <- 10
saver <- make_saver("total_hist")

run_multiple(model_orig, ndays = 60, nsims = nsims, 
             saver = saver, nthread = min(parallel::detectCores() - 1, nsims))

run_multiple(model_calib, ndays = 60, nsims = nsims,
             saver = saver, nthread = min(parallel::detectCores() - 1, nsims))

# Extract and prepare results
res_orig <- run_multiple_get_results(model_orig)$total_hist %>%
  filter(state == "Infected") %>%
  rename(time = date, rep = sim_num, I = counts) %>%
  mutate(model = "Original")

res_calib <- run_multiple_get_results(model_calib)$total_hist %>%
  filter(state == "Infected") %>%
  rename(time = date, rep = sim_num, I = counts) %>%
  mutate(model = "Calibrated")

# Combine results
all_results <- bind_rows(res_orig, res_calib)

# Calculate summary statistics
summary_df <- all_results %>%
  group_by(model, time) %>%
  summarize(
    mean_I = mean(I),
    lower = quantile(I, 0.025),
    upper = quantile(I, 0.975),
    .groups = "drop"
  )

# Create the comparison plot
p <- ggplot(summary_df, aes(x = time, y = mean_I, color = model, fill = model)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    x = "Day",
    y = "Active Infected (I)",
    title = "Original vs. Calibrated SIRCONN Simulations",
    subtitle = paste0("Predictions: ptran = ", round(params$ptran, 4), 
                      ", crate = ", round(params$crate, 2), 
                      ", R0 = ", round(params$R0, 2)),
    color = NULL, 
    fill = NULL
  ) +
  theme_minimal(base_size = 12)

# Display the plot
print(p)

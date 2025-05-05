library(reticulate)

torch <- import("torch")

model_path <- normalizePath("~/Desktop/Sims_calibrate/model4_bilstm_new.pt")

library(reticulate)

torch <- import("torch")
device <- torch$device("cpu")  # force CPU
model <- torch$load(model_path, map_location = device)

# Load required R packages
library(reticulate)

# Import Python modules
torch <- import("torch")
nn <- torch$nn
joblib <- import("joblib")
np <- import("numpy")
os <- import("os")

# Define paths (adjust to your actual paths)
model_path <- normalizePath("~/Desktop/Sims_calibrate/model4_bilstm_new.pt")
# Assume scalers are in the same directory or specify full paths
scaler_additional_path <- normalizePath("~/Desktop/Sims_calibrate/scaler_additional.pkl")
scaler_targets_path <- normalizePath("~/Desktop/Sims_calibrate/scaler_targets.pkl")

# Now define the model in Python with correct paths
py_run_string('
import torch
import torch.nn as nn
import joblib
import numpy as np
import os

# Define Model Architecture
class BiLSTMModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout):
        super(BiLSTMModel, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        
        # Bidirectional LSTM
        self.bilstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True, 
                             dropout=dropout, bidirectional=True)
        
        # Fully connected layers
        self.fc1 = nn.Linear(2 * hidden_dim + additional_dim, 64)
        self.fc2 = nn.Linear(64, output_dim)

        self.sigmoid = nn.Sigmoid()
        self.softplus = nn.Softplus()

    def forward(self, x, additional_inputs):
        _, (h_n, _) = self.bilstm(x)
        h_n = torch.cat((h_n[-2], h_n[-1]), dim=1)
        combined = torch.cat((h_n, additional_inputs), dim=1)
        x = self.fc1(combined)
        out = self.fc2(x)
        out = torch.stack([
            self.sigmoid(out[:, 0]),
            self.softplus(out[:, 1]),
            self.softplus(out[:, 2])
        ], dim=1)
        return out
')

# Pass the paths correctly to Python
py_run_string(sprintf('
# Model parameters
input_dim = 1
hidden_dim = 128
num_layers = 2
additional_dim = 2
output_dim = 3
dropout = 0.3

# Set device
device = torch.device("cpu")  # Always use CPU to avoid issues

# Initialize model
model = BiLSTMModel(input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout)

# Load model weights - make sure to use the correct path
model.load_state_dict(torch.load("%s", map_location=device))
model.to(device)
model.eval()

# Load scalers - make sure to use the correct paths
scaler_additional = joblib.load("%s")
scaler_targets = joblib.load("%s")
', model_path, scaler_additional_path, scaler_targets_path))

# Create an R function to predict with the model
predict_with_bilstm <- function(time_series, n, recov) {
  # Check inputs
  if (!is.numeric(time_series) || !is.numeric(n) || !is.numeric(recov)) {
    stop("Inputs must be numeric")
  }
  
  # Reshape time series to match expected input (batch_size, seq_len, features)
  time_series <- array(time_series, dim = c(1, length(time_series), 1))
  
  # Pass the time_series to Python environment
  py$time_series <- time_series
  
  # Run prediction using the Python model
  py_run_string(sprintf('
  # Scale additional inputs
  new_add = np.array([[%f, %f]])
  new_add_scaled = scaler_additional.transform(new_add)
  
  # Convert inputs to tensors - use the time_series variable directly
  new_seq_tensor = torch.tensor(time_series, dtype=torch.float32).to(device)
  new_add_tensor = torch.tensor(new_add_scaled, dtype=torch.float32).to(device)
  
  # Run model
  with torch.no_grad():
      pred = model(new_seq_tensor, new_add_tensor).cpu().numpy()
  
  # Inverse-transform predictions
  pred_unscaled = scaler_targets.inverse_transform(pred)
  ', n, recov))
  
  # Extract predictions
  predictions <- py$pred_unscaled[1, ]
  names(predictions) <- c("ptran", "crate", "R0")
  
  return(predictions)
}

predict_with_bilstm <- function(time_series, n, recov) {
  # Check inputs
  if (!is.numeric(time_series) || !is.numeric(n) || !is.numeric(recov)) {
    stop("Inputs must be numeric")
  }
  
  # Reshape time series to match expected input (batch_size, seq_len, features)
  time_series <- array(time_series, dim = c(1, length(time_series), 1))
  
  # Pass all data to Python
  py$time_series <- time_series
  py$n_value <- n
  py$recov_value <- recov
  
  # Run the prediction in one clean Python block
  py_run_string("
# Scale additional inputs
new_add = np.array([[n_value, recov_value]])
new_add_scaled = scaler_additional.transform(new_add)

# Convert inputs to tensors
new_seq_tensor = torch.tensor(time_series, dtype=torch.float32).to(device)
new_add_tensor = torch.tensor(new_add_scaled, dtype=torch.float32).to(device)

# Run model
with torch.no_grad():
    pred = model(new_seq_tensor, new_add_tensor).cpu().numpy()

# Inverse-transform predictions
pred_unscaled = scaler_targets.inverse_transform(pred)
")
  
  # Extract predictions
  predictions <- py$pred_unscaled[1, ]
  names(predictions) <- c("ptran", "crate", "R0")
  
  return(predictions)
}

# Example usage
set.seed(123)
example_time_series <- rnorm(60)
example_n <- 5000
example_recov <- 0.1

# Make prediction
predictions <- predict_with_bilstm(example_time_series, example_n, example_recov)
print("Predictions:")
print(predictions)

# Debug R script with detailed logging
library(reticulate)

# Updated Python helper with debugging
python_code <- "
import torch
import torch.nn as nn
import joblib
import numpy as np
from sklearn.preprocessing import MinMaxScaler

_model = None
_scaler_add = None
_scaler_tgt = None
_scaler_inc = None
_device = torch.device('cpu')

INCIDENCE_MIN = 0
INCIDENCE_MAX = 10000

class BiLSTMModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout):
        super().__init__()
        self.bilstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True, dropout=dropout, bidirectional=True)
        self.fc1 = nn.Linear(2 * hidden_dim + additional_dim, 64)
        self.fc2 = nn.Linear(64, output_dim)
        self.sigmoid = nn.Sigmoid()
        self.softplus = nn.Softplus()

    def forward(self, x, additional_inputs):
        _, (h_n, _) = self.bilstm(x)
        hid = torch.cat((h_n[-2], h_n[-1]), dim=1)
        combined = torch.cat((hid, additional_inputs), dim=1)
        
        # CRITICAL: Your R code is missing ReLU activation!
        x = torch.relu(self.fc1(combined))  # Added ReLU here
        out = self.fc2(x)
        return torch.stack([
            self.sigmoid(out[:, 0]),
            self.softplus(out[:, 1]),
            self.softplus(out[:, 2])
        ], dim=1)

def create_fixed_incidence_scaler(shape):
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.data_min_ = np.zeros(shape)
    scaler.data_max_ = np.ones(shape) * INCIDENCE_MAX
    scaler.data_range_ = scaler.data_max_ - scaler.data_min_
    scaler.scale_ = 1.0 / scaler.data_range_
    scaler.min_ = 0 - scaler.data_min_ * scaler.scale_
    return scaler

def load_model(model_path, scaler_add_path, scaler_tgt_path, scaler_inc_path=None):
    global _model, _scaler_add, _scaler_tgt, _scaler_inc

    _scaler_add = joblib.load(scaler_add_path)
    _scaler_tgt = joblib.load(scaler_tgt_path)

    # Debug scaler information
    print('=== SCALER DEBUGGING (R) ===')
    print(f'Additional scaler data_min_: {_scaler_add.data_min_}')
    print(f'Additional scaler data_max_: {_scaler_add.data_max_}')
    print(f'Additional scaler scale_: {_scaler_add.scale_}')
    print(f'Additional scaler min_: {_scaler_add.min_}')
    
    print(f'Targets scaler data_min_: {_scaler_tgt.data_min_}')
    print(f'Targets scaler data_max_: {_scaler_tgt.data_max_}')
    print(f'Targets scaler scale_: {_scaler_tgt.scale_}')
    print(f'Targets scaler min_: {_scaler_tgt.min_}')

    if scaler_inc_path:
        try:
            _scaler_inc = joblib.load(scaler_inc_path)
            print(f'Loaded incidence scaler data_min_[:5]: {_scaler_inc.data_min_[:5]}')
            print(f'Loaded incidence scaler data_max_[:5]: {_scaler_inc.data_max_[:5]}')
            print(f'Loaded incidence scaler scale_[:5]: {_scaler_inc.scale_[:5]}')
            print(f'Loaded incidence scaler min_[:5]: {_scaler_inc.min_[:5]}')
        except:
            _scaler_inc = None
            print('Failed to load incidence scaler, will create fixed one')
    else:
        _scaler_inc = None

    _model = BiLSTMModel(input_dim=1, hidden_dim=160, num_layers=3,
                         additional_dim=2, output_dim=3, dropout=0.5)
    state = torch.load(model_path, map_location=_device)
    _model.load_state_dict(state)
    _model.to(_device).eval()
    
    print(f'Model fc1 weight shape: {_model.fc1.weight.shape}')
    print(f'Model fc1 bias shape: {_model.fc1.bias.shape}')
    print(f'Model fc2 weight shape: {_model.fc2.weight.shape}')
    print(f'Model fc2 bias shape: {_model.fc2.bias.shape}')

def predict_debug(seq, additional_pair):
    global _scaler_inc
    x = np.asarray(seq, dtype=np.float32).reshape(1, -1, 1)
    
    print(f'\\n=== INPUT DEBUGGING (R) ===')
    print(f'Input shape: {x.shape}')
    print(f'Input first 10 values: {x[0, :10, 0]}')
    print(f'Input last 10 values: {x[0, -10:, 0]}')

    if _scaler_inc is None:
        _scaler_inc = create_fixed_incidence_scaler(x.shape[1])
        print(f'Created fixed scaler with data_min_[:5]: {_scaler_inc.data_min_[:5]}')
        print(f'Created fixed scaler with data_max_[:5]: {_scaler_inc.data_max_[:5]}')
        print(f'Created fixed scaler with scale_[:5]: {_scaler_inc.scale_[:5]}')
        print(f'Created fixed scaler with min_[:5]: {_scaler_inc.min_[:5]}')

    x_scaled = _scaler_inc.transform(x.reshape(1, -1)).reshape(1, -1, 1)
    print(f'\\n=== SCALING DEBUGGING (R) ===')
    print(f'Scaled input first 10 values: {x_scaled[0, :10, 0]}')
    print(f'Scaled input last 10 values: {x_scaled[0, -10:, 0]}')
    
    add_np = np.array([additional_pair], dtype=np.float32)
    add_scaled = _scaler_add.transform(add_np)
    print(f'Additional inputs (raw): {add_np}')
    print(f'Additional inputs (scaled): {add_scaled}')

    x_t = torch.tensor(x_scaled, dtype=torch.float32, device=_device)
    add_t = torch.tensor(add_scaled, dtype=torch.float32, device=_device)

    print(f'\\n=== PREDICTION DEBUGGING (R) ===')
    with torch.no_grad():
        # Get LSTM output
        _, (h_n, _) = _model.bilstm(x_t)
        hid = torch.cat((h_n[-2], h_n[-1]), dim=1)
        print(f'LSTM hidden state shape: {hid.shape}')
        print(f'LSTM hidden state first 5 values: {hid[0, :5]}')
        
        # Combine with additional inputs
        combined = torch.cat((hid, add_t), dim=1)
        print(f'Combined input shape: {combined.shape}')
        print(f'Combined input first 5 values: {combined[0, :5]}')
        
        # First linear layer with ReLU
        fc1_out = _model.fc1(combined)
        print(f'FC1 output first 5 values: {fc1_out[0, :5]}')
        
        relu_out = torch.relu(fc1_out)
        print(f'ReLU output first 5 values: {relu_out[0, :5]}')
        
        # Second linear layer
        fc2_out = _model.fc2(relu_out)
        print(f'FC2 output values: {fc2_out[0]}')
        
        # Apply activations
        out = torch.stack([
            _model.sigmoid(fc2_out[:, 0]),
            _model.softplus(fc2_out[:, 1]),
            _model.softplus(fc2_out[:, 2])
        ], dim=1)
        
        print(f'Prediction (scaled): {out[0]}')
        
        result = _scaler_tgt.inverse_transform(out.cpu().numpy())[0].tolist()
        print(f'Prediction (original scale): {result}')
        
        return result

def predict(seq, additional_pair):
    # Use the debug version for now
    return predict_debug(seq, additional_pair)
"

writeLines(python_code, "model_utils_debug.py")
source_python("model_utils_debug.py")

# Load model with debugging
base_dir <- normalizePath("~/Desktop/Sims_calibrate/LSTM_model")
model_path            <- normalizePath(file.path(base_dir, "model4_bilstm.pt"))
scaler_add_path       <- normalizePath(file.path(base_dir, "scaler_additional.pkl"))
scaler_tgt_path       <- normalizePath(file.path(base_dir, "scaler_targets.pkl"))
scaler_inc_path       <- normalizePath(file.path(base_dir, "scaler_incidence.pkl"))

load_model(
  model_path       = model_path,
  scaler_add_path  = scaler_add_path,
  scaler_tgt_path  = scaler_tgt_path,
  scaler_inc_path  = scaler_inc_path
)

# Define wrapper function
predict_with_bilstm_debug <- function(time_series, n, recov) {
  stopifnot(length(time_series) == 61, is.numeric(n), is.numeric(recov))
  time_series <- as.numeric(time_series)
  out <- predict(time_series, list(n, recov))
  names(out) <- c("ptran", "crate", "R0")
  return(out)
}

# Test prediction
incidence_vec <- c(
  103, 37, 60, 74, 108, 125, 138, 186, 215, 276, 318, 331, 414, 402, 446, 454, 405, 401, 373, 334, 285,
  241, 219, 156, 140, 108, 93, 82, 73, 78, 48, 38, 34, 22, 27, 22, 20, 11, 14, 8, 11, 14,
  8, 6, 6, 2, 0, 7, 2, 4, 4, 1, 1, 2, 2, 1, 0, 2, 1, 1, 0
)

n <- 7087
recov <- 0.203

result <- predict_with_bilstm_debug(incidence_vec, n, recov)
cat("\n=== FINAL RESULTS (R) ===\n")
print(result)

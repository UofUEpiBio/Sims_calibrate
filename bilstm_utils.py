# -*- coding: utf-8 -*-
import os
import joblib
import numpy as np
import torch
import torch.nn as nn

# Model definition (must match training architecture exactly)
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

# Best parameters from Optuna (needed to instantiate model correctly)
BEST_PARAMS = {
    'hidden_dim': 160,
    'num_layers': 3,
    'dropout': 0.5,
    'lr': 0.0002766923888203664,
    'lambda_penalty': 0.0001773491932427762
}

# Global variables to store model and scalers
model = None
scaler_additional = None
scaler_targets = None
scaler_incidence = None
device = torch.device('cpu')

def load_model_and_scalers(model_path, scaler_add_path, scaler_tar_path, scaler_inc_path):
    """Load the BiLSTM model and all necessary scalers."""
    global model, scaler_additional, scaler_targets, scaler_incidence, device
    
    # Load scalers
    scaler_additional = joblib.load(scaler_add_path)
    scaler_targets = joblib.load(scaler_tar_path)
    
    try:
        scaler_incidence = joblib.load(scaler_inc_path)
        print('Incidence scaler loaded successfully')
    except Exception as e:
        print(f'Warning: Could not load incidence scaler: {e}')
        print('Will use raw values for incidence data')
        scaler_incidence = None
    
    # Instantiate model with correct hyperparameters
    model = BiLSTMModel(
        input_dim=1,
        hidden_dim=BEST_PARAMS['hidden_dim'],
        num_layers=BEST_PARAMS['num_layers'],
        additional_dim=2,
        output_dim=3,
        dropout=BEST_PARAMS['dropout']
    ).to(device)
    
    # Load trained weights
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()  # Set to evaluation mode
    print('Model loaded successfully')
    
    return True

def predict(incidence_data, n_val, recov_val):
    """
    Make predictions using the loaded model.
    
    Args:
        incidence_data: numpy array of shape (batch, seq_len) or (seq_len,)
        n_val: population size parameter
        recov_val: recovery rate parameter
        
    Returns:
        Array of [ptran, crate, R0] in original units
    """
    global model, scaler_additional, scaler_targets, scaler_incidence, device
    
    # Ensure model is loaded
    if model is None:
        raise ValueError('Model not loaded. Call load_model_and_scalers first.')
    
    # Handle single sequence case
    if len(incidence_data.shape) == 1:
        incidence_data = incidence_data.reshape(1, -1)
    
    # Scale incidence data if scaler is available
    if scaler_incidence is not None:
        scaled_incidence = scaler_incidence.transform(incidence_data)
    else:
        scaled_incidence = incidence_data
    
    # Reshape to (batch, seq_len, features)
    scaled_incidence = scaled_incidence.reshape(incidence_data.shape[0], incidence_data.shape[1], 1)
    
    # Scale additional inputs
    additional_data = np.array([[n_val, recov_val]])
    scaled_additional = scaler_additional.transform(additional_data)
    
    # Convert to tensors
    X_tensor = torch.tensor(scaled_incidence, dtype=torch.float32).to(device)
    add_tensor = torch.tensor(scaled_additional, dtype=torch.float32).to(device)
    
    # Make prediction
    with torch.no_grad():
        pred_scaled = model(X_tensor, add_tensor).cpu().numpy()
    
    # Inverse transform to original units
    pred_original = scaler_targets.inverse_transform(pred_scaled)
    
    return pred_original[0]  # Return the first (only) prediction


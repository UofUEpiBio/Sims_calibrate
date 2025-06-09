# model_utils.py
import torch
import torch.nn as nn
import joblib
import numpy as np
import traceback
from sklearn.preprocessing import MinMaxScaler

_model = None
_scaler_add = None
_scaler_tgt = None
_scaler_inc = None
_device = torch.device('cpu')
# Fixed incidence scaling parameters
INCIDENCE_MIN = 0
INCIDENCE_MAX = 5000

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
        x = self.fc1(combined)
        out = self.fc2(x)
        return torch.stack([
            self.sigmoid(out[:,0]),
            self.softplus(out[:,1]),
            self.softplus(out[:,2])
        ], dim=1)

def create_fixed_incidence_scaler(shape):
    """Create a MinMaxScaler with fixed parameters for incidence data"""
    scaler = MinMaxScaler(feature_range=(0, 1))
    # Set fixed parameters
    scaler.data_min_ = np.zeros(shape)
    scaler.data_max_ = np.ones(shape) * INCIDENCE_MAX
    scaler.data_range_ = scaler.data_max_ - scaler.data_min_
    scaler.scale_ = 1.0 / scaler.data_range_
    scaler.min_ = 0 - scaler.data_min_ * scaler.scale_
    return scaler

def load_model(model_path, scaler_additional_path, scaler_targets_path, scaler_incidence_path=None):
    global _model, _scaler_add, _scaler_tgt, _scaler_inc
    try:
        _scaler_add = joblib.load(scaler_additional_path)
        _scaler_tgt = joblib.load(scaler_targets_path)
        
        # Try to load the incidence scaler if provided
        if scaler_incidence_path:
            try:
                _scaler_inc = joblib.load(scaler_incidence_path)
                print('Loaded existing incidence scaler')
            except Exception as e:
                print(f'Could not load incidence scaler: {e}')
                print('Will create fixed scaler during prediction')
                _scaler_inc = None
        else:
            _scaler_inc = None
            print('No incidence scaler provided, will create fixed scaler during prediction')
        
        # instantiate with training hyperparams
        _model = BiLSTMModel(input_dim=1, hidden_dim=160, num_layers=3, additional_dim=2, output_dim=3, dropout=0.5)
        state = torch.load(model_path, map_location=_device)
        _model.load_state_dict(state)
        _model.to(_device).eval()
        print('Model loaded successfully')
    except Exception as e:
        print('Error loading model/scalers:', e)
        traceback.print_exc()
        raise

def predict(seq, additional_pair):
    try:
        global _scaler_inc
        x = np.asarray(seq, dtype=np.float32)
        if len(x.shape) == 1:
            # If 1D array is passed, reshape to match expected format
            x = x.reshape(1, -1, 1)
        elif len(x.shape) == 2:
            # If 2D array, reshape to 3D
            x = x.reshape(x.shape[0], x.shape[1], 1)
        
        # Check if we have a valid incidence scaler
        if _scaler_inc is None:
            print('Creating fixed incidence scaler with range [0, 5000]')
            _scaler_inc = create_fixed_incidence_scaler(x.shape[1] if len(x.shape) > 1 else 1)
        
        # Apply fixed scaling - with proper shape handling
        if len(x.shape) == 3:
            # Flatten to 2D for transformation
            orig_shape = x.shape
            flat = x.reshape(x.shape[0], x.shape[1])
            try:
                # Try to use the scaler
                flat_scaled = _scaler_inc.transform(flat)
            except Exception as e:
                print(f'Scaling error: {e}')
                print('Using manual scaling with range [0, 5000]')
                # Manual scaling as fallback
                flat_scaled = flat / INCIDENCE_MAX
            # Reshape back to 3D
            x_scaled = flat_scaled.reshape(orig_shape)
        else:
            # Directly apply scaling for 1D/2D
            try:
                x_scaled = _scaler_inc.transform(x)
            except Exception:
                x_scaled = x / INCIDENCE_MAX
            x_scaled = x_scaled.reshape(-1, x_scaled.shape[-1], 1)
        
        # Scale additional inputs
        add_np = np.array([additional_pair], dtype=np.float32)
        add_scaled = _scaler_add.transform(add_np)
        
        # Convert to PyTorch tensors
        x_t = torch.tensor(x_scaled, dtype=torch.float32, device=_device)
        add_t = torch.tensor(add_scaled, dtype=torch.float32, device=_device)
        
        # Run prediction
        with torch.no_grad():
            out = _model(x_t, add_t).cpu().numpy()
        return _scaler_tgt.inverse_transform(out)[0].tolist()
    except Exception as e:
        print('Prediction error:', e)
        traceback.print_exc()
        raise

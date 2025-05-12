# model_utils.py
import torch
import torch.nn as nn
import joblib
import numpy as np
import traceback

_model = None
_scaler_add = None
_scaler_tgt = None
_scaler_inc = None
_device = torch.device('cpu')

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

def load_model(model_path, scaler_additional_path, scaler_targets_path, scaler_incidence_path):
    global _model, _scaler_add, _scaler_tgt, _scaler_inc
    try:
        _scaler_add = joblib.load(scaler_additional_path)
        _scaler_tgt = joblib.load(scaler_targets_path)
        _scaler_inc = joblib.load(scaler_incidence_path)
        # instantiate with training hyperparams
        _model = BiLSTMModel(input_dim=1, hidden_dim=160, num_layers=3, additional_dim=2, output_dim=3, dropout=0.5)
        state = torch.load(model_path, map_location=_device)
        _model.load_state_dict(state)
        _model.to(_device).eval()
    except Exception as e:
        print('Error loading model/scalers:', e)
        traceback.print_exc()
        raise

def predict(seq, additional_pair):
    try:
        x = np.asarray(seq, dtype=np.float32)
        flat = x.reshape(x.shape[0], x.shape[1])
        # fallback: if incidence scaler not fitted, use raw
        try:
            flat_scaled = _scaler_inc.transform(flat)
        except Exception:
            flat_scaled = flat
        x_scaled = flat_scaled.reshape(x.shape[0], x.shape[1], 1)
        add_np = np.array([additional_pair], dtype=np.float32)
        add_scaled = _scaler_add.transform(add_np)
        x_t = torch.tensor(x_scaled, dtype=torch.float32, device=_device)
        add_t = torch.tensor(add_scaled, dtype=torch.float32, device=_device)
        with torch.no_grad():
            out = _model(x_t, add_t).cpu().numpy()
        return _scaler_tgt.inverse_transform(out)[0].tolist()
    except Exception as e:
        print('Prediction error:', e)
        traceback.print_exc()
        raise

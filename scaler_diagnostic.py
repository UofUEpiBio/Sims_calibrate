# scaler_diagnostic.py
import joblib
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler

def analyze_scaler(file_path, name):
    try:
        print(f'\nAnalyzing {name}:')
        scaler = joblib.load(file_path)
        print(f'  Type: {type(scaler).__name__}')
        
        # Check attributes for different scaler types
        if isinstance(scaler, MinMaxScaler):
            print(f'  min_: {getattr(scaler, "min_", "Not found")}')
            print(f'  scale_: {getattr(scaler, "scale_", "Not found")}')
            print(f'  data_min_: {getattr(scaler, "data_min_", "Not found")}')
            print(f'  data_max_: {getattr(scaler, "data_max_", "Not found")}')
        elif isinstance(scaler, StandardScaler):
            print(f'  mean_: {getattr(scaler, "mean_", "Not found")}')
            print(f'  scale_: {getattr(scaler, "scale_", "Not found")}')
            print(f'  var_: {getattr(scaler, "var_", "Not found")}')
        
        # Try a simple transform test
        try:
            if isinstance(scaler, MinMaxScaler):
                test_data = np.array([[50.0]])
            else:
                test_data = np.array([[0.5, 0.1]])
            
            result = scaler.transform(test_data)
            print(f'  ✓ Transform test successful: {result}')
        except Exception as e:
            print(f'  ✗ Transform test failed: {e}')
        
        return scaler
    except Exception as e:
        print(f'  ✗ Error loading {name}: {e}')
        return None

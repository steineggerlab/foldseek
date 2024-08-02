import torch
import numpy as np

device = torch.device('cpu')
state = torch.load("model.pt", map_location=device)

with open("weights.bin", "wb") as f:
    for k, v in state["state_dict"].items():
        shape = list(v.shape)
        while len(shape) < 4:
            shape.append(-1)
        
        f.write(np.array(shape, dtype=np.int32).tobytes())
        f.write(v.flatten().numpy().astype(np.float32).tobytes())

import numpy as np 
def array_d(*args, **kwargs):
    # kwargs.setdefault("dtype", np.double)
    return np.array(*args, **kwargs)
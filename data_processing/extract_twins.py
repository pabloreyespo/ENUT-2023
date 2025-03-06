
import pandas as pd
import numpy as np

twins = np.load("data/raw/matriz_gemelos5.npy")
twins = pd.DataFrame(twins)
twins.to_csv(
    "data/raw/matriz_gemelos5.csv.gzip",
    header=False,
    index=False,
    compression="gzip",
    float_format='%.4f'
)
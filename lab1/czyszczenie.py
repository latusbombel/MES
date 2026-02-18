import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

with open('energie.dat', 'r') as f:
    raw_content = f.read()


blocks = raw_content.strip().split('\n\n')

data_cleaned = []
for block in blocks:
    numbers = [float(line) for line in block.split('\n') if line.strip()]
    data_cleaned.append(numbers)

df = pd.DataFrame(data_cleaned).T

df.to_dat("poukladane_dane.dat", index=False, sep=' ')
import pandas as pd
import sys
import numpy as np

cols = ["GPSD", "GPST", "lat", "lon", "height", "Q", "ns", "sdn", "sde", "sdu", "sdne", "sdeu", "sdun", "age", "ratio"]
df = pd.read_csv(sys.argv[1], names=cols, comment='%', delim_whitespace=True)

print("Lat:", df.lat.mean())
print("Lon:", df.lon.mean())
print("Alt:", df.height.mean())

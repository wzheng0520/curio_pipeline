import pandas as pd
from peak_detector import peak_seeker
import sys

# file directory as input
file_directory = sys.argv[1]

dist = int(sys.argv[2])
pro_sens = int(sys.argv[3])
pro_min = int(sys.argv[4])


df = pd.read_csv(file_directory)
data = df.x.tolist()
peaks = peak_seeker(data, dist, pro_sens, pro_min)

# this serves as the returned value, number of peaks in string format
print(len(peaks))

from read_plot import *
import matplotlib.pyplot as plt

t, w = read_plot("omega.data")

plt.scatter(t,w, s=1)
plt.show()

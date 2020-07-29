import numpy as np
import matplotlib.pyplot as plt

plt.close("all")

data = np.random.randn(100000)

fig1 = plt.figure(1)
fig1.canvas.set_window_title("histogram_example")
ax11 = plt.subplot(221)
ax11.set_title("10 bins")
ax12 = plt.subplot(222)
ax12.set_title("50 bins")
ax13 = plt.subplot(223)
ax13.set_title("100 bins")
ax14 = plt.subplot(224)
ax14.set_title("500 bins")

ax11.hist(data, bins=10)
ax12.hist(data, bins=50)
ax13.hist(data, bins=100)
ax14.hist(data, bins=500)

plt.show(block=False)

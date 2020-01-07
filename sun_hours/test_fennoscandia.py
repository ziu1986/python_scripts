import numpy as np
import matplotlib.pyplot as plt
from sun_hours import *
from mytools.met_tools import plot_month_span, plot_month_name, print_all
from cycler import cycler

station_labels = ("Jergul", "Karasjok", "Svanvik", "Pallas", "Esrange")
day_length_stations = [ daylength(doy, lat) for doy in range(1,366) for lat in (69.45,69.467,69.45,69.97,67.83)]
day_length_stations = np.array(day_length_stations).reshape(365,5)

# Plot it
plt.close("all")

fig1 = plt.figure(1, figsize=(9,6))
fig1.canvas.set_window_title("test_fennoscandia")
ax11 = plt.subplot()

plots = ax11.plot(day_length_stations)
ax11.set_xlabel("Day of Year")
ax11.set_ylabel("Day Length (hours)")
ax11.legend(plots, station_labels)
plot_month_span(ax11)
plot_month_name(ax11, ypos=24.5)

# Show it
plt.show(block=False)

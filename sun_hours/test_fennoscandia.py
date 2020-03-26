import numpy as np
import matplotlib.pyplot as plt
from sun_hours import *
from mytools.plot_tools import plot_month_span, plot_month_name, print_all
from mytools.station_info import station_location
from cycler import cycler

station_labels = ("Jergul", "Karasjok", "Svanvik", "Pallas", "Esrange", "Prestebakke")
day_length_stations = [ daylength(doy, station_location[each].lat) for doy in range(1,366) for each in station_labels]
day_length_stations = np.array(day_length_stations).reshape(365,len(station_labels))

# Plot it
plt.close("all")
custom_cycler = (cycler(color=['orange', 'orange', 'blueviolet', 'black', 'blue', 'red']) +
cycler(linestyle=['-', '--', ':', '-.', '-', '--']))

fig1 = plt.figure(1, figsize=(9,6))
fig1.canvas.set_window_title("test_fennoscandia")
ax11 = plt.subplot()
ax11.set_prop_cycle(custom_cycler)
plots = ax11.plot(day_length_stations)
ax11.set_xlabel("Time (day of year)")
ax11.set_ylabel("Day length (h)")

ax11.legend(plots, station_labels)
plot_month_span(ax11)
plot_month_name(ax11, ypos=24.5)

# Show it
plt.show(block=False)

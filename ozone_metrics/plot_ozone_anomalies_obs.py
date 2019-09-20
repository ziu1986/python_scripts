try:
    data
except NameError:
    execfile("plot_ozone_observations.py")
    plt.close('all')

def anomalies(input_data, climatology_spline):
    data = input_data
    anomaly = (data.groupby(data.index.dayofyear).mean()-climatology_spline(data.index.dayofyear.unique()))
    index = data.resample("1D").mean().index
    if (anomaly.size != index.size):
        #print(anomaly.index.size, index.size)
        new_index = []
        for each in anomaly.index:
            #print(each, index[each-1])
            new_index.append(index[each-1])
        index = new_index
    return(pd.DataFrame(index=index, data=anomaly.values))

def plot_anomalies(fig, ax, data, climat, label, color):
    anomal = [anomalies(data['%d' % (iyear)], climat) for iyear in data.index.year.unique()]
    anomal = pd.concat(anomal)
    anomal[0].plot(ax=ax, label=label, color=color)
    ax.legend()

# Plot it
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()

plot_anomalies(fig1, ax11, data_jergkara, fitSpl_dmean, "Jergul/Karasjok", "orange")

# Show it
plt.show(block=False)

plt.close('all')

# Plotting
fig1 = plt.figure(1,figsize=(16,9))
fig1.canvas.set_window_title("ozone_1997_osloctm3v1")
ax11 = plt.subplot()
data_finnmark.sel(time=slice('1997-01-01','1998-01-01')).plot(ax=ax11, label='OsloCTM3v1.0')
data_jergul['O3']['1997-01-01':'1998-01-01'].plot(ax=ax11, marker='x', ls='none', color='red', alpha=0.75, label='Jergul')
data_karasjok['O3']['1997-01-01':'1998-01-01'].plot(ax=ax11, marker='+', ls='none', color='orange', alpha=0.75, label='Karasjok')
#data_svanvik['1997-01-01':'1998-01-01'].plot(marker='v', ls='none', color='blueviolet', alpha=0.75, label='Svanvik')
ax11.set_xlabel('Time')
ax11.set_ylabel('$O_3$ (ppb)')
ax11.set_ylim(0,75)
ax11.legend()

# Show it
plt.show(block=False)

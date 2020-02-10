# Plot data
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("raw_data")
ax11 = plt.subplot(projection='3d')
#ax12 = plt.subplot(212)

for each in data_list:
    selection = each.where((each['"Cond"']>0) & (each['"Photo"']>0) & (each['"Photo"']/each['"Ci"']<0.2) & (each['"Ci"']>0)).dropna()
    ax11.scatter( 1/np.sqrt(selection['"VpdL"']), selection['"Photo"']/selection['"Ci"'], selection['"Cond"'], zorder=1)
    # Plot on y-z plane
    cset = ax11.scatter(selection['"Photo"']/selection['"Ci"'], selection['"Cond"'], zs=0, zdir='x', color='black', marker='+')
    # Plot on x-z plane
    cset = ax11.scatter(1/np.sqrt(selection['"VpdL"']), selection['"Cond"'], zs=0.1, zdir='y', color='black', marker='x')
    
    
    
ax11.set_ylabel("$A_n / [CO_2]$ $(mol\;CO_2\;m^{-2}\;s^{-1})$")
ax11.set_xlabel("$D^{-0.5}$ $(kPa^{-0.5})$")
ax11.set_zlabel("$g_{sto}$ $(mol\;H_2O\;m^{-2}\;s^{-1})$", y=0)

ax11.set_ylim(0, 0.1)
#surf = ax11.plot_surface(X, Y, stomatal_conductance(xdata, 0.08, 0.3).reshape(1115,1115), cmap='binary', zorder=4)
# Add a color bar which maps values to colors.
#fig1.colorbar(surf, shrink=0.5, aspect=5)


# Show it
plt.show(block=False)

# Plot data
plt.close('all')
fig1 = plt.figure(1)#, figsize=(16,9))
fig1.canvas.set_window_title("krekling_raw_data")
ax1 = plt.subplot(projection='3d')
#ax1.set_title("Data")

fig2 = plt.figure(2)
fig2.canvas.set_window_title("krekling_raw_data_plane_fit")
ax2 = plt.subplot(projection='3d')
#ax2.set_title("Fit")

for each in data_list:
    selection = each.where((each['"Cond"']>0) & (each['"Photo"']>0) & (each['"Photo"']/each['"Ci"']<0.2) & (each['"Ci"']>0)).dropna()
    ax1.scatter( 1/np.sqrt(selection['"VpdL"']), selection['"Photo"']/selection['"Ci"'], selection['"Cond"'], zorder=1)
    # Plot on y-z plane
    cset = ax1.scatter(selection['"Photo"']/selection['"Ci"'], selection['"Cond"'], zs=0, zdir='x', color='grey', marker='+')
    # Plot on x-z plane
    cset = ax1.scatter(1/np.sqrt(selection['"VpdL"']), selection['"Cond"'], zs=0.1, zdir='y', color='grey', marker='x')
    # Plot on x-y plane
    #cset = ax1.scatter(1/np.sqrt(selection['"VpdL"']), selection['"Photo"']/selection['"Ci"'], zs=0., zdir='z', color='grey', marker='|')
    
    
    
ax1.set_ylabel("$A_n / [CO_2]$ $(mol\;CO_2\;m^{-2}\;s^{-1})$")
ax1.set_xlabel("$D^{-0.5}$ $(kPa^{-0.5})$")
ax1.set_zlabel("$g_{sto}$ $(mol\;H_2O\;m^{-2}\;s^{-1})$", y=0)

ax1.set_ylim(0, 0.1)
ax1.set_zlim(0,1)
ax1.set_xlim(0,2.5)


X, Y = np.meshgrid(flunder(x), flunder(y))
xprobe = np.vstack((X.ravel(), Y.ravel()))
surf = ax2.plot_surface(X, Y, stomatal_conductance(xprobe, fit_params[0], fit_params[1]).reshape(1115,1115), cmap='viridis', alpha=0.2) #rstride=1, cstride=1,cmap='binary', 

# Add a color bar which maps values to colors.
#fig2.colorbar(surf, shrink=0.5, aspect=5)
ax2.set_ylabel("$A_n / [CO_2]$ $(mol\;CO_2\;m^{-2}\;s^{-1})$")
ax2.set_xlabel("$D^{-0.5}$ $(kPa^{-0.5})$")
ax2.set_zlabel("$g_{sto}$ $(mol\;H_2O\;m^{-2}\;s^{-1})$", y=0)

ax2.set_xlim(ax1.get_xbound())
ax2.set_ylim(ax1.get_ybound())
ax2.set_zlim(ax1.get_zbound())

ax2.text(0,0.1,2.5, "g0: (%1.2f +- %1.2f)\ng1: (%1.2f +- %1.2f)" % (fit_params[0], cov_mat[0,0], fit_params[1], cov_mat[1,1]), size='large')
# Show it
plt.show(block=False)

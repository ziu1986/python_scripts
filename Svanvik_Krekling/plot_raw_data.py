# Plot data
plt.close('all')
fig1 = plt.figure(1, figsize=(16,9))
fig1.canvas.set_window_title("krekling_raw_data")
ax11 = plt.subplot(121, projection='3d')
ax11.set_title("Data")
ax12 = plt.subplot(122, projection='3d')
ax12.set_title("Fit")

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

X, Y = np.meshgrid(flunder(x), flunder(y))
xprobe = np.vstack((X.ravel(), Y.ravel()))
surf = ax12.plot_surface(X, Y, stomatal_conductance(xprobe, fit_params[0], fit_params[1]).reshape(1115,1115), cmap='viridis', alpha=0.2) #rstride=1, cstride=1,cmap='binary', 

# Add a color bar which maps values to colors.
#fig2.colorbar(surf, shrink=0.5, aspect=5)
ax12.set_ylabel("$A_n / [CO_2]$ $(mol\;CO_2\;m^{-2}\;s^{-1})$")
ax12.set_xlabel("$D^{-0.5}$ $(kPa^{-0.5})$")
ax12.set_zlabel("$g_{sto}$ $(mol\;H_2O\;m^{-2}\;s^{-1})$", y=0)

ax12.set_xlim(ax11.get_xbound())
ax12.set_ylim(ax11.get_ybound())
ax12.set_zlim(ax11.get_zbound())

ax12.text(0,0.1,2.5, "g0: (%1.2f +- %1.2f)\ng1: (%1.2f +- %1.2f)" % (fit_params[0], cov_mat[0,0], fit_params[1], cov_mat[1,1]), size='large')
# Show it
plt.show(block=False)

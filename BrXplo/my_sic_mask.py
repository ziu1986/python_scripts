# Set map
plt.close('all')
try:
    meridians
except NameError:
    meridians = np.arange(-180.,181.,20.)
    parallels = np.arange(-80.,81.,20.)
    ice_range = np.arange(0,1.1,0.1)
    m_nh = Basemap(projection='npstere', boundinglat=45., lon_0=0, resolution='l', round=True)
    m_sh = Basemap(projection='spstere', boundinglat=-45., lon_0=0, resolution='l', round=True)

fig1 = plt.figure(1)
fig1.canvas.set_window_title("icecov_sum_nh_%s" % year)
nh.sum(axis=(1,2)).plot(color='blue')

fig2 = plt.figure(2)
fig2.canvas.set_window_title("icecov_sum_sh_%s" % year)
sh.sum(axis=(1,2)).plot(color='blue')

fig3 = plt.figure(3)
fig3.canvas.set_window_title("my_icecov_nh_%s" % year)
data_cyclic, lon_cyclic = addcyclic(nh.where(nh.sum(axis=(1,2))==min_nh, drop=True).mean(axis=0).data, nh.lon)
# Shift grid so lon runs from -180 to +180
data_cyclic, lon_cyclic = shiftgrid(180, data_cyclic, lon_cyclic, start=False)
# Create "D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon_cyclic, nh.lat)
# Transform lat/lon into plotting coordinates for projection
x, y = m_nh(lon2d, lat2d)
cs11 = m_nh.contourf(x, y, data_cyclic, ice_range, cmap=plt.cm.OrRd)
m_nh.drawcoastlines()
m_nh.drawmapboundary()
# draw parallels and meridians.
m_nh.drawparallels(parallels)
m_nh.drawmeridians(meridians, labels=[True,False,False,True])
plt.colorbar(cs11, aspect=30)

fig4 = plt.figure(4)
fig4.canvas.set_window_title("my_icecov_sh_%s" % year)
data_cyclic, lon_cyclic = addcyclic(sh.where(sh.sum(axis=(1,2))==min_sh, drop=True).mean(axis=0).data, sh.lon)
# Shift grid so lon runs from -180 to +180
data_cyclic, lon_cyclic = shiftgrid(180, data_cyclic, lon_cyclic, start=False)
# Create "D lat/lon arrays for Basemap
lon2d, lat2d = np.meshgrid(lon_cyclic, sh.lat)
# Transform lat/lon into plotting coordinates for projection
x, y = m_sh(lon2d, lat2d)
cs12 = m_sh.contourf(x, y, data_cyclic, ice_range, cmap=plt.cm.OrRd)
m_sh.drawcoastlines()
m_sh.drawmapboundary()
# draw parallels and meridians.
m_sh.drawparallels(parallels)
m_sh.drawmeridians(meridians, labels=[True,False,False,True])
plt.colorbar(cs12, aspect=30)

plt.show(block=False)

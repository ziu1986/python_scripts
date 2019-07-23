# Clean up
plt.close('all')

# Plot it
fig1 = plt.figure(1, figsize=(16,9))
ax11 = plt.subplot()
occu_model = []
occu_data = []
for aot in np.arange(10,65,5):
    ratio = float(data_karasjok_o3.where(data_karasjok_o3>aot).dropna().size)/data_karasjok_o3.dropna().size
    occu_data.append(ratio)
    print("data: threashold%d: %1.2f" % (aot, ratio))
    ratio = float((sel_data_o3*scale.values).where(sel_data_o3*scale.values>aot, drop=True).time.size)/(sel_data_o3*scale.values).time.size
    occu_model.append(ratio)
    print("model: threashold%d: %1.2f" % (aot, ratio))

ax11.plot(np.arange(10,65,5), np.array(occu_data)*100, marker='x', color='grey', ls=':', label="Karasjok")
ax11.plot(np.arange(10,65,5), np.array(occu_model)*100, color='blue', label="OsloCTM3v1.0")
ax11.legend()
ax11.set_xlabel("$[O_3]$ threashold (ppt)")
ax11.set_ylabel("$N_{[O_3]>threashold}/N_{tot}$ (%)")

fig2 = plt.figure(2, figsize=(16,9))
ax21 = plt.subplot()

(sel_data_o3*scale.values).groupby('time.year').apply(lambda x: x.where(x>40).count()/x.count().astype(float))['O3'].plot(ax=ax21, label="OsloCTM3v1.0")
(sel_data_o3*scale.values).where((sel_data_o3.time.dt.month>=6) & (sel_data_o3.time.dt.month<9)).groupby('time.year').apply(lambda x: x.where(x>40).count()/x.count().astype(float))['O3'].plot(ax=ax21, ls='--', color='blue', label="OsloCTM3v1.0 - summer")

data_karasjok_o3['O3'].groupby(data_karasjok_o3.index.year).apply(lambda x: x.where((x.index.month>=6) & (x.index.month<9) & (x>40)).count()/x.count().astype(float)).plot(ax=ax21, color='grey', ls=':', marker='x', label="Karasjok")
data_karasjok_o3.groupby(data_karasjok_o3.index.year).apply(lambda x: x.where(x>40).count()/x.count().astype(float))['O3'].plot(ax=ax21, color='grey', ls='--', marker='v', label="Karasjok - summer")

ax21.legend()
ax21.set_xlabel("Time (year)")
#ax21.set_ylabel("Ratio")
# Show it
plt.show(block=False)

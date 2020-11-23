#from plot_results import *

plt.close('all')

test = read_data("/data/DO3SE_results/Perennial grass results_ different temps.xlsx")

fig10 = plt.figure(10, figsize=(8,12))
fig10.canvas.set_window_title("DO3SE_results_3grass")

for sheet, iax in zip(test.sheet_names[1:-1:2], np.arange(1,10)):
    date = pd.read_excel(data, sheet, header=2)
    date.index = date.index+(date['Day'].iloc[0]-1)*24
    date = date.reindex(np.arange(1,365*24))
    ax = plt.subplot(3,3,iax)

    # Plot data
    plot_pody_gsto_o3(ax, date)

# Show it
plt.show(block=False)



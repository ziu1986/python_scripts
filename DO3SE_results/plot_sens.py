from mytools.plot_tools import print_all, plot_error_bands
from do3se_tools import *

plt.close('all')

test = read_data("/data/DO3SE_results/Perennial grass results_ different temps.xlsx")

# Plot it
stop_day = 300
start_day = 90

fig10 = plt.figure(10, figsize=(16,18))
fig10.canvas.set_window_title("DO3SE_results_3grass")

for sheet, iax, ititle in zip(sorted(test.sheet_names[1:-1:2]), np.arange(1,10), char_range('a', 'i')):
    date = pd.read_excel(test, sheet, header=2)
    date.index = date.index+(date['Day'].iloc[0]-1)*24
    date = date.reindex(np.arange(1,365*24))
    ax = plt.subplot(3,3,iax)
    ax.set_title("(%s)" % ititle)

    # Plot data
    if sheet.find('2018')>=0:
        color = 'violet'
    elif sheet.find('2019')>=0:
        color = 'purple'
    else:
        color = 'blueviolet'
    plot_pody_gsto_o3(ax, date, o3color=color)

for i in np.concatenate((np.arange(3*3), np.arange(4*3-2, 6*3-1), np.arange(6*3, 9*3))):
    fig10.axes[i].set_ylabel('')

# Show it
print_all()
#plt.show(block=False)



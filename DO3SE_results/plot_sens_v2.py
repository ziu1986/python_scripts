from mytools.plot_tools import print_all, plot_error_bands
from do3se_tools import *

plt.close('all')

# Read data
src = "/data/DO3SE_results/v2/Grassland*"
try:
    test
except NameError:
    test = {}
    for each in sorted(glob.glob(src)):
        exp_name = os.path.basename(each)[:-13]
        data = read_data(each)
        test.update({exp_name: data})

pod1 = []
pod1_keys = []

# Plot it
for key in test:
    sheet_names = test[key].sheet_names[2:-1]
    for isheet in sheet_names:
        tmp = pd.read_excel(test[key], isheet, header=2)['PODY (mmol/m^2 PLA)'].max()
        pod1.append(tmp)
        pod1_keys.append((key[key.find('_')+1:], isheet))


fig1 = plt.figure(figsize=(12,8))
ax1 = plt.subplot()
ax1.set_ylabel('$POD_1$ $(mmol\,O_3\,m^{-2}\,s^{-1})$')
ax1.set_ylim(0,35)

for i in np.arange(6, len(pod1)+1, 6):
    if pod1_keys[::-1][i-6][1].find('2018') < 0:
        icolor = 'purple'
    else:
        icolor = 'violet'
        
    ax1.plot(np.arange(i-6, i)[1::2], pod1[::-1][i-6:i][::2], marker='o', fillstyle='none', color=icolor)
    ax1.plot(np.arange(i-6, i)[::2], pod1[::-1][i-6:i][1::2], marker='o', color=icolor)

plt.show(block=False)

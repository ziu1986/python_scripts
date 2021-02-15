from mytools.plot_tools import print_all, plot_error_bands
from do3se_tools import *
from mytools.ozone_tools import flunder

plt.close('all')

def write_pod1(pod1, pod1_keys, species_name):
    outfile = open('pod1_max_%s.dat' % species_name, 'w')

    outfile.write("Type,Year,PPFD,SWP,POD1\n")
    for ikey, ival in zip(pod1_keys, pod1):
        split = ikey[1].split('_')
        if len(split) < 3:
            split.append('0')
        outfile.write("%s,%s,%s,%s,%s\n" % (ikey[0], split[0], split[1], split[2].replace('SWP', '1'), ival))
    outfile.close()


def plot_pod(results, species_name):
    fig1 = plt.figure(figsize=(12,8))
    fig1.canvas.set_window_title('do3se_results_pod_%s' % species_name)
    ax1 = plt.subplot()
    ax1.set_ylabel('$POD_1$ $(mmol\,O_3\,m^{-2})$')
    ax1.set_ylim(0,35)
    
    for ispecies in ('Boreal', 'Cold', 'MM'):
            
        results['POD1'].where((results['Year']==2018) &
                              (results['SWP']==0) &
                              (results['Type']==ispecies)).dropna().plot(color='violet', marker='o', markersize=9, label='_')
        results['POD1'].where((results['Year']==2019) &
                              (results['SWP']==0) &
                              (results['Type']==ispecies)).dropna().plot(color='purple', marker='o', markersize=9, label='_')

        results['POD1'].where((results['Year']==2018) &
                              (results['SWP']==1) &
                              (results['Type']==ispecies)).dropna().plot(color='violet', fillstyle='none', marker='o', markersize=9, label='_')
        results['POD1'].where((results['Year']==2019) &
                              (results['SWP']==1) &
                              (results['Type']==ispecies)).dropna().plot(color='purple', fillstyle='none', marker='o', markersize=9, label='_')

    # Annotate PPFD
    for idx, iPPDF, ipod in zip(results.where((results['SWP']==0)).index, results['PPFD'].where((results['SWP']==0)), results['POD1'].where((results['SWP']==0))):
        
        if (iPPDF == 'PPFD0.8'):
            annotation = '$-$'
            y = ipod+1
        elif (iPPDF == 'PPFD1.2'):
            annotation = '$+$'
            y = ipod-1.5
        else:
            annotation = ''
            y = 0
            
        ax1.text(idx, y, annotation, size='x-large')
        
    ax1.axvspan(-0.5, 11.5, color='linen')
    ax1.axvspan(23.5, 35.5, color='linen')
    ax1.set_xlim(-0.5,35.5)
    ax1.set_xticks(np.arange(6,36,12))
    ax1.set_xticklabels(('Boreal', 'Cold', 'MM'), size='xx-large')

    # Fake legend
    ax1.plot(np.arange(-1,-3,-1), color='violet', ls='-', label='2018')
    ax1.plot(np.arange(-1,-3,-1), color='purple', ls='--', label='2019')
    
    ax1.plot(np.arange(-1,-3,-1), color='black', marker='o', ls="None", label='SWP=off')
    ax1.plot(np.arange(-1,-3,-1), color='black', marker='o', fillstyle='none', ls="None", label='SWP=on')
        
    ax1.legend(loc='upper right', ncol=2)
    #for i in np.arange(6, len(pod1)+1, 6):
    #    param_type = pod1_keys[::-1][i-6][0]
    #            
    #    if pod1_keys[::-1][i-6][1].find('2018') < 0:
    #        icolor = 'purple'
    #        x_range = np.arange(i, i+6)
    #    else:
    #        icolor = 'violet'
    #        x_range = np.arange(i-12, i-6)
    #        ax1.text(i-9, 34, "%s" % param_type, size='x-large')
    #    
    #    ax1.plot(x_range[1::2], pod1[::-1][i-6:i][::2],
    #                       marker='o', fillstyle='none', color=icolor,
    #                       label='f_SW var.')
    #    ax1.plot(x_range[::2], pod1[::-1][i-6:i][1::2],
    #             marker='o', color=icolor, label='f_SW=1')

    
# Read data
species = ("Birch", "Spruce", "Grassland")
for PFT in species:
    src = "/data/DO3SE_results/v2/%s*" % PFT

    if not os.path.isfile('pod1_max_%s.dat' % PFT):
        pod1 = []
        pod1_keys = []
        for each in sorted(glob.glob(src)):
            exp_name = os.path.basename(each)[:-13]
            data = read_data(each)
            sheet_names = data.sheet_names[2:-1]
            for isheet in sheet_names:
                pod1_max = pd.read_excel(data, isheet, header=2)['PODY (mmol/m^2 PLA)'].max()
                pod1.append(pod1_max)
                pod1_keys.append((exp_name[exp_name.find('_')+1:].replace(' temp', ''), isheet.replace('_output', '').replace('-20%', 'PPFD0.8').replace('+20%', 'PPFD1.2')))

        write_pod1(pod1, pod1_keys, PFT)
    else:
        results = pd.read_csv('pod1_max_%s.dat' % PFT)

        plot_pod(results, PFT)
        # Save it
        print_all()
        # Show it
        plt.show(block=False)


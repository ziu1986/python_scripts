from mytools.plot_tools import print_all, plot_error_bands
from do3se_tools import *
from mytools.ozone_tools import flunder

plt.close('all')

def write_pod1(pod1, pod1_keys, species_name, version):
    outfile = open('%s_pod1_max_%s.dat' % (version, species_name), 'w')

    outfile.write("Type,Year,PPFD,SWP,POD1\n")
    for ikey, ival in zip(pod1_keys, pod1):
        split = ikey[1].split('_')
        if len(split) < 3:
            split.append('0')
        outfile.write("%s,%s,%s,%s,%s\n" % (ikey[0], split[0], split[1], split[2].replace('SWP', '1'), ival))
    outfile.close()


def plot_pod(results, species_name, **karg):
    crit_level = karg.pop('CL', -1)
    annotation = karg.pop('annote', None)

    # Annotate PPFD
    def print_annotation(ax, date, **karg):      
        size = karg.pop('size', 'xx-large')

        for idx, iPPDF, ipod in zip(date.where((date['SWP']==0)).index, date['PPFD'].where((date['SWP']==0)), date['POD1'].where((date['SWP']==0))):
        
            if (iPPDF == 'PPFD0.8'):
                annotation = '$-$'
                y = ipod+1
            elif (iPPDF == 'PPFD1.2'):
                annotation = '$+$'
                y = ipod-1.5
            else:
                annotation = ''
                y = 0
            
            ax.text(idx, y, annotation, size=size)

    if 'ref' in karg:
        b_ref = True
        alpha = 0.55
    else:
        b_ref = False
        
    reference = karg.pop('ref', None)
    
    fig1 = plt.figure(figsize=(12,8))
    if annotation:
        fig1.canvas.set_window_title('do3se_results_pod_%s_%s' % (annotation, species_name))
    else:
        fig1.canvas.set_window_title('do3se_results_pod_%s' % species_name)
    ax1 = plt.subplot()
    ax1.set_ylabel('$POD_1$ $(mmol\,O_3\,m^{-2})$')
    ax1.set_ylim(0,35)
    ax1.axhline(crit_level, color='grey', ls='--', linewidth=3, alpha=0.5)
    
    for ispecies in ('MM', 'Cold', 'Boreal'):
        if b_ref:
            reference['POD1'].where((reference['Year']==2018) &
                              (reference['SWP']==0) &
                              (reference['Type']==ispecies)).dropna().plot(color='violet', ls=':', marker='o', markersize=9, label='_', alpha=alpha)
            reference['POD1'].where((reference['Year']==2019) &
                              (reference['SWP']==0) &
                              (reference['Type']==ispecies)).dropna().plot(color='purple', ls=':', marker='o', markersize=9, label='_', alpha=alpha)

            reference['POD1'].where((reference['Year']==2018) &
                              (reference['SWP']==1) &
                              (reference['Type']==ispecies)).dropna().plot(color='violet', ls=':', fillstyle='none', marker='o', markersize=9, label='_', alpha=alpha)
            reference['POD1'].where((reference['Year']==2019) &
                              (reference['SWP']==1) &
                              (reference['Type']==ispecies)).dropna().plot(color='purple', ls=':', fillstyle='none', marker='o', markersize=9, label='_', alpha=alpha)

            
        results['POD1'].where((results['Year']==2018) &
                              (results['SWP']==0) &
                              (results['Type']==ispecies)).dropna().plot(color='violet', ls='-', marker='o', markersize=9, label='_')
        results['POD1'].where((results['Year']==2019) &
                              (results['SWP']==0) &
                              (results['Type']==ispecies)).dropna().plot(color='purple', ls='--', marker='o', markersize=9, label='_')

        results['POD1'].where((results['Year']==2018) &
                              (results['SWP']==1) &
                              (results['Type']==ispecies)).dropna().plot(color='violet', ls='-', fillstyle='none', marker='o', markersize=9, label='_')
        results['POD1'].where((results['Year']==2019) &
                              (results['SWP']==1) &
                              (results['Type']==ispecies)).dropna().plot(color='purple', ls='--', fillstyle='none', marker='o', markersize=9, label='_')

    # Annotate PPFD
    print_annotation(ax1, results)
    if b_ref:
        print_annotation(ax1, reference, size='small')

        
        
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
    
    ax1.plot(np.arange(-1,-3,-1), color='black', ls=':', alpha=alpha, label='GS MM')
    ax1.plot(np.arange(-1,-3,-1), color='black', ls='-', alpha=alpha, label='GS besp.')
        
    ax1.legend(loc='upper right', ncol=3)

   
        
# Read data
critical_level = {'Birch': 5.2, 'Spruce': 9.2, 'Grassland': 10.2}
species = ("Birch", "Spruce", "Grassland")
version = 'v3'
reference = 'v2'

for PFT in species:
    src = "/data/DO3SE_results/%s/*%s*" % (version, PFT)

    if not os.path.isfile('%s_pod1_max_%s.dat' % (version, PFT)):
        print('Read data')
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
        # Save to file for easier import into pandas
        write_pod1(pod1, pod1_keys, PFT, version)

    # Read data from file    
    results = pd.read_csv('%s_pod1_max_%s.dat' % (version, PFT))
    if os.path.isfile('%s_pod1_max_%s.dat' % (reference, PFT)):
        ref_results = pd.read_csv('%s_pod1_max_%s.dat' % (reference, PFT))
    # Plot it
    plot_pod(results, PFT, CL=critical_level[PFT], annote='%s' % version, ref=ref_results)
    
# Save it
print_all()
# Show it
plt.show(block=False)


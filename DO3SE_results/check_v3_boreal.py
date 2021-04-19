from mytools.plot_tools import print_all, plot_error_bands, get_month_name
from do3se_tools import *

plt.close('all')

species = ("Birch", "Spruce", "Grassland")
version = 'v2'
for PFT in species:
    src = "/data/DO3SE_results/%s/*%s*" % (version, PFT)

    data = read_data(sorted(glob.glob(src))[0])
    test = {}
    for isheet, iyear in zip((3,9),(2018, 2019)):
        test.update({iyear:pd.read_excel(data, data.sheet_names[isheet], header=2)})
        test[iyear].index = test[iyear].index+(test[iyear]['Day'].iloc[0]-1)*24

    indi = pd.date_range('2018-01-01', '2019-01-01')
    # Plot it
    fig1 = plt.figure(figsize=(10,12))
    fig1.canvas.set_window_title("check_v3_boreal_18-19_%s_%s" % (version, PFT))
    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)
    #(test[2018][u'f_temp']-test[2019][u'f_temp']).plot(ax=ax1, color='blue')
    (test[2018][u'f_SW']-test[2019][u'f_SW']).plot(ax=ax1, color='blue')
    #(test[2018][u'f_VPD']-test[2019][u'f_VPD']).plot(ax=ax2, color='red')
    (test[2018][u'f_phen']-test[2019][u'f_phen']).plot(ax=ax2, color='red')
    (test[2018][u'f_light']-test[2019][u'f_light']).plot(ax=ax3, color='black')
    (test[2018][u'PODY (mmol/m^2 PLA)']-test[2019][u'PODY (mmol/m^2 PLA)']).plot(ax=ax4, color='darkgreen')

    for ax in fig1.axes:
        ax.set_ylim(-1,1)
        ax.set_xticklabels("")
    ax4.set_xticklabels(np.round(ax.get_xticks()/24))
    for each in np.round(ax1.get_xticks()[1:-1]/24):
        ax1.text(each*24, ax.get_ylim()[-1], "%s" % get_month_name(indi.groupby(indi.dayofyear)[each].month[0], length=3), size='x-large')
    ax4.set_ylim(-4.5, 4.5)
    ax4.set_xlabel("Time (doy)")
    #ax1.set_ylabel('$\Delta f_{T}$')
    ax1.set_ylabel('$\Delta f_{SWP}$')
    #ax2.set_ylabel('$\Delta f_{VPD}$')
    ax2.set_ylabel('$\Delta f_{phen}$')
    ax3.set_ylabel('$\Delta f_{light}$')
    ax4.set_ylabel('$\Delta POD_y$ ($mmol\,m^{-2}$ PLA)')


# Show it
plt.show(block=False)
# Print it
print_all()


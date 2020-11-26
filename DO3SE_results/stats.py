import statsmodels.api as sm

def plot_stats(date, uncum_date, key_list, **karg):
    color = karg.pop('color', 'black')
    fit_color = karg.pop('fit_color', 'red')
    file_name = karg.pop('file_name', None)

    if file_name:
        fig_path = 'fit_results'
        if not os.path.isdir(fig_path):
            os.mkdir(fig_path)
        file = open(fig_path+'/'+file_name, 'w')
        
    for key, iax in zip(sorted(key_list), np.arange(len(key_list))):
        ax = plt.subplot(4,3,iax+1)
        ax.set_title("%s" % key, y=0.85, x=0.22)

        x = date.where(uncum_date['PODY']!=0)[key]
        y = uncum_date['PODY'].where(uncum_date['PODY']!=0)*1e3

        # Plot it
        ax.plot(x, y, marker='x', ls='none', color=color)
        
        if(x.dropna().size!=y.dropna().size):
            x = x[:-1]

        # Compute linear regression
        # Need to add ones to get ordinate intersception
        x = sm.add_constant(x)
        try:
            model = sm.OLS(y.dropna(), x.dropna())
        except ValueError:
            print(x.dropna().size, y.dropna().size)
            return
        results = model.fit()
        test_x = np.linspace(x.min()[key], x.max()[key])
        if (results.params.size>1):
            ax.plot(test_x, results.params[0]+test_x*results.params[1], ls='--', color=fit_color)
            ax.text(0.9, 0.9,"y=%2.3f x+%2.3f" % (results.params[1], results.params[0]), horizontalalignment='right', transform=ax.transAxes)
            if file_name:
                file.write(results.summary2().as_text())
            else:
                print(results.summary2())
    if file:
        file.close()

        
plt.close('all')
"""
for spec, figi  in zip(species, np.arange(3,6)):
    fig = plt.figure(figi, figsize=(16,6))
    fig.canvas.set_window_title("DO3SE_results_pody_gsto_o3_%s" % spec.replace(' ', '_'))

    data = data_list[spec]
    for iax, sheet, color, ititle in zip(np.arange(1,10), data.sheet_names[1::2][:3], ('violet', 'purple', 'blueviolet'), char_range('a', 'c')):
        ax = plt.subplot(1,3,iax)
        ax.set_title("(%s)" % ititle)

        date = pd.read_excel(data, sheet, header=2)
        date.index = date.index+(date['Day'].iloc[0]-1)*24
        date = date.reindex(np.arange(1,365*24))
        # Plot data
        plot_pody_gsto_o3(ax, date, o3color=color)

    for ax in fig.axes[1:-2]:
        ax.set_ylabel("")
        ax.set_yticklabels("")
"""

for spec in species:
    data = data_list[spec]
    for sheet, color, ifig in zip(data.sheet_names[1::2][:3], ('violet','purple','blueviolet'), np.arange(7,10)):
        fig = plt.figure(ifig, figsize=(16,10))
        fig.canvas.set_window_title("DO3SE_results_%s" % sheet.replace(' ', '_'))
        date = pd.read_excel(data, sheet, header=2)
        uncum_date = pd.DataFrame({'PODY':uncumsum(date,'PODY (mmol/m^2 PLA)'), 'Month':date['Month'], 'Doy':date['Day']})
        plot_stats(date, uncum_date, key_list, color=color, file_name='fit_results_%s_%s.csv' % (sheet.replace(' ', '_'), 'clim'))
        for ax in fig.axes:
            ax.set_ylim(0, np.round(uncum_date['PODY'].max()*1e3))
        fig.axes[6].set_xlabel("PODY ($10^{-6} mol m^{-2}$)")
    print_all()
    plt.close('all')
# Show it
#plt.show(block=False)

import statsmodels.api as sm

def plot_stats(date, uncum_date, key_list, **karg):
    color = karg.pop('color', 'black')
    fit_color = karg.pop('fit_color', 'red')
    for key, iax in zip(sorted(key_list), np.arange(len(key_list))):
        ax = plt.subplot(4,3,iax+1)
        ax.set_title("%s" % key, y=0.85, x=0.22)

        x = date.where(uncum_date['PODY']!=0)[key]
        y = uncum_date['PODY'].where(uncum_date['PODY']!=0)

        # Plot it
        ax.plot(x, y, marker='x', ls='none', color=color)
        # Compute linear regression
        # Need to add ones to get ordinate intersception
        x = sm.add_constant(x)
        if(x.dropna().size!=y.dropna().size):
            x = x[:-1]
        
        model = sm.OLS(y.dropna(), x.dropna())
        results = model.fit()
        test_x = np.linspace(x.min()[key], x.max()[key])
        if (results.params.size>1):
            ax.plot(test_x, results.params[0]+test_x*results.params[1], ls='--', color=fit_color)
        print(results.summary2())

plt.close('all')

fig7 = plt.figure(7, figsize=(16,10))
plot_stats(date_clim, uncum_date_clim, key_list, color='blueviolet')

fig8 = plt.figure(8, figsize=(16,10))
plot_stats(date, uncum_date, key_list, color='purple')

    
# Show it
plt.show(block=False)

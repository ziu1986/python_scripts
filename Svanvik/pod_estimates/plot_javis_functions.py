from func_def_plot_javis import *

# Data sources
src_svanvik = os.environ['DATA']+'/astra/observations/metdata_svanvik/Svanvik_temp_relhum_wind_*.csv'
src_svanvik_rad = os.environ['DATA']+'/astra/observations/metdata_svanvik/svanvik_glob_rad_*.csv'

# main
alpha_var_decr_20 = {'coniferous':0.0075, 'deciduous':0.00525, 'grassland':0.01375}
alpha_var_incr_20 = {'coniferous':0.005, 'deciduous':0.0035, 'grassland':0.0092}

# Set up the different species
# Coniferous
coniferous = JavisModel('coniferous', Tmin=0, Tmax=200, Topt=20, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
coniferous_subarctic = JavisModel('coniferous_subarctic', Tmin=0, Tmax=200, Topt=10, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
coniferous_cold = JavisModel('coniferous_cold', Tmin=0, Tmax=200, Topt=15, fmin=0.1, Dmax=0.8, Dmin=2.8, alpha=0.006, gmax=125)
# Deciduous
deciduous = JavisModel('deciduous', Tmin=5, Tmax=200, Topt=20, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
deciduous_subarctic = JavisModel('deciduous_subarctic', Tmin=0, Tmax=200, Topt=10, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
deciduous_cold = JavisModel('deciduous_cold', Tmin=5, Tmax=200, Topt=15, fmin=0.1, Dmax=0.5, Dmin=2.7, alpha=0.0042, gmax=240)
# Grassland
grassland = JavisModel('grassland', Tmin=10, Tmax=36, Topt=24, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)
grassland_subarctic = JavisModel('grassland_subarctic', Tmin=0, Tmax=24, Topt=10, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)
grassland_cold = JavisModel('grassland_cold', Tmin=5, Tmax=36, Topt=15, fmin=0.1, Dmax=1.75, Dmin=4.5, alpha=0.011, gmax=190)


# Load data
# Convert Wm^-2 to mumolm^-2s^-1
conv_rad = 2.77e18 / 6.022141e23 * 1e6
data_temp = import_data(src_svanvik)
data_rad = import_data(src_svanvik_rad)*conv_rad
# Comput vpd
vpd = VPD(data_temp.iloc[:,1], data_temp.iloc[:,0])/kilo

# Clean up
plt.close('all')

# Loop through species and plot them
#for species, i in zip((coniferous, deciduous, grassland, coniferous_subarctic, deciduous_subarctic, grassland_subarctic, coniferous_cold, deciduous_cold, grassland_cold), np.arange(1,10)):
#    plot_f_functions(species, i)

#print_all()
#plt.close('all')

for species in ((coniferous, coniferous_subarctic, coniferous_cold), (deciduous, deciduous_subarctic, deciduous_cold), (grassland, grassland_subarctic, grassland_cold)):
    plot_temperature_histogram(data_temp.iloc[:,0], javis_model=species)
   
for species in ([coniferous,], [deciduous,], [grassland,]):
    for kind in ('coniferous', 'deciduous', 'grassland'):
        if kind in species[0].name:
            for alpha_var, alpha_mode in zip((alpha_var_incr_20[kind], alpha_var_decr_20[kind]), ('PPFD1.2', 'PPFD0.8')):
                temp = copy.deepcopy(species[0])
                temp.name = temp.name + "_" + alpha_mode
                temp.alpha_light = alpha_var
                species.append(temp)
    plot_histogram(data_rad.iloc[:,0], javis_model=species, var='light')

for species in ([coniferous,], [deciduous,], [grassland,]):
    plot_histogram(vpd, javis_model=species, var='vpd')

#print_all()
#plt.close('all')
#for species in (coniferous, coniferous_subarctic, coniferous_cold, deciduous, deciduous_subarctic, deciduous_cold, grassland, grassland_subarctic, grassland_cold):
#    print(species.name)
#    result = get_f_function(species)
#    temp = result[-1].where((result[-1].index.month>=5)&(result[-1].index.month<9)).resample('1Y').sum()['2000':]
#    temp_mean = temp.mean()
#    temp_std = temp.std()
#    print(temp.apply(lambda x: (x-temp_mean)/temp_std))


"""
list = []
err_list = []
midnight_sun_list = []
# Loop through species and print variance
for species in (deciduous_subarctic, deciduous_cold, deciduous, coniferous_subarctic, coniferous_cold, coniferous, grassland_subarctic, grassland_cold, grassland):
    for kind in ( 'deciduous', 'coniferous','grassland'):
        if kind in species.name:
            print(species.name)
            result = get_f_function(species, data_temp, vpd, data_rad)
            list.append(get_stats(result[-1], type='noon', stats='mean'))
            list.append(get_stats(result[-1], type='morning', stats='mean'))
            err_list.append(get_stats(result[-1], type='noon', stats='std'))
            err_list.append(get_stats(result[-1], type='morning', stats='std'))
            test_res =result[2].where(result[2].index.hour==0)
            midnight_sun_list.append(test_res.groupby([test_res.index.month,test_res.index.day]))

            tmp = copy.deepcopy(species)
            
            tmp.name = tmp.name[:-3] + "-20"
            tmp.alpha_light = alpha_var_decr_20[kind]
            result = get_f_function(tmp, data_temp, vpd, data_rad)
            print(tmp.name)
            list.append(get_stats(result[-1], type='noon', stats='mean'))
            list.append(get_stats(result[-1], type='morning', stats='mean'))
            err_list.append(get_stats(result[-1], type='noon', stats='std'))
            err_list.append(get_stats(result[-1], type='morning', stats='std'))
            test_res =result[2].where(result[2].index.hour==0)
            midnight_sun_list.append(test_res.groupby([test_res.index.month,test_res.index.day]))
           
            tmp.name = tmp.name + "+20"
            tmp.alpha_light = alpha_var_incr_20[kind]
            result = get_f_function(tmp, data_temp, vpd, data_rad)
            print(tmp.name)
            list.append(get_stats(result[-1], type='noon', stats='mean'))
            list.append(get_stats(result[-1], type='morning', stats='mean'))
            err_list.append(get_stats(result[-1], type='noon', stats='std'))
            err_list.append(get_stats(result[-1], type='morning', stats='std'))
            test_res =result[2].where(result[2].index.hour==0)
            midnight_sun_list.append(test_res.groupby([test_res.index.month,test_res.index.day]))
            
    
plot_optimal(list, err=err_list, stats='mean')   
"""
# Show it
plt.show(block=False)

print_all()


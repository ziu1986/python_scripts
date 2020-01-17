import pandas as pd
import numpy as np

#creating some dummy data
n_years = 3
n_hours = 3
st_times = ['01-01-198{0} 00:00'.format(i) for i in range(n_years)]
nd_times = ['01-01-198{0} 0{1}:00'.format(i,n_hours-1) for i in range(n_years)]

indx_list = []
for s, e in zip(st_times, nd_times):
    indx = pd.date_range(start=s, end=e, freq='H')
    indx_list.append(indx.values)
index = pd.DatetimeIndex(np.concatenate(indx_list,axis=0))

data = pd.DataFrame({'rainfall': list(range(n_years*n_hours)),
              'rainfall_1': list(reversed(range(n_years*n_hours)))
             }, index=index)

#creating the hour, day, month, & day columns
data.loc[:,'hour'] = data.index.hour.values
data.loc[:,'day'] = data.index.day.values
data.loc[:,'month'] = data.index.month.values

#create groups and calculate the mean of each group
print(data.groupby(['month','day','hour']).mean())

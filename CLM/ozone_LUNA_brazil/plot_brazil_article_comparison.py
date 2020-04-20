import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

plt.close('all')

fig_clm = pkl.load(open('brazil_test_vcmax_jmax_ratio.pkl','rb'))
#fig_meta = pkl.load(open('../../plant_model/ozone_response_vcmax_jmax_ratio.pkl','rb'))

# Show it
plt.show(block=False)


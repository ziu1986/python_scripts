import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

plt.close('all')

fig_clm = pkl.load(open('brazil_test_vcmax_jmax_ratio.pkl','rb'))
# Python2->python3 comp issue!
with open('../../plant_model/ozone_response_vcmax_jmax_ratio.pkl','rb') as f:
    fig_meta = pkl.load(f, encoding='bytes')

# Show it
plt.show(block=False)


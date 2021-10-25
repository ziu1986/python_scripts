import random
import numpy as np
import pandas as pd

n = 22
kategories = 22
probe = 10000
penalty = 4*n

test_src = pd.read_csv("preferences.txt")
test = [test_src['Pref%d' % i].to_numpy() for i in range(1,4)]


if n > len(test[0]):
    n = len(test[0])
    print("Change sample size: %d" % n)


for i in range(probe):
    weight = np.repeat(4, n)
    assign = np.array(random.sample(range(1,kategories+1), k=n))
    for i, itest in zip(np.arange(3,0,-1), test[::-1]):
        weight[np.where(itest==assign)] = i
    if weight.sum() < penalty:
        penalty = weight.sum()
        result = assign

if penalty <= 3*n:
    print("OK")
    
print(penalty, result)


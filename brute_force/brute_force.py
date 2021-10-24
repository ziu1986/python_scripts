import random
import numpy as np 

n = 22
probe = 1000000
penalty = 4*n

test = [np.array(random.sample(range(1,n+1), k=n)) for i in range (3)]


for i in range(probe):
    weight = np.repeat(4, n)
    assign = np.array(random.sample(range(1,n+1), k=n))
    for i, itest in zip(np.arange(3,0,-1), test[::-1]):
        weight[np.where(itest==assign)] = i
    if weight.sum() < penalty:
        penalty = weight.sum()
        result = assign

if penalty <= 3*n:
    print("OK")
    
print(penalty, result)


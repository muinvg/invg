# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Inverse G

# <codecell>

import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

# <codecell>

a = 3**3
print(a)

# <codecell>

x=np.linspace(0,5,21)

# <codecell>

print(x)

# <codecell>

plt.plot(x,np.sin(x))

# <codecell>

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# <codecell>

def func(x, t, k, m):
    return [x[2]/m, 
            x[3]/m, 
            -x[0]*x[1]*x[1]*k,
            -x[0]*x[0]*x[1]]

# <codecell>

m=1
k=1
x0 = [0.0, 0.0, -1, -0.9]
t = np.arange(0, 100, 0.01)

# <codecell>

x = odeint(func, x0, t, args=(k,m))

# <codecell>

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:, 0], x[:, 1], x[:, 2])

# <codecell>

h = np.column_stack([x[:,0],x[:,1]])
h.shape

# <codecell>

U = x[:,0]*x[:,0]*x[:,1]*x[:,1]/2 
K =  (x[:,2]*x[:,2] + x[:,3]*x[:,3] )/(2*m)
H = U + K
plt.plot(t,np.column_stack([H,U,K]))

# <codecell>

plt.plot(x[:,0],x[:,1])

# <codecell>


# <codecell>



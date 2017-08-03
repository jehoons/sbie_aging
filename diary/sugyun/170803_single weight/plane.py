
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-3.9, 2.1, 0.1)
Y = np.arange(-3.9, 2.1, 0.1)
# print(Y,'ddd')
X, Y = np.meshgrid(X, Y)
# print(X.shape)
# # print(np.meshgrid(X, Y))
# print(len(Y))
R = np.sqrt(X**2 + Y**2)
# print(R)
list= []
f= open('output.txt','r')
for i,line in enumerate(f):
    line = line.strip().split('\t')
    list.append(float(line[2]))
    # print(line[2])
Z = np.array(list)
Z =Z.reshape(60,60)
print(X.shape,Y.shape,Z.shape)
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
# ax.set_zlim(-0.01, 50.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
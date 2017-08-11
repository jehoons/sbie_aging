import numpy as np
from sklearn.manifold import TSNE
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import numpy as np
import matplotlib.pyplot as plt



f=open('output_all2.txt','r')
X = []
C = []
# s = np.random.choice(range(5**9),10000)
# print(s)
for n, line in enumerate(f):
    if n > 10000:
        continue
    line = line.strip().split('\t')
    x = [float(x) for x in line]
    # if x[-1] == 0:
    #     print(line)
    #     X.append(x[0:-1])
    #     C.append(x[-1])
    # elif n not in s:
    #     continue
    # x=[float(x) for x in line]
    X.append(x[0:-1])
    C.append(x[-1])
X = np.array(X)
model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100)
np.set_printoptions(suppress=True)
X2 =model.fit_transform(X)
# print(X2)
x, y = X2.T

# colors = cm.rainbow(np.linspace(0, 1, len(C)))
plt.scatter(x, y, c=C, cmap=cm.bwr)
# plt.scatter(x,y,c=C)

plt.show()
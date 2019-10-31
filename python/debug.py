import h5py
import matplotlib.pyplot as plt
import numpy as np
import sys

file = h5py.File(sys.argv[1], 'r')
path = sys.argv[2].split('/')
op = sys.argv[3]

dataset = file
for group in path:
    dataset = dataset[group]
y = np.array(dataset)

if op == "real":
    y = np.real(y)
if op == "imag":
    y = np.imag(y)
if op == "angle":
    y = np.angle(y)
if op == "norm2":
    y = np.abs(y)**2

fig,ax = plt.subplots()
ax.plot(y,label=sys.argv[2])
ax.set(xlim=[0, y.shape[0]])

# x = x_name = None
# y_name = sys.argv[-1]
# y = np.array(file["Debug"][y_name])

# if(len(sys.argv) == 4):
#     x_name = sys.argv[-2]
#     x = np.array(file["Debug"][x_name])

# print(y)

# fig,ax = plt.subplots()
# if x is None:
#     ax.plot([i for i in range(len(y))], y)
#     ax.set(ylabel=y_name)
# else:
#     ax.plot(x, y)
#     ax.set(xlabel=x_name, ylabel=y_name)

ax.legend(*ax.get_legend_handles_labels())
plt.grid(which="both")
plt.show()






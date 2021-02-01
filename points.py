from matplotlib import pyplot as plt

# ps = [(0,5), (3,4), (5,4), (2,3), (4,2), (1, 2), (4, 1), (5, 0)]
# plt.scatter([p[0] for p in ps], [p[1] for p in ps])
# plt.show()
xs = [4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 0, 0, 0]
ys = [0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 4, 3, 2]

plt.scatter(xs, ys)
plt.show()

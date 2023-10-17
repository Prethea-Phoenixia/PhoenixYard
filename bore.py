import csv

with open("file.csv", newline="\n") as csvfile:
    reader = csv.DictReader(csvfile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
    data = [line for line in reader]


print(data)

import matplotlib.pyplot as plt


ax = plt.subplot()

for line in data:
    Pb = float(line["Breech Pres."])
    Ps = float(line["Shot Pres."])
    x = float(line["Travel"])

    ax.scatter([x, 0], [Ps, Pb], c=["r", "k"])

    xs = []
    ys = []
    for i in range(100):
        k = i / 100
        xs.append(x * k)
        ys.append(Pb + (Ps - Pb) * (k) ** 2)

    ax.plot(xs, ys, c="k", ls="-", alpha=0.5)

ax.set_ylim(0, None)
ax.set_xlim(0, None)
plt.show()

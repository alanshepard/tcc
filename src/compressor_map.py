import matplotlib.pyplot as plt

from compressor import Compressor

fig, ax = plt.subplots(1,1, figsize=(4,3))
c=Compressor()
c.plot_map(ax, samples=144)

plt.savefig('compressor_map.pdf')


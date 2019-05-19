import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from compressor import Compressor

plt.figure(figsize=(6,8),dpi=150)
Compressor().plot_map(plt.gca())
plt.savefig('map.pdf', dpi=300)

plt.show()

import matplotlib.pyplot as plt

from turbine import TurbineExtendedMap

fig, ax = plt.subplots(1,1, figsize=(4,3))
t=TurbineExtendedMap()
t.geom['A_ratio']=4
t.geom['R_ratio']=2
t.plot_map(ax, extent=(0,0.55,1,3), choke_line=False, samples=144)

plt.savefig('turbine_map_typical.pdf')


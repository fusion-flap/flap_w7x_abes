import os
import flap
from matplotlib import pyplot as plt

folder = "/DATA/repos/flap/modules/flap_spade/recon"
filename = '20180823.011_light_orig_1696426019.hdf5'
file = os.path.join(folder, filename)
data = flap.load(file)
for light in data.data:
    plt.plot(light)
    plt.show()
    plt.pause(1)
    plt.clf()
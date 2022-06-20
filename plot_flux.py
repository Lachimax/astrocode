# interactive plot for a spectral cube
# Author: C. Derkenne. Date: 3/2019. 


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot_flux(cube_in):
""" Takes in a 3D data cube, with 0th axis wavelength.
Returns the chosen x and y dimension values in pixel units """

	flatcube = np.sum(cube_in, axis = 0)
	b, a = flatcube.shape

	fig, ax = plt.subplots()
	plt.subplots_adjust(bottom=0.25)

	l = ax.imshow(flatcube, interpolation=None, origin='lower', cmap='gray')

	axcolor = 'lightgoldenrodyellow'
	axcropx = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
	axcropy = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

	slcropx = Slider(axcropx, 'Crop X (px)', 0, a, valinit=a, valstep=2)
	slcropy = Slider(axcropy, 'Crop Y (px)', 0, b, valinit=b, valstep=2)

	def update(val):
	    x = slcropx.val
	    y = slcropy.val
	    l = ax.imshow(flatcube[int(b/2 - y/2):int(b/2 + y/2),int(a/2 - x/2):int(a/2 + x/2)],
	    			  interpolation=None, origin='lower', cmap='gray')
	    fig.canvas.draw()

	slcropx.on_changed(update)
	slcropy.on_changed(update)

	plt.show()

	return slcropx.val, slcropy.val
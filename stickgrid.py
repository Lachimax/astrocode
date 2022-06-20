from random import random
import numpy as np
import matplotlib.pyplot as plt

def stickOnLine(x, y, length, angle):
    """
    A stick is thrown onto a uniform grid with gridlines spaced 1 unit apart.
    This function determines if the line crosses a gridline.
    """
    x2 = abs((x + l*np.cos(angle))//1)
    y2 = abs((y + l*np.sin(angle))//1)
    return x2+y2==1


if __name__=='__main__':
	iterations=100
	l_vals = []
	p_vals = []
	
	# simulate line lengths from 0 to length of the diagonal of the grid squares
	for l in np.arange(0, 5**0.5, 0.02):
	    c=0
	    for i in range(iterations):
	        if (stickOnLine(random(), random(), l, random()*np.pi**2)):
                    c+=1
	    print(l, float(c)/iterations)
	    l_vals.append(l)
	    p_vals.append(float(c)/iterations)
	    
	    
	plt.plot(l_vals, p_vals)
	plt.title('Chance of stick of length L landing on a gridline')
	plt.xlabel('length')
	plt.ylabel('probability')
	plt.show()

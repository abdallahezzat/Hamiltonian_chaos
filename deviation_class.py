import pickle
import numpy as np
import matplotlib.pyplot as plt

class deviation():
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Constructor ::

	def __init__(self, x1, y1, x2, y2):
		self.x1 = x1
		self.y1 = y1
		self.x2 = x2
		self.y2 = y2

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Functions - Distance ::

	def distance(self):
		dist_list = []
		ite_list  = []
		for i in range(len(self.x1) - 20):
			d = np.sqrt( ( np.abs(self.x2[i*10] - self.x1[i]) )**2 + ( np.abs(self.y2[i*10] - self.y1[i]) )**2 )
			dist_list.append(d)
			ite_list.append(i)

		plt.plot(ite_list, dist_list)
		plt.title("Difference in trajectories as a function of the iteration", fontname = 'serif', size = 13)
		plt.xlabel("i", fontname = 'serif', size = 13)
		plt.ylabel("Distance", fontname = 'serif', size = 13)
		plt.show()

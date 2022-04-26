import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import statistics
from mpl_toolkits.mplot3d import Axes3D

class plot() :
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Constructor ::

	def __init__(self, N, x, y, z, impactx, impacty, impactz, poincarex, poincarevx, time, energy, angmom):
		self.N          = N
		self.x          = x
		self.y          = y
		self.z          = z
		self.impactx    = impactx
		self.impacty    = impacty
		self.impactz    = impactz
		self.poincarex  = poincarex
		self.poincarevx = poincarevx
		self.time       = time
		self.energy     = energy
		self.angmom     = angmom

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Functions - Plot ::

	def TrajectoryPlot(self): #Plot of the trajectory
		fig = plt.figure()
		ax  = plt.axes(projection='3d')
		for i in range(self.N) :
			ax.plot3D(self.x[i], self.y[i], self.z[i], label = "Trajectory", linewidth = 1)
			ax.scatter3D(self.impactx[i], self.impacty[i], self.impactz[i], c="red", label = "Intersection points")
			# ax.set_title("Trajectory of star", fontname = 'serif', size = 13)
			ax.set_xlabel("x", fontname = 'serif', size = 13)
			ax.set_ylabel("y", fontname = 'serif', size = 13)
			ax.set_zlabel("z", fontname = 'serif', size = 13)
			#ax.legend()
		plt.savefig("TrajectoryPlot.pdf")
		plt.show()

	def PoincarePlot(self): #Plot of the poincare section
		for i in range(self.N):
			plt.scatter(self.poincarex[i], self.poincarevx[i], s = 1)
			plt.title("Poincare map", fontname = 'serif', size = 13)
			plt.xlabel("Position along x", fontname = 'serif', size = 13)
			plt.ylabel("Velocity along x", fontname = 'serif', size = 13)
		plt.savefig("PoincarePlot.pdf")
		plt.show()

	def EnergyPlot(self): #Plot of the energy as a function of time
		for i in range(self.N):
			print(statistics.stdev(self.energy[i]))
			plt.plot(self.time[i], self.energy[i])
			#plt.title("Energy as a function of time", fontname = 'serif', size = 13)
			plt.xlabel("Time", fontname = 'serif', size = 13)
			plt.ylabel("Energy", fontname = 'serif', size = 13)
			#plt.ylim([0., 0.1])
		plt.savefig("EnergyPlot.pdf")
		plt.show()

	def AngMomPlot(self): #Plot of the angular momentum as a function of time
		for i in range(self.N):
			plt.plot(self.time[i], self.angmom[i])
			plt.title("Angular momentum as a function of time", fontname = 'serif', size = 13)
			plt.xlabel("Time", fontname = 'serif', size = 13)
			plt.ylabel("Energy", fontname = 'serif', size = 13)
		plt.show()
		plt.savefig("AngMomPlot.pdf")



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Functions - Animation ::

	def func(self, num, dataSet, line, redDots):
		line.set_data(dataSet[0:2, :num])
		line.set_3d_properties(dataSet[2, :num])
		redDots.set_data(dataSet[0:2, :num])
		redDots.set_3d_properties(dataSet[2, :num])
		return line

	def animation(self):
		X = []                   #Lists that will contain 1 point out of a 100 to simplify the animation
		Y = []
		Z = []
		i = 0

		while i!=len(self.x)/10: #Only plots a tenth of the total trajectory
			X.append(self.x[i])
			Y.append(self.y[i])
			Z.append(self.z[i])
			i+=100

		dataSet = np.array([X, Y, Z])
		numDataPoints = len(X)

		fig = plt.figure()
		ax = Axes3D(fig)
		redDots = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=2, c='r', marker='o')[0] # For scatter plot
		line = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=2, c='g')[0] # For line plot

		ax.set_xlabel("x", fontname = 'serif', size = 13)
		ax.set_ylabel("y", fontname = 'serif', size = 13)
		ax.set_zlabel("z", fontname = 'serif', size = 13)
		ax.set_title("This is the 3D plot of the trajectory in \n 3D space of our 3D world", fontname = 'serif', size = 13)

		line_ani = animation.FuncAnimation(fig, self.func, frames=numDataPoints, fargs=(dataSet,line,redDots), interval=1, blit=False)

		plt.show()

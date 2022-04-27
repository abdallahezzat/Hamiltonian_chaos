import numpy as np
import matplotlib.pyplot as plt
import random 
from star_class import star
from plot_class import plot
from deviation_class import deviation

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Load data from data.txt ::

input     = open('data.txt','r')
data_file = input.readlines()
N  = int(data_file[0].split()[2])     #Number of stars
dt = float(data_file[1].split()[2])   #Time step
T  = float(data_file[2].split()[2])   #Max time
b = []								  #b parameter for each stars
c = []								  #c parameter for each stars

input.close()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Lists ::

star_list      = []
map_list       = []
xyz            = [[], [], []]  #List that will contain the coordinates for the trajectories
impactxyz      = [[], [], []]  #List that will contain the coordinates of the points of the trajectory crossing the plan y=0
poincare       = [[], []]      #List that will contain the coordinates for the poincare map
time           = []
energy         = []
angmom         = []

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Enter the initial conditions manually for all the stars ::

init_condition = [[0.2, 0., 0., 0.3, 0.3, 0.1]]*N
# init_condition = [[]]
# init_condition[0] = [0.2, 0.7, 0, 0.3, 0.3, 0.1]
# 
# init_condition.append([0.2, 0.7, 0, 0.3, 0.3, 0.1])
# init_condition.append([0.2, 0.7, 0, 0.3, 0.3, 0.1])
# init_condition.append([0.2, 0.7, 0, 0.3, 0.3, 0.1])

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Generating random values for the parameters b and c and writing them in a text file ::
prm = open("parameters.txt", "x")
prm.write("Star{}b{}c \n".format("    ","      "))
for i in range(1, N+1):
    b.append(round(random.uniform(0.01,1.0), 2))
    c.append(round(random.uniform(0.01,b[i-1]),2))
    assert c[i-1] <= b[i-1]
    assert b[i-1] <= 1.0
    prm.write("%s     %s \n" % (b[i-2], c[i-2]))
prm.close()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Loop over the N stars ::

for i in range(N) :
	print(init_condition[i])

	f = star(b[i], c[i], dt)                    #Creates a star object
	map_list.append(f)
	map_list[i].rk4(init_condition[i], T) #Choose the integrator

	xyz[0].append(map_list[i].GetPosition()[0])
	xyz[1].append(map_list[i].GetPosition()[1])
	xyz[2].append(map_list[i].GetPosition()[2])

	impactxyz[0].append(map_list[i].GetVisualizer()[0])
	impactxyz[1].append(map_list[i].GetVisualizer()[1])
	impactxyz[2].append(map_list[i].GetVisualizer()[2])

	poincare[0].append(map_list[i].GetXmap())
	poincare[1].append(map_list[i].GetVxmap())

	time.append(map_list[i].GetTime())

	energy.append(map_list[i].GetEnergy())

	angmom.append(map_list[i].GetAngMom())

	print("Star", i+1, b[i], c[i])


p = plot(N, xyz[0], xyz[1], xyz[2], impactxyz[0], impactxyz[1], impactxyz[2], poincare[0], poincare[1], time, energy, angmom)
p.TrajectoryPlot()  #Plot of the trajectory
p.PoincarePlot()	#Plot of the poincare section
p.EnergyPlot()      #Plot of the energy as a function of time
#p.AngMomPlot()     #Plot of the angular momentum as a function of time

#p.animation() #Only use for one star

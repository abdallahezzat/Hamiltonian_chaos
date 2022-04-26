import numpy as np

class star() :
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Constructor ::

	def __init__(self, b, c, dt):
		self.b              = b
		self.c              = c
		self.dt             = dt
		self.tlist          = []
		self.xmap           = []
		self.vxmap          = []
		self.visualizer     = [[],[],[]]
		self.energy_list    = []
		self.angmom_list    = []
		self.position       = [[],[],[]]
		self.velocity       = [[],[],[]]
		self.v_0            = 0
		self.init_parameter = [0,0,0,0,0,0]

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Initial parameter ::

	def set_init_parameter(self, p0):
		self.init_parameter[0] = p0[0]
		self.init_parameter[1] = p0[1]
		self.init_parameter[2] = p0[2]
		self.init_parameter[3] = p0[3]
		self.init_parameter[4] = p0[4]
		self.init_parameter[5] = p0[5]
		self.v_0 = (self.init_parameter[3]**2 + self.init_parameter[4]**2 + self.init_parameter[5]**2)/2

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Getters ::

	def GetInitParameter(self):
		return self.init_parameter

	def GetPosition(self):
		return self.position

	def GetVelocity(self):
		return self.velocity

	def GetXmap(self):
		return self.xmap

	def GetVxmap(self):
		return self.vxmap

	def GetVisualizer(self):
		return self.visualizer

	def GetEnergy(self):
		return self.energy_list

	def GetAngMom(self):
		return self.angmom_list

	def GetTime(self):
		return self.tlist

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Functions - Potentials ::

	def grad_gravityPotential(self, x, y, z):  #To use this potential with Runge-Kutta it is necessary to comment line 81 and uncomment line 80
		d_phi_x = -x*(x**2 + y**2 + z**2)**(-3/2)
		d_phi_y = -y*(x**2 + y**2 + z**2)**(-3/2)
		d_phi_z = -z*(x**2 + y**2 + z**2)**(-3/2)

		return np.array([d_phi_x, d_phi_y, d_phi_z])

	def grad_logPotential(self, x, y, z):  #To use this potential with Runge-Kutta it is necessary to comment line 80 and uncomment line 81
		d_phi_x = -((self.v_0)*x*2)/((y/self.b)**2+(z/self.c)**2+x**2)
		d_phi_y = -((self.v_0)*y*2)/(self.b**2*((y/self.b)**2+(z/self.c)**2+x**2))
		d_phi_z = -((self.v_0)*z*2)/(self.c**2*((y/self.b)**2+(z/self.c)**2+x**2))

		return np.array([d_phi_x, d_phi_y, d_phi_z])

	def F(self, u):
		v, r = u
		#return np.array((self.grad_gravityPotential(r[0], r[1], r[2]), np.array([v[0], v[1], v[2]])))
		return np.array((self.grad_logPotential(r[0], r[1], r[2]), np.array([v[0], v[1], v[2]])))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Functions - Integrators ::

	def rk4(self, p0, T): #Runge-Kutta 4th integrator
		self.set_init_parameter(p0)
		t, i, n = 0, 0, int(T/self.dt)
		X = np.array((n+1)*[np.array([self.init_parameter[:3], self.init_parameter[3:]])]) #initial conditions

		while t < T: #RK algoritm
			k1 = self.dt * self.F(X[i])
			k2 = self.dt * self.F(X[i] + 0.5 * k1)
			k3 = self.dt * self.F(X[i] + 0.5 * k2)
			k4 = self.dt * self.F(X[i] + k3)

			if i < len(X)-1:
				X[i+1] = X[i] + (k1 + 2*(k2 + k3 ) + k4) / 6

				self.position[0].append(X[i+1][0][0])
				self.position[1].append(X[i+1][0][1])
				self.position[2].append(X[i+1][0][2])
				self.velocity[0].append(X[i+1][1][0])
				self.velocity[1].append(X[i+1][1][1])
				self.velocity[2].append(X[i+1][1][2])

				if i > 1:
					self.poincare_map(self.position[1][i-1], self.position[1][i], self.position[2][i], self.velocity[0][i-1],
                                		self.velocity[0][i], self.position[0][i-1], self.position[0][i], self.velocity[1][i-1], self.velocity[1][i])
					self.visualizer_fill(self.position[0][i], self.position[1][i], self.position[2][i], self.position[1][i-1])
					self.energy(self.position[0][i], self.position[1][i], self.position[2][i], self.velocity[0][i], self.velocity[1][i], self.velocity[2][i], t)
					self.angmom(self.position[0][i], self.position[1][i], self.position[2][i], self.velocity[0][i], self.velocity[1][i], self.velocity[2][i])

			t += self.dt
			i += 1

	def euler(self, p0, T): #Euler integrator
		self.set_init_parameter(p0)
		t, i, n = 0, 0, int(T/self.dt)
		X = np.array((n+1)*[np.array([self.init_parameter[:3], self.init_parameter[3:]])]) #Initial conditions

		while t < T: #Euler algorithm
			if i < len(X)-1:
				X[i+1][0][0] = X[i][0][0] + self.dt*X[i][1][0]
				X[i+1][0][1] = X[i][0][1] + self.dt*X[i][1][1]
				X[i+1][0][2] = X[i][0][2] + self.dt*X[i][1][2]

				d_phi_x, d_phi_y, d_phi_z = self.grad_logPotential(X[i+1][0][0], X[i+1][0][1], X[i+1][0][2])

				X[i+1][1][0] = X[i][1][0] + self.dt*d_phi_x
				X[i+1][1][1] = X[i][1][1] + self.dt*d_phi_y
				X[i+1][1][2] = X[i][1][2] + self.dt*d_phi_z

				self.position[0].append(X[i+1][0][0])
				self.position[1].append(X[i+1][0][1])
				self.position[2].append(X[i+1][0][2])
				self.velocity[0].append(X[i+1][1][0])
				self.velocity[1].append(X[i+1][1][1])
				self.velocity[2].append(X[i+1][1][2])

				if i > 1:
					self.poincare_map(self.position[1][i-1], self.position[1][i], self.position[2][i], self.velocity[0][i-1],
                                      self.velocity[0][i], self.position[0][i-1], self.position[0][i], self.velocity[1][i-1], self.velocity[1][i])
					self.visualizer_fill(self.position[0][i], self.position[1][i], self.position[2][i], self.position[1][i-1])
					self.energy(self.position[0][i], self.position[1][i], self.position[2][i], self.velocity[0][i], self.velocity[1][i], self.velocity[2][i], t)
					self.angmom(self.position[0][i], self.position[1][i], self.position[2][i], self.velocity[0][i], self.velocity[1][i], self.velocity[2][i])

			t += self.dt
			i += 1

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#						:: Functions - Other ::

	def poincare_map(self, y1, y0, z0, vx1, vx0, x1, x0, vy1, vy0): #This function allows to store the values for the plot of the Poincare section
		if y1 * y0 < 0:
			self.xmap.append((x1 + x0) / 2)
			self.vxmap.append((vx1 + vx0) / 2)

	def visualizer_fill(self, x, y, z, y1): #This function store the coordinates of the "impact" points with the y=0 plane
		if y1 * y < 0:
			self.visualizer[0].append(x)
			self.visualizer[1].append(y)
			self.visualizer[2].append(z)

	def energy(self, x, y, z, vx, vy, vz, t): #This plot calculates the energy for every coordinates
		EK   = 0.5*(vx**2 + vy**2 + vz**2)
		EPot = self.v_0*np.log(x**2 + (y/self.b)**2 + (z/self.c)**2)
		self.energy_list.append(EK + EPot)
		self.tlist.append(t)

	def angmom(self, x, y, z, vx, vy, vz): #This plot calculates the angular momentum for every coordinates
		Lx = y*vz - z*vy
		Ly = z*vx - x*vz
		Lz = x*vy - y*vx
		self.angmom_list.append(np.sqrt(Lx**2 + Ly**2 + Lz**2))

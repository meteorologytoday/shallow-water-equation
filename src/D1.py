import numpy as np

class D1:

	def __init__(self, dt, N, L, g=10.0, S_fun=None):
		self.dt = dt
		self.N  = N
		self.L  = L
		self.u0 = np.array((N,), dtype=float)
		self.h0 = 
		self.g  = g
		self.S_fun = S_fun


		self.u = np.array(self.u0)
		self.h = np.array(self.h0)



	def run(self):
		















import numpy as np


class RK4:
	
	def __init__(self, N, shape, dt, dysdt_func):
		# N      : number of variables
		# shape  : field shape 
		
		self.total_shape = [N]
		self.total_shape.extend(shape)

		self.ks = [np.zeros(self.total_shape, dtype=float) for _ in range(4)]

		self.dt = dt
		self.dt_2 = dt * 2.0
		self.dysdt_func = dysdt_func

		self.wkspace = np.zeros(self.total_shape, dtype=float)


	def update(self, ys_in, ys_out):
		# [ys_in] and [ys_out] should have dimension (N, shape)
		#
		# 	This function take [ys_int] as initial value ad write
		# new values into [ys_out].
		#
		# 	This function calls [self.dysdt] as provided by user
		# to calcualte k1, k2, k3, and k4.


		# k1
		self.dysdt_func(ys_in, self.ks[0])

		self.ks[0] *= self.dt
		ys_out[:] = self.ks[0]
		self.ks[0] /= 2.0

		# k2
		self.wkspace[:]  = ys_in
		self.wkspace[:] += self.ks[0]

		self.dysdt_func(self.wkspace, self.ks[1])

		self.ks[1] *= self.dt_2
		ys_out += self.ks[1]
		self.ks[1] /= 4.0

		# k3
		self.wkspace[:]  = ys_in
		self.wkspace[:] += self.ks[1]

		self.dysdt_func(self.wkspace, self.ks[2])

		self.ks[2] *= self.dt_2
		ys_out += self.ks[2]
		self.ks[2] /= 2.0
	
		# k4
		self.wkspace[:]  = ys_in
		self.wkspace[:] += self.ks[2]

		self.dysdt_func(self.wkspace, self.ks[3])

		self.ks[3] *= self.dt
		ys_out += self.ks[3]
	
	
		# ys_out
		ys_out /= 6.0
		ys_out += ys_in	


class D1:

	def __init__(self, dt, N, L, g=10.0, S_fun=None):
		self.dt = dt
		self.N  = N
		self.L  = L
		self.y0s = np.array((2, N), dtype=float)
		self.g  = g
		self.S_fun = S_fun

		self.dx   = L / float(N)
		self.dx_2 = self.dx * 2.0

		self.ys = np.zeros((2, N), dtype=float)

		self.dudx = np.zeros((N,), dtype=float)
		self.dhdx = np.zeros_like(self.dudx)

		self.ddudxdx = np.zeros((N,), dtype=float)
		self.ddhdxdx = np.zeros((N,), dtype=float)
		
		self.dysdt = np.zeros_like(self.ys)

		self.RK4 = RK4(2, (N,), dt, self.dysdt_func)



	def ddx(self, Q, dQdx):
		dQdx[:] = (np.roll(Q, -1) - np.roll(Q, 1)) / (self.dx_2)

	def dddxdx(self, Q, ddQdxdx):
		ddQdxdx[:] = (np.roll(Q, -1) + np.roll(Q, 1) - Q * 2.0) / (self.dx ** 2.0)


	def update(self):
		self.u += self.dt * self.dudt
		self.h += self.dt * self.dhdt


	def dysdt_func(self, ys_in, dysdt):
		# 0: u, 1: h
		self.ddx(ys_in[0], self.dudx)
		self.ddx(ys_in[1], self.dhdx)
		
		self.dddxdx(ys_in[0], self.ddudxdx)
		self.dddxdx(ys_in[1], self.ddhdxdx)

		dysdt[0,:] = - ys_in[0] * self.dudx - self.g * self.dhdx   + self.ddudxdx * .005
		dysdt[1,:] = - ys_in[0] * self.dhdx - ys_in[1] * self.dudx + self.ddhdxdx * .005
		



	def run(self, steps, record_interval, u0, h0):

		records = 1 + int(float(steps) / float(record_interval))
		self.u_rec = np.zeros((records, self.N), dtype=float)
		self.h_rec = np.zeros_like(self.u_rec)


		self.u_rec[0] = u0
		self.h_rec[0] = h0


		ys_in  = np.zeros_like(self.ys)
		ys_out = ys_in.copy()

		ys_in[0, :] = u0
		ys_in[1, :] = h0
 
		t = 0

		for rec in range(1, records):
			for step in range(record_interval):
				print("t = %.2f" % t)
				
				self.RK4.update(ys_in, ys_out)
				ys_in[:] = ys_out
				

				t += self.dt

			print("Save array.")
			self.u_rec[rec] = ys_out[0]
			self.h_rec[rec] = ys_out[1]
	













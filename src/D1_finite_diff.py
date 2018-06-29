import numpy as np

class D1:

	def __init__(self, dt, N, L, g=10.0, S_fun=None):
		self.dt = dt
		self.N  = N
		self.L  = L
		self.u0 = np.array((N,), dtype=float)
		self.h0 = np.zeros_like(self.u0) 
		self.g  = g
		self.S_fun = S_fun

		self.dx   = L / float(N)
		self.dx_2 = self.dx * 2.0

		self.u = np.zeros((N,), dtype=float)
		self.h = np.zeros_like(self.u)

		self.dudx = np.zeros_like(self.u)
		self.dhdx = np.zeros_like(self.u)

		self.dudt = np.zeros_like(self.u)
		self.dhdt = np.zeros_like(self.u)

	def ddx(self, Q, dQdx):
		dQdx[:] = (np.roll(Q, -1) - np.roll(Q, 1)) / (self.dx_2)

	def calDudt(self):
		self.dudt[:] = - self.u * self.dudx - self.g * self.dhdx

	def calDhdt(self):
		self.dhdt[:] = - self.u * self.dhdx - self.h * self.dudx

	def update(self):
		self.u += self.dt * self.dudt
		self.h += self.dt * self.dhdt


	def run(self, steps, record_interval, u0, h0):

		records = 1 + int(float(steps) / float(record_interval))
		self.u_rec = np.zeros((records, self.N), dtype=float)
		self.h_rec = np.zeros_like(self.u_rec)


		self.u_rec[0] = u0
		self.h_rec[0] = h0



		self.u[:] = u0
		self.h[:] = h0

		t = 0

		for rec in range(1, records):
			for step in range(record_interval):
				print("t = %.2f" % t)
				# cal dudx, dhdx
				self.ddx(self.u, self.dudx)
				self.ddx(self.h, self.dhdx)
			
				self.calDudt()
				self.calDhdt()

				self.update()
				t += self.dt

			print("Save array.")
			self.u_rec[rec] = self.u
			self.h_rec[rec] = self.h
	













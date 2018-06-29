import numpy as np
from D1_finite_diff import D1

dt = 0.05
N  = 100
L  = 10.0

steps = 10000
record_interval = 10

x = np.arange(N, dtype=float) * L/N

u0 = np.zeros_like(x) * 0.0
h0 = np.exp(-(x - L/2.0)**2.0 / 1.0**2.0 / 2.0) / (2.0 * np.pi)**0.5 / 1.0 / 10

sandbox = D1(
	dt = dt,
	N  = N,
	L  = L
)

sandbox.run(
	steps=steps,
	record_interval = record_interval,
	u0 = u0,
	h0 = h0
)




import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize=(8,6))

ax1.plot(x, sandbox.h_rec[0])
ax1.plot(x, sandbox.h_rec[1])

ax2.plot(x, sandbox.u_rec[0])
ax2.plot(x, sandbox.u_rec[1])


for i in range(0, 5):
	ax1.plot(x, sandbox.h_rec[i])
	ax2.plot(x, sandbox.u_rec[i])


plt.show()

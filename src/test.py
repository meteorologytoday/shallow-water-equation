import numpy as np
from D1_finite_diff import D1

dt = 0.02
N  = 500
L  = 10.0

steps = 10000
record_interval = 10

x = np.arange(N, dtype=float) * L/N

u0 = np.sin((x-2.0) * np.pi * 4.0 / L) * 0.1 + (np.random.rand(N) - 0.5) * 0.2
h0 =  0.01 + np.exp(-(x - L/2.0)**2.0 / 1.0**2.0 / 2.0) / (2.0 * np.pi)**0.5 / 1.0 / 10 + np.sin((x - 2.0) * np.pi * 4.0 / L) * 0.01


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
import matplotlib.animation as animation

fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True, figsize=(8,6))

"""
ax1.plot(x, sandbox.h_rec[0])
ax1.plot(x, sandbox.h_rec[1])

ax2.plot(x, sandbox.u_rec[0])
ax2.plot(x, sandbox.u_rec[1])
"""

cnt = 0

ax1.set_ylim([0,np.amax(h0) * 2.0])
ax2.set_ylim([-1, 1])

ax2.set_xlim([0, L])

line1, = ax1.plot(x, sandbox.h_rec[0])
line2, = ax2.plot(x, sandbox.u_rec[0])

def update(data):
    global ax1, ax2, cnt, line1
    line1.set_ydata(sandbox.h_rec[cnt])
    line2.set_ydata(sandbox.u_rec[cnt])

    cnt = (cnt + 1) % len(sandbox.h_rec)

    return line1,

def data_gen():
    while True:
        yield None

    print("data_gen ends")

ani = animation.FuncAnimation(fig, update, data_gen, interval=50)

plt.show()



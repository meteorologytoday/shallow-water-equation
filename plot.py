import numpy as np
import matplotlib as mplt

mplt.use("Agg")

import matplotlib.pyplot as plt

from matplotlib import rc

#from mycolormap import import cmap_vorticity

nx = 768
ny = 768

Lx = 600000.0
Ly = 600000.0

x_range = np.array([0.0, Lx])
y_range = np.array([0.0, Ly])

x_vec = np.linspace(0, Lx, nx)
y_vec = np.linspace(0, Ly, ny)

t_range = np.array([0, 12.0])

dt = 3.0
record_step = 1
total_steps = int(12 * 60 * 60 / dt)

# 'streamline', 'barb', or 'none'
draw_wind = 'barb'
barb_skip=48



rc('font', **{'family':'sans-serif', 'serif': 'Bitstream Vera Serif', 'sans-serif': 'MS Reference Sans Serif', 'size': 20.0, 'weight' : 100})
rc('axes', **{'labelsize': 20.0, 'labelweight': 100})
rc('mathtext', **{'fontset':'stixsans'})

default_linewidth = 2.0
default_ticksize = 10.0

mplt.rcParams['lines.linewidth'] =   default_linewidth
mplt.rcParams['axes.linewidth'] =    default_linewidth
mplt.rcParams['xtick.major.size'] =  default_ticksize
mplt.rcParams['xtick.major.width'] = default_linewidth
mplt.rcParams['ytick.major.size'] =  default_ticksize
mplt.rcParams['ytick.major.width'] = default_linewidth
mplt.rcParams['xtick.major.pad']='12'
mplt.rcParams['ytick.major.pad']='12'

in_dir  = "output"
out_dir = "output_img"

x_vec = x_vec / 1000.0
y_vec = y_vec / 1000.0


figdpi = 100

t_cnt = 0.0

for step in range(0, total_steps, record_step):
    print("Step: %d" % step)
    #i = int(step / record_step)
    try:
        vort = np.fromfile("%s/vort_step_%04d.bin" % (in_dir, step), dtype='<f4', count=(nx*ny)).reshape((nx, ny)).transpose()
        divg = np.fromfile("%s/divg_step_%04d.bin" % (in_dir, step), dtype='<f4', count=(nx*ny)).reshape((nx, ny)).transpose()
        geop = np.fromfile("%s/geop_step_%04d.bin" % (in_dir, step), dtype='<f4', count=(nx*ny)).reshape((nx, ny)).transpose()
        u = np.fromfile("%s/u_step_%04d.bin" % (in_dir, step), dtype='<f4', count=(nx*ny)).reshape((nx, ny)).transpose()
        v = np.fromfile("%s/v_step_%04d.bin" % (in_dir, step), dtype='<f4', count=(nx*ny)).reshape((nx, ny)).transpose()
    
    except IOError as e:
        print(e)
        print(e.strerror)
        raise Exception("Error when loading files... quit program.")
    

    t_cnt = step * dt

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(12, 6));

    ax1.set_title(r'$\zeta$ [$\times\,10^{-3}\,\mathrm{s}^{-1}$]', fontsize=20)
    ax2.set_title(r'$\delta$ [$\mathrm{s}^{-1}$]', fontsize=20)
    ax3.set_title(r'$\Phi$ [$\mathrm{m}^2 \, \mathrm{s}^{-2}$]', fontsize=20)
    
    ax1.text(10, 10, "%02d:%02d:%02d" % (int(t_cnt/3600), int(t_cnt/60) % 60, t_cnt % 60))

    # VORTICITY
    ax1.set_xlim([x_vec[0], x_vec[-1]])
    ax1.set_ylim([y_vec[0], y_vec[-1]])
    ax1.set_xlabel(r'x [$\mathrm{km}$]')
    ax1.set_ylabel(r'y [$\mathrm{km}$]')
    ax1.set_aspect(1)

    cmap = plt.get_cmap("jet")
    cbar_mappable = ax1.contourf(x_vec, y_vec, vort * 1000.0, cmap=cmap)
    cbar = fig.colorbar(cbar_mappable, ax=ax1, orientation='horizontal')

    ax1.barbs( x_vec[::barb_skip],
              y_vec[::barb_skip],
              u[::barb_skip,::barb_skip] * 0.5144,
              v[::barb_skip,::barb_skip] * 0.5144,
              length=8
    )

    
    # DIVERGENCE
    ax2.set_xlim([x_vec[0], x_vec[-1]])
    ax2.set_ylim([y_vec[0], y_vec[-1]])
    ax2.set_xlabel(r'x [$\mathrm{km}$]')
    ax2.set_ylabel(r'y [$\mathrm{km}$]')
    ax2.set_aspect(1)

    cmap = plt.get_cmap("jet")
    cbar_mappable = ax2.contourf(x_vec, y_vec, divg, cmap=cmap)
    cbar = fig.colorbar(cbar_mappable, ax=ax2, orientation='horizontal')


    # GEOPOTENTIAL HEIGHT
    ax3.set_xlim([x_vec[0], x_vec[-1]])
    ax3.set_ylim([y_vec[0], y_vec[-1]])
    ax3.set_xlabel(r'x [$\mathrm{km}$]')
    ax3.set_ylabel(r'y [$\mathrm{km}$]')
    ax3.set_aspect(1)

    cmap = plt.get_cmap("jet")
    cbar_mappable = ax3.contourf(x_vec, y_vec, geop, cmap=cmap)
    cbar = fig.colorbar(cbar_mappable, ax=ax3, orientation='horizontal')


    file_out = "%s/step_%04d.png" % (out_dir, step) 
    fig.savefig(file_out, dpi=figdpi, format='png')
    print("Output image: %s" % (file_out,))

    plt.close(fig)




#==============================================
# plotmodel.py
#==============================================
import os
import sys
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# Add ISSM functions to current path:
sys.path.append(os.getenv('ISSM_DIR') + '/bin')
sys.path.append(os.getenv('ISSM_DIR') + '/lib')
# sys.path.append(os.getenv('ISSM_DIR') + '/share')
# sys.path.append(os.getenv('ISSM_DIR') + '/share/proj')

#----------------------------------------------------
os.chdir("/Users/dora/GitHub/axel-thk")
#
from loadmodel import loadmodel
from plotmodel import plotmodel



rname = 'axeltest'
# sys.exit(0)

#===============================================================================
#
# REGION / MESHES
#
#===============================================================================

print('Load initial mesh')
mname = "models/" + rname + '_initmesh.nc'
md = loadmodel(mname)

# === Initial mesh =====#
fname = "figs/" + rname + "_initmesh.png"
plotmodel(md,'data','mesh','title','initial mesh')
plt.savefig(fname, dpi=300)



print('Load refined mesh')
mname = "models/" + rname + '_mesh.nc'
md = loadmodel(mname)

#=== Refined mesh =====#
fname = "figs/" + rname + "_mesh.png"
plotmodel(md,'data','mesh','title','mesh')
plt.savefig(fname, dpi=300)

# # sys.exit(0)

#===============================================================================
#
# MODEL PARAMS
#
#===============================================================================

print('Load model params')
mname = "models/" + rname + '_input.nc'
md = loadmodel(mname)



#=== Velocity obs =====#
fname = "figs/" + rname + "_velobs.png"
plotmodel(md, 'data', md.initialization.vel, 'caxis',[0, 200], 'title','velocities [m/yr]')
plt.savefig(fname, dpi=300)

#=== Surface elevation =====#
fname = "figs/" + rname + "_surf.png"
plotmodel(md, 'data', md.geometry.surface, 'caxis',[0, 1800], 'title', 'ice surface elev [m]')
plt.savefig(fname, dpi=300)

#=== Ice thickness =====#
fname = "figs/" + rname + "_thk.png"
plotmodel(md, 'data', md.geometry.thickness, 'caxis',[0, 950], 'title', 'ice thickness [m]')
plt.savefig(fname, dpi=300)

#=== SMB =====#
fname = "figs/" + rname + "_smb.png"
plotmodel(md, 'data', md.smb.mass_balance,'caxis',[-2.6, 0.6],'title', 'smb [m/yr]')
plt.savefig(fname, dpi=300)

#=== dh/dt =====#
fname = "figs/" + rname + "_dhdt.png"
plotmodel(md, 'data', md.balancethickness.thickening_rate, 'caxis',[-5, 5], 'title', 'dhdt [m/yr]')
plt.savefig(fname, dpi=300)

#=== Thickness constraints =====#
fname = "figs/" + rname + "_thk_constraints.png"
plotmodel(md, 'data', md.balancethickness.spcthickness, 'caxis', [0, 650], 'title', 'thickness constraints [m]')
plt.savefig(fname, dpi=300)


# ===============================================================================
#
# PLOT RESULTS # plt.show()
#
# ===============================================================================

print('Load results model')
mname = "models/" + rname + '_donesolve.nc'
md = loadmodel(mname)


#=== Balance thickness results =====#
fname = "figs/" + rname + "_results_thk_full.png"
plotmodel(md, 'data', md.results.BalancethicknessSolution.Thickness, 'title', 'balance thickness solution m')
plt.savefig(fname, dpi=300)

fname = "figs/" + rname + "_results_thk_lim.png"
plotmodel(md, 'data', md.results.BalancethicknessSolution.Thickness, 'title', 'balance thickness solution m', 'caxis', [0,1000])
plt.savefig(fname, dpi=300)

fname = "figs/" + rname + "_boundary_conditions.png"
plotmodel(md, 'data', 'BC')
plt.savefig(fname, dpi=300)

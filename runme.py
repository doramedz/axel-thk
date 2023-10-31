#==============================================
# runme.py
#==============================================
import os
import sys
import numpy as np
import netCDF4 as nc
# from scipy.io import loadmat
# from osgeo import gdal
import matplotlib.pyplot as plt

# Add ISSM functions to current path:
sys.path.append(os.getenv('ISSM_DIR') + '/bin')
sys.path.append(os.getenv('ISSM_DIR') + '/lib')
# sys.path.append(os.getenv('ISSM_DIR') + '/share')
# sys.path.append(os.getenv('ISSM_DIR') + '/share/proj')

#----------------------------------------------------
os.chdir("/Users/dora/GitHub/axel-thk")

rname = 'axeltest'
# 
from model import *
from triangle import triangle
from bamg import bamg
from InterpFromGridToMesh import InterpFromGridToMesh
from export_netCDF_DM import export_netCDF    # edited file in "./bin_edit/"
from loadmodel import loadmodel
from plotmodel import plotmodel
from verbose import verbose
from socket import gethostname
from generic import generic
from solve import solve
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from SetIceSheetBC import SetIceSheetBC

# Not use for now
# from cuffey import cuffey
# from paterson import paterson
# from read_netCDF import netCDFRead
# from m1qn3inversion import m1qn3inversion

# For meshing radar tracks??... not work
# from ContourToNodes import ContourToNodes
# from ContourToMesh import ContourToMesh
# from BamgTriangulate import BamgTriangulate
# from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d


# Steps: 1=mesh, 2=parameterize, 3=solve
steps = [3]



#===============================================================================
#
# REGION / MESH: STEP 1
#
#===============================================================================

if 1 in steps:
  # Step 1: Mesh creation {{{
  print('   Step 1: Make mesh')

  #Generate observations
  domain = 'data/domain_rgiX.exp'
  hinit = 500 		# element size for the initial mesh
  hmax = 5000 		# maximum element size of the final mesh (low res in regions of slow flow)
  hmin = 100 		# minimum element size of the final mesh (high res in regions of fast flow)
  gradation = 1.5 	# maximum size ratio between two neighboring elements
  errmax = 5 			# maximum error between interpolated and control field

  # #===== TO DO: FOLLOW TRACKS ==========================
  # # Load thickness observations (point measurements)
  # tracks = 'data/thk_obs_axeltest.exp' # exp, csv

  # Generate an initial uniform mesh (resolution = hinit m)
  md = bamg(model(),'domain',domain,'hmax',hinit) ## ,'tracks',tracks)
  #Name and Coordinate system
  md.miscellaneous.name = rname
  md.mesh.epsg = 3413


  mname = "models/" + rname + '_initmesh.nc'
  export_netCDF(md, mname)

  # sys.exit(0)

#=============================================x
#===== VELOCITY DATA (m/yr) ===================
  # Load velocities: ITSLIVE
  print('       ----Load velocities')
  ncdata = nc.Dataset('data/velxy_NCAA02_120m.nc', mode='r') # CAN_G0120_0000.nc

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get velocities
  vx = np.squeeze(ncdata.variables['vx'][:].data)
  vy = np.squeeze(ncdata.variables['vy'][:].data)
  ncdata.close()

  # Interpolate velocities onto coarse mesh
  vx_obs = InterpFromGridToMesh(x, np.flipud(y), np.flipud(vx), md.mesh.x, md.mesh.y, 0)
  vy_obs = InterpFromGridToMesh(x, np.flipud(y), np.flipud(vy), md.mesh.x, md.mesh.y, 0)
  vel_obs = np.sqrt(vx_obs ** 2 + vy_obs ** 2)
  del vx, vy, x, y

#=============================================x
#===== WIDE SWATH RADAR DATA (m) ==============
  # Load wide swath thickness: IceBridge
  print('       ----Load radar data')
  ncdata = nc.Dataset('data/rds_wide_thk_2014.nc', mode='r')

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get thickness
  dat = np.squeeze(ncdata.variables['thk'][:].data)
  ncdata.close()

  # Interpolate onto coarse mesh
  thk_obs = InterpFromGridToMesh(x, np.flipud(y), np.flipud(dat), md.mesh.x, md.mesh.y, 0)
  del x, y, dat



  #===== MAKE BAMG MESH =====
  # Build combined velocity + thickness fields
  field = np.array([vel_obs, thk_obs]).T

  # Combine errors
  err = errmax * np.ones(md.mesh.numberofvertices)
  err = np.array([err, err]).T


  # Adapt the mesh to minimize error in velocity AND thickness interpolation
  print('       ----Adapt mesh')
  md = bamg(md,'hmax',hmax,'hmin',hmin,'gradation',gradation,'field',field,'err',err)


  print('       ----Export mesh model')
  mname = "models/" + rname + '_mesh.nc'
  export_netCDF(md, mname)



#===============================================================================
#
# PARAMETERISE: STEP 2
#
#===============================================================================

if 2 in steps:
  # Step 2: Parameterise {{{
  print('   Step 2: Load data, Set params')
  mname = "models/" + rname + '_mesh.nc'

  md = loadmodel(mname)

#=============================================x
#===== VELOCITY (m/yr) ========================
  # Load velocities: ITSLIVE
  print('       ----Load velocities')
  ncdata = nc.Dataset('data/velxy_NCAA02_120m.nc', mode='r') # CAN_G0120_0000.nc

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get velocities
  vx = np.squeeze(ncdata.variables['vx'][:].data)
  vy = np.squeeze(ncdata.variables['vy'][:].data)
  ncdata.close()

  # Interpolate onto mesh
  vx_obs = InterpFromGridToMesh(x, np.flipud(y), np.flipud(vx), md.mesh.x, md.mesh.y, 0)
  vy_obs = InterpFromGridToMesh(x, np.flipud(y), np.flipud(vy), md.mesh.x, md.mesh.y, 0)
  vel_obs = np.sqrt(vx_obs ** 2 + vy_obs ** 2)
  del x, y, vx, vy

#=============================================x
#===== SURFACE ELEVATION (m) ==================
  # Load dem
  print('       ----Load surf elev')
  ncdata = nc.Dataset('data/arcticdem_NCAA02_2f100m.nc', mode='r')

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get elevation
  dat = np.squeeze(ncdata.variables['elev'][:].data)
  ncdata.close()

  # Interpolate onto mesh
  surf = InterpFromGridToMesh(x, np.flipud(y), np.flipud(dat), md.mesh.x, md.mesh.y, 0)
  del x, y, dat

#=============================================x
#===== ICE THICKNESS (m) ======================
  # Load ice thickness
  print('       ----Load ice thickness')
  ncdata = nc.Dataset('data/thickness_NCAA02_100m.nc', mode='r')

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get thickness
  dat = np.squeeze(ncdata.variables['thick'][:].data)
  ncdata.close()

  # Interpolate onto mesh
  thk = InterpFromGridToMesh(x, np.flipud(y), np.flipud(dat), md.mesh.x, md.mesh.y, 0)
  del x, y, dat

#=============================================x
#===== SMB (mm/yr w.e.) =======================
  # Load SMB
  print('       ----Load SMB')
  # ncdata = nc.Dataset('SMB_rec_1996-2015.CAA_North_1km.YYmean.nc', mode='r')
  ncdata = nc.Dataset('data/SMB_NCAA_2f1km.nc', mode='r')

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get SMB
  dat = np.squeeze(ncdata.variables['SMB_rec'][:].data)
  ncdata.close()

  # Interpolate onto mesh
  smb = InterpFromGridToMesh(x, np.flipud(y), np.flipud(dat), md.mesh.x, md.mesh.y, 0)
  del x, y, dat

#=============================================x
#===== DH / DT (m/yr) =========================
  # Load dh/dt
  print('       ----Load DH/DT')
  ncdata = nc.Dataset('data/dhdt_NCAA02_2010-2020_100m.nc', mode='r')

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get dh/dt
  dat = np.squeeze(ncdata.variables['dhdt'][:].data)
  ncdata.close()

  # Interpolate onto mesh
  dhdt = InterpFromGridToMesh(x, np.flipud(y), np.flipud(dat), md.mesh.x, md.mesh.y, 0)
  del x, y, dat

#=============================================x
#===== WIDE SWATH RADAR DATA (m) ==============
  # Load wide swath thickness: IceBridge
  print('       ----Load radar data')
  ncdata = nc.Dataset('data/rds_wide_thk_2014.nc', mode='r')

  # Get coords
  x = np.squeeze(ncdata.variables['easting'][:].data)
  y = np.squeeze(ncdata.variables['northing'][:].data)
  # Get thickness
  dat = np.squeeze(ncdata.variables['thk'][:].data)
  ncdata.close()

  # Interpolate onto coarse mesh
  thk_obs = InterpFromGridToMesh(x, np.flipud(y), np.flipud(dat), md.mesh.x, md.mesh.y, np.nan)
  del x, y, dat


#=============================================x
#===== SET PARAMS =============================
  # Set params from file
  # md = parameterize(md, 'X.py')

  print('       ----Parameterize')
#===== Initialise velocities
  md.initialization.vx = vx_obs
  md.initialization.vy = vy_obs
  md.initialization.vz = np.zeros((md.mesh.numberofvertices))
  # Remove zero vel
  md.initialization.vx[np.where(md.initialization.vx == 0)] = 1
  md.initialization.vy[np.where(md.initialization.vy == 0)] = 1
  # Norm
  md.initialization.vel = np.sqrt(md.initialization.vx ** 2 + md.initialization.vy ** 2)


  # md.inversion.vx_obs = md.initialization.vx
  # md.inversion.vy_obs = md.initialization.vy
  # md.inversion.vel_obs = md.initialization.vel

#===== Set geometry
  md.geometry.surface = surf
  md.geometry.thickness = thk
  # Set min thickness
  md.geometry.thickness[np.where(md.geometry.thickness <= 0)] = 1

  md.geometry.base = md.geometry.surface - md.geometry.thickness
  md.geometry.bed = md.geometry.surface - md.geometry.thickness

#===== Set mass balance
  # Convert to m/yr ice: materials.rho_water and rho_ice (pre-defined)
  md.smb.mass_balance = smb/1000 # m/yr w.e.
  md.smb.mass_balance = md.smb.mass_balance * md.materials.rho_water / md.materials.rho_ice

#===== Set other boundary conditions
  # Ocean mask (grounded > 0, floating < 0, coast/grounding line == 0)
  md.mask.ocean_levelset = np.ones((md.mesh.numberofvertices)) # grounded == 1)
  # Ice mask (ice < 0, rock > 0, ice front == 0)
  md.mask.ice_levelset = -1 * np.ones((md.mesh.numberofvertices)) # ice == -1)
  # Sets ice front at boundary (ice front == 0)
  md.mask.ice_levelset[np.nonzero(md.mesh.vertexonboundary == 1)] = 0

  md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))  # (positive if melting)
  md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))  # (positive if melting)

  # md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
  # md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
  # md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

  #=== Params for balancethickness solution
  # stabilization # 0: None, 1: SU, 2: SSA's artificial diffusivity, 3: discontinuous Galerkin (DG?)
  md.balancethickness.stabilization = 1
  md.balancethickness.thickening_rate = dhdt

  # thickness constraints: set minimum thickness at domain boundary
  md.balancethickness.spcthickness = np.nan * np.ones((md.mesh.numberofvertices))
  md.balancethickness.spcthickness[np.nonzero(md.mesh.vertexonboundary == 1)] = 1
  # fill in radar thickness measurements
  pos = np.where(~np.isnan(thk_obs))[0]
  md.balancethickness.spcthickness[pos] = thk_obs[pos]


  print('       ----Export model')
  mname = "models/" + rname + '_input.nc'
  export_netCDF(md, mname)

# ===============================================================================
#
# DO SOMETHING: STEP 3
#
#===============================================================================

if 3 in steps:
  # Step 3: Do something {{{
  print('   Step 3: make solve')
  mname = "models/" + rname + '_input.nc'

  md = loadmodel(mname)

# Set cluster #md.cluster (set only the name and number of process)
  md.cluster = generic('name', gethostname(), 'np', 2)
  md.verbose = verbose(5)
# Solve #help solve
  md = solve(md, 'Balancethickness') # Stressbalance

  print('       ----Export results model')
  mname = "models/" + rname + '_donesolve.nc'
  export_netCDF(md, mname)

  print('   Ok done!')

#===============================================================================
# ...
#===============================================================================



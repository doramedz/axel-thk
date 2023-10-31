# axel-thk
 issm test: axel ice thickness

Testing ISSM for modelling ice thickness over single ice mass on southern Axel Heiberg Island, Canadian Arctic.


 ## System and installation
 ISSM 4.22 (Release 2022-10-27) installed from source on Intel-based MacBookPro with OS version 12.6 (Monterey) and Python 3.7.16+ (dev).

 Standard installation with recommended packages: autotools, cmake, petsc, triangle, chaco, m1qn3


 #### Configuration (Python interface)
 configure.sh
 ```
 {
 	./configure \
	--prefix="${ISSM_DIR}" \
	--with-python-dir="${HOME}/.pyenv/versions/3.7-dev" \
	--with-python-numpy-dir="${HOME}/.local/lib/python3.7/site-packages/numpy" \
	--with-fortran-lib="-L/usr/local/Cellar/gcc/13.1.0/lib/gcc/13 -lgfortran" \
	--with-mpi-include="${ISSM_DIR}/externalpackages/petsc/install/include" \
	--with-mpi-libflags="-L${ISSM_DIR}/externalpackages/petsc/install/lib -lmpi -lmpicxx -lmpifort" \
	--with-metis-dir="${ISSM_DIR}/externalpackages/petsc/install" \
	--with-parmetis-dir="${ISSM_DIR}/externalpackages/petsc/install" \
	--with-blas-lapack-dir="${ISSM_DIR}/externalpackages/petsc/install" \
	--with-scalapack-dir="${ISSM_DIR}/externalpackages/petsc/install" \
	--with-mumps-dir="${ISSM_DIR}/externalpackages/petsc/install" \
	--with-petsc-dir="${ISSM_DIR}/externalpackages/petsc/install" \
	--with-triangle-dir="${ISSM_DIR}/externalpackages/triangle/install" \
	--with-chaco-dir="${ISSM_DIR}/externalpackages/chaco/install" \
	--with-m1qn3-dir="${ISSM_DIR}/externalpackages/m1qn3/install"
 }
 ```

## Input datasets

#### Gridded data
Ice surface velocities: [ITSLIVE](https://its-live.jpl.nasa.gov/#about) Regional mosaic. 1985–2018 average m/yr, 120 m

Ice surface elevation: [ArcticDEM](https://www.pgc.umn.edu/data/arcticdem/) Mosaic version 4.1 (2023 release), 100 m

Ice thickness: [Millan et al. (2022)](https://www.sedoo.fr/theia-publication-products/?uuid=55acbdd5-3982-4eac-89b2-46703557938c) Global ice thickness dataset version 1.1. 2010-2020, 50 m (downscaled to 100 m)

Surface mass balance: [Noël (2017)](https://doi.pangaea.de/10.1594/PANGAEA.881315?format=html#download) Surface mass balance components for the Northern Canadian Arctic Archipelago. 1996-2015 average mm w.e. per year, downscaled from RACMO2.3 to 1 km resolution

Ice surface elevation change: [Hugonnet et al. (2021)](https://www.sedoo.fr/theia-publication-products/?uuid=c428c5b9-df8f-4f86-9b75-e04c778e29b9) Global glacier elevation change dataset. 2010-2020 average m/yr, 100 m


#### Glacier outlines
Glacier outlines: [Randolph Glacier Inventory 7.0](https://www.glims.org/RGI/) (RGIv7) Version 9 (2023 Release)
ARGUS files:
domain_rgiX.exp contains domain outline and nunataks (holes) as closed contours
domain_rgiX_outline.exp and domain_rgiX_holes.ext split outline and hole contours in separate files


#### Radar data
Ice thickness observations: Operation Icebridge MCoRDS ice thickness measurements, ~15 m point spacing
ARGUS, CSV files:
thk_obs_axeltest.exp contains XY coordinates of point measurements along flight tracks (split into segments)
thk_obs_axeltest.csv also includes measured ice thickness (as single segment)

Wide-swath OIB MCoRDS ice thickness data, ~15 m along-track and ~15-130 m across-track spacing point spacing
Gridded at 25 m across a 1.5-3 km wide swath along flight track


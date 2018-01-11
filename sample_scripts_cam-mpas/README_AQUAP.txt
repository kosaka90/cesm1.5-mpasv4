README file for aquaplanet experiments using the CESM-MPAS code (cesm2.0 beta05)
Contributed from Jihyeon Jang, Koichi Sakaguchi, and Bryce Harrop

Aquaplanet configuration in CESM2 (Thanks to Brian Medeiros, Jim Benedict, and Colin!)  

Available from cesm2_0_beta06 (But.. CAM-MPAS version: cesm2_0_beta05) 

compset: FC5AQUAP 
- in env_run.xml, CAM_NML_USE_CASE set to 'aquaplanet_cam5' 

- no topography (topo=0), water surface (land fraction=0) in the MPAS initial file as well as the topo file (details will be added)

- set both atmo and ocean resolution to be mpas grid ('-res mp120a-mp120a' for create_newcase)

- Perpetual symmetric forcing with the prescribed solar irradiance   -  eccentricity and obliquity are set to zero 
(user_nl_cpl) 
orb_eccen = 0. 
orb_obliq = 0. 
orb_mvelp = 0. 
orb_mode = "fixed_parameters"


- No aerosol-radiation interaction, No aerosol-cloud interaction  -  Aerosol emissions files are not used. 
(user_nl_atm) 
ext_frc_specifier = ''
srf_emis_specifier = ''
tracer_cnst_specifier = ''
prescribed_aero_datapath='/global/project/projectdirs/m1867/beharrop/MPAS/inputdata/APE_aerosols/'
prescribed_aero_file='mam3_1.9x2.5_L30_1850clim_aquap_equinoctial_c171108.nc'
prescribed_aero_cycle_yr = 1850
prescribed_aero_type = 'CYCLICAL'

-  A constant drop number in the microphysics   
(user_nl_atm) 
micro_mg_nccons = .true.,
micro_mg_nicons = .true.

-Prescribed Aerosol File: 'mam3_1.9x2.5_L30_1850clim_aquap_equinoctial_c171108.nc'
 
The file was generated from aerosol file mam3_1.9x2.5_L30_1850clim_c130319.nc, which is seasonally varying and has realistic land-sea geometry.  First, all aerosol variables are averaged zonally and temporally.  Next, the aerosol variables are made to be symmetric about the equator by averaging the northern and southern hemispheres together and applying the result to both hemispheres.  Finally, a smoother is applied in the meridional direction to smooth any remaining signatures of the original land geometry.  The smoother is a Savitsky Golay filter from the scipy cookbook: http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay.  Each variable is broadcast back to its original coordinate dimensions to match the original file so that no new code is needed to read the file into MPAS.  The surface pressure from the original file could not be used because of the influence of mountain ranges.  Therefore, surface pressure in this aerosol file is taken from a seasonally varying 240 km MPAS aquaplanet simulation.  The same averaging, making symmetric about the equator, and smoothing algorithm is applied to the surface pressure.
 
Since the surface pressure comes from a simulation and not the original data, there may be a mismatch in the aerosol abundances and the mass of the atmosphere.  It is unclear what the impact of this discrepancy would be.  It is also unclear whether all aerosol species from the original aerosol file should exist in an aquaplanet simulation.  This aerosol file was generated under the assumption that the aquaplanet simulation would be run with cloud-aerosol interactions turned off.  If the user desires to examine aerosol interactions in greater detail, careful thought is suggested on whether this aerosol file is appropriate.



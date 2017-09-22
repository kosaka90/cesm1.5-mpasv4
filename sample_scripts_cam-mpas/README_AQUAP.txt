Aquaplanet configuration in CESM2 (Thanks to Brian Medeiros, Jim Benedict, and Colin!)  

Available from cesm2_0_beta06 (But.. CAM-MPAS version: cesm2_0_beta05) 

compset: FC5AQUAP 
- in env_run.xml, CAM_NML_USE_CASE set to ‘aquaplanet_cam5’  

- no topography (topoà0), water surface (land fractionà0) in the MPAS initial file as well as the topo file (details will be added)


- Perpetual symmetric forcing with the prescribed solar irradiance   -  eccentricity and obliquity are set to zero 
(user_nl_cpl) 
orb_eccen = 0. 
orb_obliq = 0. 
orb_mvelp = 0. 
orb_mode = “fixed_parameters”


- No aerosol-radiation interaction, No aerosol-cloud interaction  -  Aerosol emissions files are not used. 
(user_nl_atm) 
ext_frc_specifier = ''
srf_emis_specifier = ''
tracer_cnst_specifier = ''


-  A constant drop number in the microphysics   
(user_nl_atm) 
micro_mg_nccons = .true.,
micro_mg_nicons = .true.


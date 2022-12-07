import lava2d as sim
#
visc_melt_vent_dense = 200
vesic = 0.25
phi_max = 0.6
phi_inf = 0.58
max_cryst_rate = 0.1e-4
yield_strength_crust = 4e4
sim_num = 1
#
#-------------------------------------------------------------------------------
sim.set_topo( # set DEM info
    path_to_dem_file    = ('ned_10m_dem.tif'),
    Lon_SRC             = -155.531147, # source longitude
    Lat_SRC             = 19.513278,   # source latitude
    Lon_LowerLeft       = -155.541, # bounding box: lower-left longitude
    Lat_LowerLeft       = 19.511,   # bounding box: lower-left latitude
    Lon_UpperRight      = -155.44, # bounding box: upper-right longitude
    Lat_UpperRight      = 19.71,   # bounding box: upper-right latitude
    fill_level          = 0.0,    # set sea level to fill up to
    dx_desired          = 40,     # meters
    smooth_n_times      = 1
    )

#
#-------------------------------------------------------------------------------
sim.set_output( # where to store out.nc?
    path_out =  ('sim_{:04d}'.format(sim_num))
    )

#
#-------------------------------------------------------------------------------
sim.set_source( # set vent/fissure info: where is vent_nn.txt located?
    path_to_vent_files      = ('..//example_vents//vent_01.txt')
    )

#-------------------------------------------------------------------------------
sim.set_init( # set initialization type and file
    init_type = None,  # 'prior_model' or None currently supported
    init_file = None
    )

#
#-------------------------------------------------------------------------------
sim.set_vent_props( # set lava properties @ vents
    temperature_vent_C  = 1150, # deg C
    viscosity_melt_vent = visc_melt_vent_dense*(1 - vesic/.64)**-2.5, # Pa s
    cryst_vent          = 0.0  # erupted crystal fraction
    )

#
#-------------------------------------------------------------------------------
sim.set_lava_props( # set lava properties throughout
    liquid_density      = 2700,   # kg m-3
    porosity            = vesic,
    lava_specific_heat  = 1500,   # J kg-1 K-1
    lava_diffusivity    = 5e-7, # m2 s-1
    lava_conductivity   = None,    # W m-1 K-1
    lava_emissivity     = 0.95
    )

#
#-------------------------------------------------------------------------------
# using Avrami n = 4
sim.set_rheo( # set rheological properties
    phi_max                 = phi_max, # max crystal packing fraction (could be 0.6 e.g., Marsh 1981; Pinkerton and Stevenson 1992), # could be higher (Cashman et al., 1999)
    phi_inf                 = phi_inf,
    max_cryst_rate          = max_cryst_rate, # s-1 # max crystalization rate: max d(phi)/dt
    yield_strength_crust    = yield_strength_crust, # Pa
    glass_temperature       = None, # K
    T_core_T_vent_equal     = True  # core temperature equals vent temperature
    )

#
#-------------------------------------------------------------------------------
sim.set_ambient( # set ambient properties
    atm_temperature     = 300, # K
    h_conv              = 50,   # W m-2 K-1
    rainfall            = 0,   # m s-1
    ground_temperature  = 300  # K
    )

#
#-------------------------------------------------------------------------------
sim.set_numerics( # set numerical method details
    method = 'OHA',
    efficiency_min      = 0, # minmum allowable model-clock ratio
    efficiency_max      = 10000,  # maximum allowable model-clock ratio
    cfl_max             = 0.5,
    dt_max              = 5.,
    fraction_to_freeze  = 1.0,  # fraction of freezable lava
    tiny_flow           = 0.1,  # min thickness of lava (m)
    )

#
#-------------------------------------------------------------------------------
sim.set_runtime(
    max_iter = None, #one of them can be none, so the default value for both is none, if both none run till killed
    max_time_hr = 1*24, # hr
    out_times = [i*24. for i in range(1,1)], # hr, list of intermediate output times
    run_to_ocean = True
    )
#
#-------------------------------------------------------------------------------
#start simulation
sim.run()
#
#-------------------------------------------------------------------------------
#

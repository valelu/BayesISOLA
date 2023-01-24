#! /usr/bin/env python

import BayesISOLA

#event="IT-2013-0005"
event="EMSC-20160824_0000006"
inputs=BayesISOLA.load_italy(outdir=event,
        logfile='$outdir/log.txt',
        output_mkdir=True) #here all are optional
inputs.load_from_esm(eventID=event,max_distance=150.e3,invert_Z_only=False)#,proc_type="CV,MP",data_type="ACC",save_to="query.h5")
#inputs.load_asdf(file='query.h5',max_distance=150.e3,min_distance=None,invert_Z_only=False)
#inputs.load_itaca(dir='query/',max_distance=150.e3)
inputs.read_crust('crustal.dat', output='green/crustal.dat') #MUST COME after station read

grid=BayesISOLA.grid( inputs,
        location_unc = 0, # m
	depth_unc = 2000, # m
	time_unc = 0, # s
	step_x = 500, # m
	step_z = 500, # m
	max_points = 225,
        grid_radius= 0., #m
        grid_min_depth = 0.,#m
        grid_max_depth = 0.,#m
        grid_min_time = 0.,#s
        grid_max_time = 0.,#s
	circle_shape = True,
	rupture_velocity = 1000 # m/s
	)


data = BayesISOLA.process_italy( inputs,
	grid,
        s_velocity = 3000, #m/s
	threads = 8,
        invert_displacement = True, 
	use_precalculated_Green = 'auto', #checks if GFs exist and have the same grid
        correct_data = False,
        set_parameters = True, 	
        fmax = 0.08,
	fmin = 0.02, 
        min_depth =0,#m
        multichannel=True,
        skip_short_records=True,
        suspect_mouse=True,
        use_stainfo=True,
#        noise=False, #TODO: check for necessity of this parameter
        calculate_or_verify_Green = True,
        trim_filter_data = True,
        decimate_shift = True)

cova=BayesISOLA.covariance_matrix(data) #Must be initialized

solution = BayesISOLA.resolve_MT(data, 
           cova = cova, 
           deviatoric = True,
           decompose = True, #from here default values, not necessary to be here
           run_inversion = True,
           find_best_grid_point = True,
           save_seismo = False,
           VR_of_components = True,
           print_solution = True,
           print_fault_planes=True)


plot = BayesISOLA.plot(solution,
        maps = True,
        slices = True,
        maps_sum = True,
        MT = True,
 #       uncertainty = 3000, #number of samples for sampling PPDF
        seismo = True, #default is False
        seismo_sharey = False,#True,
        seismo_cova = True, #default is True
        noise = False, #default is True
        spectra =  False, #default is True, we do not have noise spectra
        stations = True,
        covariance_matrix = False, #default is True 
        covariance_function = False,
        pygmt=True)


plot.html_log(h1= inputs.event['name']+":"+event,
        backlink=False,
        plot_MT='auto',
        plot_uncertainty= 'auto',
        plot_stations='auto',
        plot_seismo_cova='auto',
        plot_seismo='auto',
        plot_seismo_sharey = 'auto',
        mouse_figures="mouse/", #we are not using mousetrap
        plot_spectra='auto',
        plot_noise='auto',
        plot_covariance_function='auto',
        plot_covariance_matrix='auto',
        plot_maps='auto',
        plot_slices='auto',
        plot_maps_sum='auto')

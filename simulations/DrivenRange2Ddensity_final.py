{

    'run_time': 12000, # ms
    'dt': 0.1, # ms

    'Populations' : {
        'drive' : {
            # 'n' : 25*25,
            # 'n' : 50*50,
            'n' : 100*100,
            'type': sim.SpikeSourcePoisson,
            'cellparams' : {
                'start':0.0,
                'rate':4.,
                'duration': 120000.0
            }
        },

       'py' : {
            'n': 100*100, # units
            'type': sim.EIF_cond_alpha_isfa_ista,
            'structure' : Grid2D(aspect_ratio=1, dx=1.0, dy=1.0, fill_order='sequential', rng=sim.NumpyRNG(seed=2**32-1)),
            'cellparams': {
                'e_rev_I'    : -80,   # mV, reversal potential of excitatory synapses
                'e_rev_E'    : 0,     # mV, reversal potential of inhibitory synapses
                'tau_syn_E'  : 5.0,   # ms, time constant of excitatory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_syn_I'  : 5.0,   # ms, time constant of excitatory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_refrac' : 5.0,   # ms, refractory period
                'v_reset'    : -65.0, # mV, reset after spike
                'v_thresh'   : -50.0, # mV, spike threshold (modified by adaptation)
                'delta_T'    : 2.,    # mV, steepness of exponential approach to threshold
                'cm'         : 0.150, # nF, tot membrane capacitance
                'a'          : 4.,    # nS, conductance of adaptation variable
                'tau_m'      : 15.0,  # ms, time constant of leak conductance (cm/gl)
                'v_rest'     : -65.0, # mV, resting potential E_leak
                'tau_w'      : 500.0, # ms, time constant of adaptation variable
                'b'          : .02,   # nA, increment to adaptation variable
            },
        },
        'inh' : {
            'n': 50*50, #{'ref':'py','ratio':0.25},
            'type': sim.EIF_cond_alpha_isfa_ista,
            'structure' : Grid2D(aspect_ratio=1, dx=2.0, dy=2.0, fill_order='sequential', rng=sim.NumpyRNG(seed=2**32-1)),
            'cellparams': {
                'e_rev_I'    : -80,   # mV, reversal potential of excitatory synapses
                'e_rev_E'    : 0,     # mV, reversal potential of inhibitory synapses
                'tau_syn_E'  : 5.0,   # ms, time constant of excitatory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_syn_I'  : 5.0,   # ms, time constant of inhibitory synaptic short-term plasticity, YgerBoustaniDestexheFregnac2011
                'tau_refrac' : 5.0,   # ms, refractory period
                'v_reset'    : -65.0, # mV, reset after spike
                'v_thresh'   : -50.0, # mV, spike threshold (modified by adaptation)
                'delta_T'    : 0.5,   # mV, steepness of exponential approach to threshold
                'cm'         : 0.150, # nF, tot membrane capacitance
                'a'          : 0.0,   # nS, conductance of adaptation variable
                'tau_m'      : 15.0,  # ms, time constant of leak conductance (cm/gl)
                'v_rest'     : -65.0, # mV, resting potential E_leak
                'tau_w'      : 500.0, # ms, time constant of adaptation variable
                'b'          : 0.0,   # nA, increment to adaptation variable
            },
        },
    },

    'Projections' : {
        'drive_py' : {
            'pre' : 'drive',
            'post' : 'py',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.FixedProbabilityConnector(.01, rng=sim.random.NumpyRNG(2**32-1)),
            'synapse_type' : sim.StaticSynapse(),
            # 'weight' : .001, # uS # 25*25 *1000 *.008 = 5000
            # 'weight' : .002, # uS # 50*50 *1000 *.002 = 5000
            'weight' : .0005, # uS # 100*100 *1000 *0.005 = 5000
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },
        'drive_inh' : {
            'pre' : 'drive',
            'post' : 'inh',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.FixedProbabilityConnector(.01, rng=sim.random.NumpyRNG(2**32-1)),
            'synapse_type' : sim.StaticSynapse(),
            # 'weight' : .001, # uS # 25*25 *1000 *.008 = 5000
            # 'weight' : .002, # uS # 50*50 *1000 *.002 = 5000
            'weight' : .0005, # uS # 100*100 *1000 *0.005 = 5000
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },

        'py_py' : {
            'pre' : 'py',
            'post' : 'py',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            # 'connector' : sim.DistanceDependentProbabilityConnector("14*exp(-(0.5*numpy.random.randn()+1.2)*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # 200um
            'connector' : sim.DistanceDependentProbabilityConnector("14*exp(-1.2*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # 200um
            'weight' : .002, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms
            'receptor_type' : 'excitatory',
            'save_connections':True,
            'print_statistics':True,
        },
        'py_inh' : {
            'pre' : 'py',
            'post' : 'inh',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("14*exp(-1.2*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # 200um
            'weight' : {'ref':'py_py'}, # µS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms,
            'receptor_type' : 'excitatory',
            'save_connections':False,
            'print_statistics':False,
        },
        'inh_inh' : {
            'pre' : 'inh',
            'post' : 'inh',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("24*exp(-1.5*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # 150um
            'weight' : .005, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms,
            'receptor_type' : 'inhibitory',
            'save_connections':False,
            'print_statistics':False,
        },
        'inh_py' : {
            'pre' : 'inh',
            'post' : 'py',
            'space' :  sim.Space(periodic_boundaries=((0,100), (0,100), None)), # torus
            'connector' : sim.DistanceDependentProbabilityConnector("24*exp(-1.5*d)", allow_self_connections=False, rng=sim.NumpyRNG(2**32-1)), # 150um
            'weight' : {'ref':'inh_inh'}, # µS
            # 'weight' : .005, # uS
            'synapse_type' : sim.StaticSynapse(),
            'delay' : .5, # ms,
            'receptor_type' : 'inhibitory',
            'save_connections':False,
            'print_statistics':False,
        },

    },


    'Recorders' : {
        'py' : {
            'spikes' : 'all',
            # 'v' : {
            #     'start' : 0,
            #     'end' : 100,
            # },
            'gsyn_exc' : {
                'start' : 0,
                'end' : 10,
            },
            'gsyn_inh' : {
                'start' : 0,
                'end' : 10,
            },
        },
        'inh' : {
            'spikes' : 'all',
            # 'v' : {
            #     'start' : 0,
            #     'end' : 100,
            # },
            'gsyn_exc' : {
                'start' : 0,
                'end' : 10,
            },
            'gsyn_inh' : {
                'start' : 0,
                'end' : 10,
            },
        },
    },


    'Modifiers' :{
    },

    'Injections' : {
    },

    'Analysis' : {
        'scores' : ['py'],
        # 'scores' : ['py','inh'],
        #'subsampling' : 1000, # units
        'transient' : 1000, # ms

        'Structure': {
            'py':{
                'conns': 'py_py',
                'shape': (100**2, 100**2),
            },
        },

        'Events_Clustering': {
            'py':{
                # 'add': ['inh'],
                'limits': [(0,100),(0,100)], # only central ones
                'bin': 5, # ms
                # 'bin': 30, # ms (at 2photon resolution)
                'trials': ['default'], # for which trials the analysis has to be computed
                'ylim': [0,20],
                'print2Dmap': True,
            },
        },

        'ConductanceBalance' : {
            'py':{
                'trials': ['default'], # for which trials the analysis has to be computed
            },
            'inh':{
                'trials': ['default'], # for which trials the analysis has to be computed
            },
        },

        'FiringRate' : {
            'bin': 5, # ms
            'py':{
                'firing': [0,20],
            },
            'inh':{
                'firing': [0,20],
            },
        },

        'Rasterplot' : {
            'py':{
                'limits': [(0,100),(0,100)], # coords: [(from x, to x), (from y, to y)]
                'color': 'red',
            },
            'inh':{
                'limits': [(0,100),(0,100)], # coords: [(from x, to x), (from y, to y)]
                'color': 'blue',
            },
            'type': '.png',
            # 'type': '.svg',
            'interval': False, # all
            # 'interval': [2000.,3000.], # ms # from 2s to 3s
            'dpi':800,
        },

        # 'ISI' : {
        #     'py':{
        #         'bin': 50, # ms, 20 per second
        #         # 'limits': [(0,100),(0,100)], # coords: [(from x, to x), (from y, to y)]
        #         'limits': [(10,50),(10,50)], # only central ones
        #     },
        # },
    },

}

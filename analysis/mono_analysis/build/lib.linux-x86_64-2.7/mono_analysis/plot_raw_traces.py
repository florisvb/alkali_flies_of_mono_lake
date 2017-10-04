from optparse import OptionParser
import sys, os
import h5py
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.signal

import fly_plot_lib.text as flytext
import fly_plot_lib.plot as fpl
import flystat
from multi_tracker_analysis import bag2hdf5

import figurefirst

import mono_paper_locations

def make_figure(path, solvent='', oddonly=False, evenonly=False, calc='work', flynum=None, normalize_times=True, axes=None, solution_to_plot='none', make_ax_spines_clean=True, fill_work_out=False, ylim=None):
    
    if axes is None:
        fig = plt.figure(figsize=(8,3))
        ax_before = fig.add_axes([0.15,0.15,0.05,0.8])
        ax_in = fig.add_axes([0.21,0.15,0.3,0.8])
        ax_underwater = fig.add_axes([0.52,0.15,0.05,0.8])
        ax_out = fig.add_axes([0.58,0.15,0.3,0.8])
        ax_after = fig.add_axes([0.89,0.15,0.05,0.8])
        axes = {'before': ax_before,
                'in': ax_in,
                'underwater': ax_underwater,
                'out': ax_out,
                'after': ax_after,
                }
            
    flyweight = 1
    calibrationmass = 0.041e-3 # 41 mg = .041 g = .041e-3 kg
    
    experiments_to_solutions = {'mono_concentration': ['double', 'half', 'normal'],
                                'salt_mono_base': ['salt', 'mono', 'base'],
                                'ephydra_2016_mono_sunscreen_160mg_2hands': ['mono_2016', 'sunscreen_2016'],
                                'ephydra_2016_mono_sunscreen_160mg_2hands_control': ['mono_2016', 'sunscreen_2016'],
                                'ephydra_2016_sunscreen_160mg_2hands_diwater': ['diwater_2016', 'sunscreen_2016'],
                                'ephydra_2016_mono_sunscreen_160mg_2hands_15min': ['mono_2016', 'sunscreen_2016'],
                                'fresh_diwater_monowater_hexane': ['diwater_2016', 'diwater_hexane', 'monowater_2016', 'monowater_hexane'],
                                'live_from_the_field_ephydra_acetonedips': ['water'],
                                'frozen_diwater_monowater_hexane': ['diwater_2016', 'monowater_2016', 'monowater_hexane'],
                                '_ph': ['bicarb', 'mix', 'carb'],
                                'diwater_monowater_carb': ['diwater_2016', 'monowater_2016', 'carb'],
                                'ephydra_2016_neutralized_monowater': ['monowater_2016', 'neutralizedmono_2016'],
                                'ephydra_2016_naoh': ['diwater_2016', 'naoh_2016'],
                                'ephydra_2016_co3_tests': ['neutralcarb', '3rdcarb', 'carb_2016'],
                                'calibration_units_check': ['check_2016',],
                               }
    experiments_to_colors = {'mono_concentration': {'double': (0.001, 0.001, 0.001),
                                                    'normal': 'blue',
                                                    'half': 'teal',
                                                    },
                             'salt_mono_base': {'salt': 'magenta',
                                                    'mono': 'blue',
                                                    'base': 'orangered',
                                                    },
                             '_ph': {            'bicarb': 'orange',
                                                'mix': 'orangered',
                                                'carb': 'red',
                                                    },
                             'ephydra_2016_mono_sunscreen_160mg_2hands': {'mono_2016': 'blue', 
                                                                          'sunscreen_2016': 'chocolate'},
                              'ephydra_2016_mono_sunscreen_160mg_2hands_15min': {'mono_2016': 'blue', 
                                                                          'sunscreen_2016': 'chocolate'},
                            'ephydra_2016_mono_sunscreen_160mg_2hands_control': {'mono_2016': 'blue', 
                                                                          'sunscreen_2016': 'indianred'},
                            'ephydra_2016_sunscreen_160mg_2hands_diwater': {'diwater_2016': 'green', 
                                                                          'sunscreen_2016': 'chocolate'},
                              'live_from_the_field_ephydra_acetonedips': {'water': 'blue'},
                             'fresh_diwater_monowater_hexane':    {'diwater_2016': 'green',
                                                             'monowater_2016': 'blue',
                                                             'diwater_hexane': 'orange',
                                                             'monowater_hexane': 'purple',
                                                            },
                            'frozen_diwater_monowater_hexane':    {'diwater_2016': 'green',
                                                             'monowater_2016': 'blue',
                                                             'monowater_hexane': 'purple',
                                                            },
                            'diwater_monowater_carb':    {  'diwater_2016': 'green',
                                                             'monowater_2016': 'blue',
                                                             'carb': 'red',
                                                            },
                            'ephydra_2016_neutralized_monowater': {  'monowater_2016': 'blue',
                                                             'neutralizedmono_2016': 'purple',
                                                            },
                            'ephydra_2016_naoh': {          'diwater_2016': 'green',
                                                            'naoh_2016': 'red',
                                                            },
                            'ephydra_2016_co3_tests': {     'neutralcarb': 'green',
                                                            '3rdcarb': 'orange',
                                                            'carb_2016': 'red',
                                                            },
                            'calibration_units_check': {'check_2016': 'red',},
                            }
    for key, solutions in experiments_to_solutions.items():
        if key in path: 
            print 'found experiment from path: ', key
            experiment_name = key
            break
    print experiment_name
        
    if flynum is None:
        if oddonly:
            numbers = [i for i in range(1,40,2)]
        elif evenonly:
            numbers = [i for i in range(2,40,2)]
        else:
            numbers = [i for i in range(1,40)]
    else:
        numbers = [flynum]
        
    basenames = []
    for solution in solutions:
        basenames.extend(['fly' + str(i) + '_' + solution for i in numbers])
    
    max_force = -1000
    min_force = 1000
    
    colors = []
    f_all = None
        
    summary_data_out = {solution: {'in': {}, 'out': {}} for solution in solutions}
    for solution in solutions:
        summary_data_out[solution]['distance_travelled'] = {'in': {}, 'out': {}}
    
    tags = ['before', 'in', 'underwater', 'out', 'after']
    interpolated_t = {tag:[] for tag in tags}
    interpolated_f = {tag:[] for tag in tags}
        
    for basename in basenames:  
        print basename  
        do_not_include = ['calibration_2016']
        if 'hexane' not in basename:
            do_not_include.append('hexane')
        filenames = get_all_filenames(path, basename, ['.bag'], do_not_include)
        print filenames
        if len(filenames) != 1:
            continue
        for filename in filenames:
            t, f, l = get_forces_for_plunge(filename, 0, calibrationmass)
        
        for solution in experiments_to_solutions[experiment_name]:
            if solution in os.path.basename(filename):
                break 
            
        color = experiments_to_colors[experiment_name][solution]
        
        for tag in tags:
            print basename, tag 
            print filenames
            print t[tag]
            print f[tag]
            print
            ax = axes[tag]
            
            #if len(basenames) > 1:
            linewidth=0.25
            #else:
            #    linewidth=1
            
            if solution_to_plot in basename:
                ax.plot(t[tag]/np.max(t[tag]), f[tag], color=color, linewidth=linewidth)
                if fill_work_out:
                    if tag == 'out':
                        ax.fill_between(t[tag]/np.max(t[tag]), f[tag], np.zeros_like(f[tag]), facecolor='blue', edgecolor='none', alpha=0.3)

                interp_t = np.linspace(0, 1, 100)
                interp_f = np.interp(interp_t, t[tag]/np.max(t[tag]), f[tag])
                interpolated_t[tag].append( interp_t)
                interpolated_f[tag].append( interp_f)
        
        flynum_ = basename.split('_')[0]
        
        if calc == 'work':
            results = scipy.stats.linregress(t['in'], l['in'])
            dip_speed_in = np.abs(results[0])
            results = scipy.stats.linregress(t['out'], l['out'])
            dip_speed_out = np.abs(results[0])
        
            total_time_out = t['out'][-1] - t['out'][0]
            seconds_per_frame_out = total_time_out / float(len(t['out']))
            total_time_in = t['in'][-1] - t['in'][0]
            seconds_per_frame_in = total_time_in / float(len(t['in']))
            
            distance_travelled = (l['out'][-1] - l['out'][0])*1e3
            
            summary_data_out[solution]['out'][flynum_] = np.sum( f['out']-np.mean(f['after']) )*dip_speed_out*1e3*seconds_per_frame_out # dip_speed_out*1e3 in mm / sec, so work is in mN/frame*mm/sec*sec/frame = uJ
            summary_data_out[solution]['in'][flynum_] = np.sum( f['in']-np.mean(f['before']) )*dip_speed_in*1e3*seconds_per_frame_in
            
            print summary_data_out[solution]['distance_travelled']['out']
            summary_data_out[solution]['distance_travelled']['out'][flynum_] = distance_travelled
            
            print 'DIP SPEED: ', dip_speed_out*1e3
            print 'F after: ', np.mean(f['after'])
            
        elif calc == 'peak':
            summary_data_out[solution]['out'][flynum_] = np.max( f['out'] )-np.mean(f['after']) 
            summary_data_out[solution]['in'][flynum_] = np.max( f['in'] )-np.mean(f['before']) 
        elif calc == 'meanforce':
            summary_data_out[solution]['out'][flynum_] = np.mean( f['out'] -np.mean(f['after']) )
            summary_data_out[solution]['in'][flynum_] = np.mean( f['in'] -np.mean(f['before']) )
            
    try:
        print len(interpolated_t)
        for tag in tags:
            ax = axes[tag]
            if flynum is None:
                t = np.mean(np.vstack(interpolated_t[tag]), axis=0)
                f = np.mean(np.vstack(interpolated_f[tag]), axis=0)
                color = experiments_to_colors[experiment_name][solution_to_plot]
                ax.plot(t, f, color=color, linewidth=2)
                
        for tag in tags:
            ax = axes[tag]
            if ylim is None:
                ax.set_ylim(-0.0001, 0.00125)
            else:
                ax.set_ylim(ylim[0], ylim[1])
            if make_ax_spines_clean:
                figurefirst.mpl_functions.adjust_spines(ax, [])
            else:
                if tag == 'before':
                    axes['before'].set_yticks([0,0.001])
                    figurefirst.mpl_functions.adjust_spines(ax, ['left'], yticks=[0,0.001], spine_locations={'left': 5})
                    ax.set_yticklabels(['0', '1'])
                    flytext.set_fontsize(ax.figure, 6) 
                else:
                    figurefirst.mpl_functions.adjust_spines(ax, [])
    except:
        print 'COULD NOT MAKE PLOT'
        pass
        
    return summary_data_out, experiments_to_solutions, experiments_to_colors
    
def full_concentration_plot(in_or_out='out', ax=None, baselinesubtract=False):
    sdo, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_fresh_diwater_monowater_hexane', solvent='')
    keys = sdo['diwater_2016'][in_or_out].keys()
    sdo_solution = np.array([sdo['diwater_2016'][in_or_out][key] for key in keys])
    sdo_monosolution = np.array([sdo['monowater_2016'][in_or_out][key] for key in keys])
    diwater = sdo_solution #- sdo_monosolution + np.mean(sdo_monosolution)
    if baselinesubtract:
        diwater = diwater - sdo_monosolution + np.mean(sdo_monosolution)
    monowater_1 = sdo_monosolution# - np.mean(sdo_monosolution) + np.mean(sdo_monosolution)
    if baselinesubtract:
        monowater_1 = monowater_1 - np.mean(sdo_monosolution) + np.mean(sdo_monosolution)
        
    sdo, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_concentration', solvent='')
    
    keys = sdo['normal'][in_or_out].keys()
    sdo_monosolution = np.array([sdo['normal'][in_or_out][key] for key in keys])
    half = np.array([sdo['half'][in_or_out][key] for key in keys]) #- sdo_monosolution + np.mean(sdo_monosolution)
    double = np.array([sdo['double'][in_or_out][key] for key in keys]) #- sdo_monosolution + np.mean(sdo_monosolution)
    monowater_2 = sdo_monosolution #- np.mean(sdo_monosolution) + np.mean(sdo_monosolution)
    
    if baselinesubtract:
        half = half - sdo_monosolution + np.mean(sdo_monosolution)
        double = double - sdo_monosolution + np.mean(sdo_monosolution)
        monowater_2 = monowater_2 - np.mean(sdo_monosolution) + np.mean(sdo_monosolution)
    
    monowater = np.hstack((monowater_1, monowater_2))
    
    LR_x = np.array([1]*len(diwater) + [2]*len(half) + [3]*len(monowater) + [4]*len(double))
    LR_y = np.hstack((diwater, half, monowater, double))
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    color='green'
    fpl.scatter_box(ax, 1, diwater, xwidth=0.4, ywidth=0.1, color=color, edgecolor=color, flipxy=False, shading='95conf', markersize=.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
    color='teal'
    fpl.scatter_box(ax, 2, half, xwidth=0.4, ywidth=0.1, color=color, edgecolor=color, flipxy=False, shading='95conf', markersize=.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
    color='blue'
    fpl.scatter_box(ax, 3, monowater, xwidth=0.4, ywidth=0.1, color=color, edgecolor=color, flipxy=False, shading='95conf', markersize=.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
    color=(0.001, 0.001, 0.001)
    fpl.scatter_box(ax, 4, double, xwidth=0.4, ywidth=0.1, color=color, edgecolor=color, flipxy=False, shading='95conf', markersize=.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
    
    r = scipy.stats.linregress(LR_x, LR_y)
    slope, intercept, rval, pval, stderr = r
    print 'slope, intercept, rval, pval, stderr'
    print slope, intercept, rval, pval, stderr
    
    if 1:#pval < 0.05:
        xs = np.linspace(np.min(LR_x)-1, np.max(LR_x)+1, 20)
        ys = slope*xs + intercept
        ax.plot(xs, ys, color='red',zorder=-10)
    
    
    # null hypothesis density line
    water_densities = [1, 1.03, 1.07, 1.14]
    estimated_forces = np.array(water_densities)*np.mean(diwater)
    r_null = scipy.stats.linregress(water_densities, estimated_forces)
    slope, intercept, rval, pval, stderr = r_null
    
    xs = [1,2,3,4]
    r_xs = scipy.stats.linregress(xs, water_densities)
    slope_xs, intercept_xs, rval_xs, pval_xs, stderr_xs = r_xs
    
    xs = [0,1,2,3,4,5]
    ys = [(x*slope_xs + intercept_xs)*slope + intercept for x in xs]
    
    ax.plot(xs, ys, color='black',zorder=-10)
    
    
    # now calculate corrected r; subtracting out the density effect
    LR_x = np.array([1]*len(diwater) + [2]*len(half) + [3]*len(monowater) + [4]*len(double))
    LR_y = np.hstack((diwater-estimated_forces[0], half-estimated_forces[1], monowater-estimated_forces[2], double-estimated_forces[3]))
    
    r = scipy.stats.linregress(LR_x, LR_y)
    slope, intercept, rval, pval, stderr = r
    print 'slope, intercept, rval, pval, stderr'
    print slope, intercept, rval, pval, stderr
    
    return r
    

def scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=None, show_regression=True):
    sdo, experiments_to_solutions, experiments_to_colors = make_figure(path, solvent='', oddonly=oddonly, evenonly=evenonly)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    
    for key, solutions in experiments_to_solutions.items():
        if key in path: 
            print 'found experiment from path: ', key
            experiment_name = key
    print experiment_name
    
    experiment_name_to_mono_solution = {'mono_concentration': 'normal',
                                        'fresh_diwater_monowater_hexane': 'monowater_2016',
                                        'frozen_diwater_monowater_hexane': 'monowater_2016',
                                        'salt_mono_base': 'mono',
                                        '_ph': 'mix',
                                        'diwater_monowater_carb': 'monowater_2016',
                                        'ephydra_2016_mono_sunscreen_160mg_2hands': 'mono_2016',
                                        'ephydra_2016_mono_sunscreen_160mg_2hands_15min': 'mono_2016',
                                        'ephydra_2016_mono_sunscreen_160mg_2hands_control': 'mono_2016',
                                        'ephydra_2016_sunscreen_160mg_2hands_diwater': 'diwater_2016',
                                        'ephydra_2016_neutralized_monowater': 'monowater_2016',
                                        'ephydra_2016_naoh': 'diwater_2016',
                                        'ephydra_2016_co3_tests': '3rdcarb',
                                        'calibration_units_check': 'check_2016',
                                        
                                        }
    monosolution = experiment_name_to_mono_solution[experiment_name]
    
    n = 0
    
    LR_x = []
    LR_y = []
    
    solution_data = {}
    for solution in experiments_to_solutions[experiment_name]:
        n += 1
        color = experiments_to_colors[experiment_name][solution]
        
        keys = sdo[solution][in_or_out].keys()
        
        sdo_solution = np.array([sdo[solution][in_or_out][key] for key in keys])
        sdo_monosolution = np.array([sdo[monosolution][in_or_out][key] for key in keys])
        
        if baselinesubtract:
            if solution != monosolution:        
                sdo_val = sdo_solution - sdo_monosolution
            else:
                sdo_val = sdo_solution - np.mean(sdo_solution)
            
            sdo_val += np.mean(sdo[monosolution][in_or_out].values()) # add the total mono solution mean back in
            
        else:
            sdo_val = sdo_solution
            
        if experiment_name == 'frozen_diwater_monowater_hexane' and solution == 'monowater_2016':
            sdo_val_stuck = sdo_val[sdo_val<0]
            fpl.scatter_box(ax, n, sdo_val_stuck, xwidth=0.4, ywidth=0.1, color=(0.001, 0.001, 0.001), edgecolor=(0.001, 0.001, 0.001), flipxy=False, shading='95conf', markersize=0.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
            
            sdo_val_notstuck = sdo_val[sdo_val>0]
            fpl.scatter_box(ax, n, sdo_val_notstuck, xwidth=0.4, ywidth=0.1, color=color, edgecolor=color, flipxy=False, shading='95conf', markersize=0.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
            
        else:
            fpl.scatter_box(ax, n, sdo_val, xwidth=0.4, ywidth=0.1, color=color, edgecolor=color, flipxy=False, shading='95conf', markersize=0.7, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1)
    
        LR_y.extend(sdo_val.tolist())
        LR_x.extend( (n*np.ones_like(sdo_val)).tolist() )
        
        solution_data[solution] = sdo_val
    
    r = scipy.stats.linregress(LR_x, LR_y)
    slope, intercept, rval, pval, stderr = r
    print 'slope, intercept, rval, pval, stderr'
    print slope, intercept, rval, pval, stderr
    
    if 'hexane' not in path and 'diwater_monowater_carb' not in experiment_name:
        if pval < 0.05 and show_regression:
            xs = np.linspace(np.min(LR_x), np.max(LR_x), 20)
            ys = slope*xs + intercept
            ax.plot(xs, ys, color='red',zorder=-10)
    
    
    return r, solution_data
    
def get_forces_for_plunge(filename, objectmass, calibrationmass):
    fname = os.path.splitext(filename)[0]
    output_fname = fname + '.hdf5'
           
    try: 
        if not os.path.exists(output_fname):
            if not os.path.exists(filename):
                print >> sys.stderr, 'No file %s' % filename
                sys.exit(1)
            bag2hdf5.bag2hdf5(filename,
                     output_fname,
                     max_strlen=200,
                     topics=['/phidgets_daq/force','/phidgets_daq/lvdt'])  
    except:
        if not os.path.exists(output_fname):
            if not os.path.exists(filename):
                print >> sys.stderr, 'No file %s' % filename
                sys.exit(1)
            bag2hdf5.bag2hdf5(filename,
                     output_fname,
                     max_strlen=200,
                     topics=['/phidgets_daq/force'])            
        
    
    t, force, lvdt = get_time_and_force_traces(output_fname)  
    t -= t[0]
    
    calibration = get_calibration(filename, calibrationmass)

    try:
        metadata_filename = fname +  '_calibration.pickle'
        metadata = get_metadata(metadata_filename)
    except:
        metadata_filename = fname +  '_metadata.pickle'
        metadata = get_metadata(metadata_filename)
    
    indices_before = [metadata['1']['indices'][0]-100, metadata['1']['indices'][0]]
    indices_going_in = metadata['1']['indices']
    indices_underwater = metadata['2']['indices']
    indices_coming_out = metadata['3']['indices']
    indices_after = [metadata['3']['indices'][-1], metadata['3']['indices'][-1]+100]

    start = np.mean(force[indices_before[0]:indices_before[1]])
    underwater = np.mean(force[indices_underwater[0]:indices_underwater[1]])
    finish = np.mean(force[indices_after[0]:indices_after[1]])

    force_N = ((force - start)*calibration)# - objectmass*9.81
    
    # Calibration
    
    #
    
    # before
    t_before = t[indices_before[0]: indices_before[-1]]
    t_before -= t_before[0]
    if lvdt is not None:
        lvdt_before = lvdt[indices_before[0]: indices_before[-1]]
        lvdt_before -= lvdt_before[0]
        lvdt_before /= 360218.97262842144/4.
    else:
        lvdt_before = t_before / (np.mean(np.diff(t_before)))
        lvdt_before *= 7.2e-7
    force_before = force_N[indices_before[0]: indices_before[-1]]
    
    # in
    t_in = t[indices_going_in[0]: indices_going_in[-1]]
    t_in -= t_in[0]
    if lvdt is not None:
        lvdt_in = lvdt[indices_going_in[0]: indices_going_in[-1]]
        lvdt_in -= lvdt_in[0]
        lvdt_in /= 360218.97262842144/4.
    else:
        lvdt_in = t_in / (np.mean(np.diff(t_in)))
        lvdt_in *= 7.2e-7
    force_in = force_N[indices_going_in[0]: indices_going_in[-1]]
    
    # out
    t_out = t[indices_coming_out[0]: indices_coming_out[-1]]
    t_out -= t_out[0]
    if lvdt is not None:
        lvdt_out = lvdt[indices_coming_out[0]: indices_coming_out[-1]]
        lvdt_out -= lvdt_out[0]
        lvdt_out /= 360218.97262842144/4.
        lvdt_out *= -1
    else:
        lvdt_out = t_out / (np.mean(np.diff(t_out)))
        lvdt_out *= 7.2e-7
    force_out = force_N[indices_coming_out[0]: indices_coming_out[-1]]
    
    # underwater
    t_underwater = t[indices_underwater[0]: indices_underwater[0]+100]
    t_underwater -= t_underwater[0]
    if lvdt is not None:
        lvdt_underwater = lvdt[indices_underwater[0]: indices_underwater[0]+100]
        lvdt_underwater -= lvdt_underwater[0]
        lvdt_underwater /= 360218.97262842144/4.
    else:
        lvdt_underwater = t_underwater / (np.mean(np.diff(t_underwater)))
        lvdt_underwater *= 7.2e-3
    force_underwater = force_N[indices_underwater[0]: indices_underwater[0]+100]

    # after
    t_after = t[indices_after[0]: indices_after[-1]]
    t_after -= t_after[0]
    if lvdt is not None:
        lvdt_after = lvdt[indices_after[0]: indices_after[-1]]
        lvdt_after -= lvdt_after[0]
        lvdt_after /= 360218.97262842144/4.
    else:
        lvdt_after = t_after / (np.mean(np.diff(t_after)))
        lvdt_after *= 7.2e-7
    force_after = force_N[indices_after[0]: indices_after[-1]]
    
    # work
    #work_in = np.sum( (force_in_total[1:] + objectmass*9.81)*np.diff(lvdt_in) ) 
    #work_out = np.sum( (force_out_total[1:] + objectmass*9.81)*np.diff(lvdt_out) ) 
    
    #water_stuck_to_fly = (start - finish)*calibration
    #bouyant_force = (underwater - start)*calibration
    
    times = {'before': t_before,
             'in': t_in,
             'out': t_out,
             'underwater': t_underwater,
             'after': t_after,
             }
    
    forces = {  'before': force_before,
                'in': force_in,
                'out': force_out,
                'underwater': force_underwater,
                'after': force_after,
                }
    
    lvdt = {    'before': lvdt_before,
                'in': lvdt_in,
                'out': lvdt_out,
                'underwater': lvdt_underwater,
                'after': lvdt_after,
                }
                
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(force_N)
    #ax.plot(forces['before'])
             
    return times, forces, lvdt
    
def get_time_and_force_traces(filename):
    f = h5py.File(filename, 'r')  
            
    t = f['phidgets_daq']['force']['t']
    force = f['phidgets_daq']['force']['value']
    
    try:
        lvdt = f['phidgets_daq']['lvdt']['value']
    except:
        lvdt = None
        
    a,b = scipy.signal.filter_design.butter(3,0.01)
    force_filtered = scipy.signal.filtfilt(a,b,force)

    return t, force_filtered, lvdt
    
def get_metadata(filename):
    if not os.path.exists(filename):
        print >> sys.stderr, 'No file %s' % filename
        sys.exit(1)
    
    f = open(filename)
    metadata = pickle.load(f)
    
    return metadata
    
def get_calibration(filename_of_data_bag, calibrationmass, splitstr='_'):
    path = os.path.dirname(filename_of_data_bag)
    identifier = os.path.basename(filename_of_data_bag).split(splitstr)[0]
    calibration_filename_startswith = identifier + '_calibration'
    possible_filenames = get_all_filenames(path, calibration_filename_startswith, ['.pickle'], '~')
    if len(possible_filenames) != 1:
        raise ValueError('too many, or too few, possible filenames')
    
    calibration_filename = possible_filenames[0]
    
    
    calibration_filename = os.path.join(path, calibration_filename)
    cal = open(calibration_filename)
    calibration = pickle.load(cal)
    
    factor1 = calibration['1']['mean'] - calibration['2']['mean']
    factor2 = calibration['3']['mean'] - calibration['2']['mean']
    factor = np.mean([factor1, factor2])
    
    # F = ma
    weight = calibrationmass*9.81 
    gain = weight / factor
    
    return gain
    
def get_all_filenames(basepath, startswith, contains, does_not_contain):
    
    cmd = 'ls ' + basepath
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass

    filelist = []
    for i, filename in enumerate(all_filelist):
        fileisgood = True
        if filename.startswith(startswith):
            for c in contains:
                if c not in filename:
                    fileisgood = False
            for d in does_not_contain:
                if d in filename:
                    fileisgood = False
            if fileisgood:
                filelist.append( os.path.join(basepath, filename) )
    
    return filelist
    

def plot_sunscreen_mortality(ax):
    
    mg_neutrogena = [0,2,8,20,43,97]
    n_flies_stuck_0 = np.array([1,2,3,9,10,10]) # out of 10 
    n_flies_stuck_1 = np.array([0,0,5,9,10,10])
    n_flies_stuck_2 = np.array([0,0,6,10,10,10])
    n_flies_stuck = n_flies_stuck_0 + n_flies_stuck_1 + n_flies_stuck_2
    
    n_flies_stuck_expanded = []
    for n in n_flies_stuck:
        a = np.array([1]*n + [0]*(30-n))
        n_flies_stuck_expanded.append(a)
        
    nth_data_pt = 0
    for i in range(len(mg_neutrogena)):
        nth_data_pt += 1
        fpl.scatter_box(ax, nth_data_pt, n_flies_stuck_expanded[i], xwidth=0.4, ywidth=0.1, color='chocolate', edgecolor='none', flipxy=False, shading='95conf', markersize=1, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)

    if 0:
        LR_x = [[mgn]*10 for mgn in mg_neutrogena]
        LR_x = np.hstack(LR_x)
        LR_y = np.hstack(n_flies_stuck_expanded)
        
        r = scipy.stats.linregress(LR_x, LR_y)
        slope, intercept, rval, pval, stderr = r
        print 'slope, intercept, rval, pval, stderr'
        print slope, intercept, rval, pval, stderr
        if pval < 0.05:
            xs = np.linspace(np.min(LR_x), np.max(LR_x), 20)
            ys = slope*xs + intercept
            ax.plot(np.log(xs), ys, color='red',zorder=-10)

def plot_sunscreen_type_mortality(ax):
    
    sunscreens = ['neutrogena', 'banana_spray', 'baby', 'badger', 'copper', 'bullfrog', 'control']
    n_flies_stuck_0 = np.array([6, 5, 2, 2, 0, 0, 0])
    n_flies_stuck_1 = np.array([8, 6, 0, 0, 0, 0, 1])
    n_flies_stuck_2 = np.array([9, 3, 3, 0, 0, 0, 0])
    n_flies_stuck = n_flies_stuck_0 + n_flies_stuck_1 + n_flies_stuck_2
                    
    n_flies_stuck_expanded = []
    for n in n_flies_stuck:
        a = np.array([1]*n + [0]*(30-n))
        n_flies_stuck_expanded.append(a)
        
    for i, name in enumerate(sunscreens):
        fpl.scatter_box(ax, i+1, n_flies_stuck_expanded[i], xwidth=0.4, ywidth=0.1, color=(0.001, 0.001, 0.001), edgecolor='none', flipxy=False, shading='95conf', markersize=1, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)


def plot_dimethicone_mortality(ax):
    
    volumes = [0, 2, 10] # uL
    n_flies_stuck_0 = np.array([1,10,10])
    n_flies_stuck_1 = np.array([0,5,10])
    n_flies_stuck_2 = np.array([2,4,9])
    n_flies_stuck = n_flies_stuck_0 + n_flies_stuck_1 + n_flies_stuck_2
                    
    n_flies_stuck_expanded = []
    for n in n_flies_stuck:
        a = np.array([1]*n + [0]*(30-n))
        n_flies_stuck_expanded.append(a)
        
    for i, name in enumerate(volumes):
        fpl.scatter_box(ax, i+1, n_flies_stuck_expanded[i], xwidth=0.4, ywidth=0.1, color='lime', edgecolor='none', flipxy=False, shading='95conf', markersize=1, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)



############################

def get_paper_layout(figurename):
    svg = figurename
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

'''
def make_concentration_plot():
    layout = get_paper_layout()
    
    ax_out = layout.axes[('concentration', 'out')]
    ax_out.set_ylim(-2.5, 2.5)
    full_concentration_plot(in_or_out='out', ax=ax_out)
    figurefirst.mpl_functions.adjust_spines(ax_out, ['left'], yticks=[-1,0,1])
    
    ax_in = layout.axes[('concentration', 'in')]
    ax_in.set_ylim(0, 5)
    full_concentration_plot(in_or_out='in', ax=ax_in)
    figurefirst.mpl_functions.adjust_spines(ax_in, ['left'], yticks=[0, 5])
    
    flytext.set_fontsize(ax_in.figure, 10)
    
    layout.append_figure_to_layer(layout.figures['concentration'], 'concentration', cleartarget=True)
    layout.write_svg('figure_output.svg' )
    

def make_ph_plot():
    layout = get_paper_layout()
    
    ax_out = layout.axes[('ph', 'out')]
    ax_out.set_ylim(-2.5, 2.5)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_ph'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=True, in_or_out='out', ax=ax_out)
    figurefirst.mpl_functions.adjust_spines(ax_out, ['left'], yticks=[-1,0,1])
    
    ax_in = layout.axes[('ph', 'in')]
    ax_in.set_ylim(0, 5)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_ph'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=True, in_or_out='in', ax=ax_in)
    figurefirst.mpl_functions.adjust_spines(ax_in, ['left'], yticks=[0, 5])
    
    flytext.set_fontsize(ax_in.figure, 10)
    
    layout.append_figure_to_layer(layout.figures['ph'], 'ph', cleartarget=True)
    layout.write_svg('figure_output.svg' )
    
def make_salt_mono_base_plot():
    layout = get_paper_layout()
    
    ax_out = layout.axes[('salt_mono_base', 'out')]
    ax_out.set_ylim(-2.5, 2.5)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_salt_mono_base'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=True, in_or_out='out', ax=ax_out)
    figurefirst.mpl_functions.adjust_spines(ax_out, ['left'], yticks=[-1,0,1])
    
    ax_in = layout.axes[('salt_mono_base', 'in')]
    ax_in.set_ylim(0, 5)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_salt_mono_base'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=True, in_or_out='in', ax=ax_in)
    figurefirst.mpl_functions.adjust_spines(ax_in, ['left'], yticks=[0, 5])
    
    flytext.set_fontsize(ax_in.figure, 10)
    
    layout.append_figure_to_layer(layout.figures['salt_mono_base'], 'salt_mono_base', cleartarget=True)
    layout.write_svg('figure_output.svg' )
'''

def make_hexane_plot(baselinesubtract=False):
    layout = get_paper_layout()
    
    ax_out = layout.axes[('hexane', 'out')]
    ax_out.set_ylim(-2, 2)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_fresh_diwater_monowater_hexane'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out)
    figurefirst.mpl_functions.adjust_spines(ax_out, ['left'], yticks=[-1,0,1])
    
    ax_in = layout.axes[('hexane', 'in')]
    ax_in.set_ylim(0, 5)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_fresh_diwater_monowater_hexane'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='in', ax=ax_in)
    figurefirst.mpl_functions.adjust_spines(ax_in, ['left'], yticks=[0, 5])
    
    flytext.set_fontsize(ax_in.figure, 10)
    
    layout.append_figure_to_layer(layout.figures['hexane'], 'hexane', cleartarget=True)
    layout.write_svg('figure_output.svg' )
    
def make_species_di_mono_carb_plot():
    layout = get_paper_layout(mono_paper_locations.figure_template_locations.figure3)
    
    if 1:
        scaling = {'ephydra': 6,      # approximate millimeters, length
                   'blue_kelp': 6,
                   'black_kelp': 7,
                   'oil': 2.5,
                   'virilis': 3.5,
                   'melanogaster': 2.5,
                   'santa_ana_sp': 1.5,
                   }
                   
    else:
        scaling = {'ephydra': 1,      
                   'blue_kelp': 1,
                   'black_kelp': 1,
                   'oil': 1,
                   'virilis': 1,
                   'melanogaster': 1,
                   'santa_ana_sp': 1,
                   }
        
    def make_scale_bar(ax, layout, species):
        ymax = ax.get_ylim()[1]
        ymid = ymax/2.
        ax.vlines(0.01,0+ymid*.1, 0+ymid*.1+ymid, color='black',linewidth=0.5)
        svgitem = 'yscale_' + species
        layout.svgitems[svgitem].style['font-size'] = 6
        d = (ymax-ymid)*1e3
        layout.svgitems[svgitem].text = '{0:.2f}'.format(d) + ' uJ'
        
    scale_power_factor = 1
    ylim = [-.0015, .002]
    
    ax = layout.axes[('species', 'black_kelp')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/blackkelp_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['black_kelp']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['black_kelp']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'black_kelp')
    
    ax = layout.axes[('species', 'melanogaster')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/melanogaster_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['melanogaster']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['melanogaster']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'melanogaster')
    
    ax = layout.axes[('species', 'ephydra')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['ephydra']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['ephydra']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'ephydra')
    
    ax = layout.axes[('species', 'virilis')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/virilis_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['virilis']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['virilis']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'virilis')
    
    ax = layout.axes[('species', 'blue_kelp')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/bluekelp_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['blue_kelp']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['blue_kelp']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'blue_kelp')
    
    ax = layout.axes[('species', 'oil')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/oilflies_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['oil']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['oil']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'oil')
    
    ax = layout.axes[('species', 'santa_ana_sp')]
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/santa_ana_2016_diwater_monowater_carb'
    scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax)
    ax.set_ylim(ylim[0]*scaling['santa_ana_sp']**scale_power_factor/float(6**scale_power_factor), ylim[1]*scaling['santa_ana_sp']**scale_power_factor/float(6**scale_power_factor))
    figurefirst.mpl_functions.adjust_spines(ax, [])
    ax.hlines(0, 0, 4, linewidth=0.5, color='black')
    make_scale_bar(ax, layout, 'santa_ana_sp')
    
    flytext.set_fontsize(ax.figure, 10)
    
    layout.append_figure_to_layer(layout.figures['species'], 'species', cleartarget=True)
    layout.apply_svg_attrs()
    layout.write_svg(mono_paper_locations.figure_template_locations.figure3 )
    
def make_frozen_hexane_plot(baselinesubtract=False):
    layout = get_paper_layout(mono_paper_locations.figure_template_locations.figure3_chc)
    
    if 0:
        ax_out = layout.axes[('hexane_frozen', 'frozen')]
        ax_out.set_ylim(-.002, .002)
        ax_out.set_xlim(0,4)
        path = mono_paper_locations.data_locations.frozen_hexane_data
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        ax_out.hlines(0, 0, 5, linewidth=0.5, color='black')
    
    ax_out = layout.axes[('hexane_frozen', 'hexane')]
    ax_out.set_ylim(-.002, .002)
    ax_out.set_xlim(0,5)
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_fresh_diwater_monowater_hexane'
    r, solution_data = scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out)
    figurefirst.mpl_functions.adjust_spines(ax_out, [])
    #figurefirst.mpl_functions.adjust_spines(ax_out, ['left'], xticks=[0,0.001])
    ax_out.hlines(0, 0, 5, linewidth=0.5, color='black')

    print solution_data.keys()
    if 1:
        print
        pval, actual_diff, diffs = permutation_test(solution_data['monowater_2016'], solution_data['monowater_hexane'], permutations=1000)
        print pval
        if pval < 0.001:
            pval = 0.001
        layout.svgitems['pval_hexane_mono'].style['font-size'] = 6
        layout.svgitems['pval_hexane_mono'].text = 'p={0:.3f}'.format(pval)

        pval, actual_diff, diffs = permutation_test(solution_data['diwater_2016'], solution_data['diwater_hexane'], permutations=1000)
        print pval
        if pval < 0.001:
            pval = 0.001
        layout.svgitems['pval_hexane_diwater'].style['font-size'] = 6
        layout.svgitems['pval_hexane_diwater'].text = 'p={0:.3f}'.format(pval)

    layout.apply_svg_attrs(['pval_hexane_mono', 'pval_hexane_diwater'])
    layout.append_figure_to_layer(layout.figures['hexane_frozen'], 'hexane_frozen', cleartarget=True)
    layout.write_svg(mono_paper_locations.figure_template_locations.figure3_chc )
    
    
def permutation_test(A, B, permutations=1000):
    mixture = np.hstack((A,B))
    diffs = []
    for p in range(permutations):
        np.random.shuffle(mixture)
        A_p = mixture[0:len(A)]
        B_p = mixture[len(A):len(A)+len(B)]
        A_m = np.mean(A_p)
        B_m = np.mean(B_p)
        diffs.append(A_m-B_m)
    diffs.sort()
    diffs = np.array(diffs)
    actual_diff = np.mean(A)-np.mean(B)
    idx = np.argmin( np.abs(diffs-actual_diff))
    L = len(diffs)
    pval = (L/2. - np.abs(L/2. - idx))/L*2
    return pval, actual_diff, diffs

def update_pvals_for_figure_2(layout, baselinesubtract=True):
    
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_salt_mono_base'
    r, solution_data = scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=None, show_regression=False)
    pval, actual_diff, diffs = permutation_test(solution_data['salt'], solution_data['base'], permutations=1000)
    if pval < 0.001:
        pval = 0.001
    layout.svgitems['pval_salt_base'].style['font-size'] = 6
    layout.svgitems['pval_salt_base'].text = 'p<{0:.3f}'.format(pval)
    
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_ph'
    r, solution_data = scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=None, show_regression=False)
    pval, actual_diff, diffs = permutation_test(solution_data['bicarb'], solution_data['carb'], permutations=1000)
    if pval < 0.001:
        pval = 0.001
    layout.svgitems['pval_carbonate'].style['font-size'] = 6
    layout.svgitems['pval_carbonate'].text = 'p<{0:.3f}'.format(pval)
    
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_co3_tests'
    r, solution_data = scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=None, show_regression=False)
    pval, actual_diff, diffs = permutation_test(solution_data['neutralcarb'], solution_data['carb_2016'], permutations=1000)
    if pval < 0.001:
        pval = 0.001
    layout.svgitems['pval_neutral_carbonate'].style['font-size'] = 6
    layout.svgitems['pval_neutral_carbonate'].text = 'p<{0:.3f}'.format(pval)
    
    path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_co3_tests'
    r, solution_data = scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=None, show_regression=False)
    pval, actual_diff, diffs = permutation_test(solution_data['3rdcarb'], solution_data['carb_2016'], permutations=1000)
    if pval < 0.001:
        pval = 0.001
    layout.svgitems['pval_carbonate_concentration'].style['font-size'] = 6
    layout.svgitems['pval_carbonate_concentration'].text = 'p<{0:.3f}'.format(pval)
    
    #layout.apply_svg_attrs()
    #layout.write_svg(mono_paper_locations.figure_template_locations.figure2 )
    
    
def make_figure_2(baselinesubtract=False):
    layout = get_paper_layout(mono_paper_locations.figure_template_locations.figure2)
    
    axes = {}
    tags = ['before', 'in', 'underwater', 'out', 'after']
    for tag in tags:
        axes[tag] = layout.axes[('force_trace', tag)]
    
    summary_data_out, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_concentration', calc='work', axes=axes, solution_to_plot='normal', make_ax_spines_clean=False )
    layout.append_figure_to_layer(layout.figures['force_trace'], 'force_trace', cleartarget=True)
    
    if 1:
        ax_out = layout.axes[('concentration', 'concentration')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0,5)
        r = full_concentration_plot(in_or_out='out', ax=ax_out, baselinesubtract=baselinesubtract)
        slope, intercept, rval, pval, stderr = r
        ax_out.hlines(0, 0, 5, linewidth=0.5, color='black') 
        ax_out.set_xlim(0,5)
        #figurefirst.mpl_functions.adjust_spines(ax_out, [])
        figurefirst.mpl_functions.adjust_spines(ax_out, ['left'], yticks=[0,0.002], spine_locations={'left': 5})
        ax_out.set_yticklabels(['0', '2'])
        flytext.set_fontsize(ax_out.figure, 6) 
        layout.append_figure_to_layer(layout.figures['concentration'], 'concentration', cleartarget=True)
        slope, intercept, rval, pval, stderr = r
        layout.svgitems['pval_concentration_regression'].style['font-size'] = 6
        if pval < 0.001:
            pval = 0.001
            layout.svgitems['pval_concentration_regression'].text = 'p<{0:.3f}'.format(pval) + ', r2={0:.2f}'.format(rval**2)
        else:
            layout.svgitems['pval_concentration_regression'].text = 'p={0:.3f}'.format(pval) + ', r2={0:.2f}'.format(rval**2)
            
        ax_out = layout.axes[('salt_mono_base', 'salt_mono_base')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0,4)
        ax_out.hlines(0, 0, 4, linewidth=0.5, color='black')
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_salt_mono_base'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        layout.append_figure_to_layer(layout.figures['salt_mono_base'], 'salt_mono_base', cleartarget=True)
        
        ax_out = layout.axes[('ph', 'ph')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0,4)
        ax_out.hlines(0, 0, 4, linewidth=0.5, color='black')
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_ph'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        layout.append_figure_to_layer(layout.figures['ph'], 'ph', cleartarget=True)
        
        ax_out = layout.axes[('ph_odd', 'ph_odd')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0,4)
        ax_out.hlines(0, 0, 4, linewidth=0.5, color='black')
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_ph'
        scatter_plot_summary_data_out(path, oddonly=True, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        layout.append_figure_to_layer(layout.figures['ph_odd'], 'ph_odd', cleartarget=True)

        ax_out = layout.axes[('neutral_mono', 'neutral_mono')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0, 3)
        ax_out.hlines(0, 0, 4, linewidth=0.5, color='black')
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_neutralized_monowater'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        layout.append_figure_to_layer(layout.figures['neutral_mono'], 'neutral_mono', cleartarget=True)
        
        ax_out = layout.axes[('naoh', 'naoh')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0, 3)
        ax_out.hlines(0, 0, 4, linewidth=0.5, color='black')
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_naoh'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        layout.append_figure_to_layer(layout.figures['naoh'], 'naoh', cleartarget=True)

        ax_out = layout.axes[('co3_tests', 'co3_tests')]
        ax_out.set_ylim(-.001, .002)
        ax_out.set_xlim(0,4)
        ax_out.hlines(0, 0, 4, linewidth=0.5, color='black')
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_co3_tests'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        layout.append_figure_to_layer(layout.figures['co3_tests'], 'co3_tests', cleartarget=True)
        
    
    
    # single force trace
    axes = {}
    tags = ['before', 'in', 'underwater', 'out', 'after']
    for tag in tags:
        axes[tag] = layout.axes[('work_example', tag)]
    
    flynum = 10
    
    summary_data_out, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_concentration', calc='work', axes=axes, solution_to_plot='normal', flynum=flynum, fill_work_out=True, make_ax_spines_clean=True)
    layout.append_figure_to_layer(layout.figures['work_example'], 'work_example', cleartarget=True)
    
    svgitem = 'distance_travelled'
    layout.svgitems[svgitem].style['font-size'] = 8
    d = summary_data_out['normal']['distance_travelled']['out']['fly10']
    layout.svgitems[svgitem].text = '{0:.2f}'.format(d) + ' mm'
    
    update_pvals_for_figure_2(layout, baselinesubtract=baselinesubtract)
    
    layout.apply_svg_attrs()
    layout.write_svg(mono_paper_locations.figure_template_locations.figure2 )
    
def make_angle_plot_for_field_flies():

    layout = get_paper_layout('figure_s1.svg')
    
    axes = {}
    tags = ['before', 'in', 'underwater', 'out', 'after']
    for tag in tags:
        axes[tag] = layout.axes[('field_exps', tag)]
    
    summary_data_out, experiments_to_solutions, experiments_to_colors = make_figure('/media/Orchard2/monolake/liveflies/live_from_the_field_ephydra_acetonedips', calc='peak', axes=axes, solution_to_plot='water', make_ax_spines_clean=True, ylim=[-.0006, .00125] )
    
    fly_angles_20150729 = {  'fly2': 10,
                             'fly3': 45,
                             'fly4': 30,
                             'fly5': 100,
                             'fly6': 50,
                             'fly8': 120,
                             'fly9': 45,
                             'fly10': 15,
                             'fly11': 120,
                             'fly12': 45,
                             'fly13': 90,
                             'fly14': 100,
                             'fly16': 100,
                             'fly17': 5,
                             'fly18': 5,
                             'fly19': 40,
                             'fly20': 135,
                             'fly21': 100,
                             }
    
    
    ax = layout.axes[('field_exps', 'angle_in')]
    keys = fly_angles_20150729.keys()
    ax.plot( [fly_angles_20150729[key] for key in keys], [summary_data_out['water']['in'][key] for key in keys], 'o', markersize=1, markeredgecolor='none', markerfacecolor='blue' )
    LR_x, LR_y = [fly_angles_20150729[key] for key in keys], [summary_data_out['water']['in'][key] for key in keys]
    r = scipy.stats.linregress(LR_x, LR_y)
    slope, intercept, rval, pval, stderr = r
    print 'slope, intercept, rval, pval, stderr'
    print slope, intercept, rval, pval, stderr
    if pval < 0.05:
        xs = np.linspace(np.min(LR_x), np.max(LR_x), 20)
        ys = slope*xs + intercept
        ax.plot(xs, ys, color='red',zorder=-10)
    ax.set_ylim(-0.0005, 0.00125)
    ax.set_xlim(0,150)
    figurefirst.mpl_functions.adjust_spines(ax, [])
    
    if 1:
        ax = layout.axes[('field_exps', 'angle_out')]
        keys = fly_angles_20150729.keys()
        a = np.array([summary_data_out['water']['out'][key] for key in keys])
        ax.plot( [fly_angles_20150729[key] for key in keys], a, 'o', markersize=1, markeredgecolor='none', markerfacecolor='blue' )
        LR_x, LR_y = [fly_angles_20150729[key] for key in keys], [summary_data_out['water']['out'][key] for key in keys]
        r = scipy.stats.linregress(LR_x, LR_y)
        slope, intercept, rval, pval, stderr = r
        print 'slope, intercept, rval, pval, stderr'
        print slope, intercept, rval, pval, stderr
        if pval < 0.05:
            xs = np.linspace(np.min(LR_x), np.max(LR_x), 20)
            ys = slope*xs + intercept
            ax.plot(xs, ys, color='red',zorder=-10)
        ax.set_ylim(-0.0005, 0.0012)
        ax.set_xlim(0,150)
        figurefirst.mpl_functions.adjust_spines(ax, [])
    
    
    layout.append_figure_to_layer(layout.figures['field_exps'], 'field_exps', cleartarget=True)
    layout.write_svg('figure_s1.svg' )
    
    
def make_figure_5(baselinesubtract=False, redo_wordcloud=False):
    layout = get_paper_layout(mono_paper_locations.figure_template_locations.figure5)
    
    axes = {}
    tags = ['before', 'in', 'underwater', 'out', 'after']
    for tag in tags:
        axes[tag] = layout.axes[('force_trace', tag)]
    summary_data_out, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_sunscreen_160mg_2hands', calc='work', axes=axes, solution_to_plot='mono_2016', ylim=None, make_ax_spines_clean=True )
    layout.append_figure_to_layer(layout.figures['force_trace'], 'force_trace', cleartarget=True)

    axes = {}
    tags = ['before', 'in', 'underwater', 'out', 'after']
    for tag in tags:
        axes[tag] = layout.axes[('sunscreen', tag)]
    summary_data_out, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_sunscreen_160mg_2hands', calc='work', axes=axes, solution_to_plot='sunscreen_2016', ylim=[-0.0005,0.00135-0.0005], make_ax_spines_clean=True )
    layout.append_figure_to_layer(layout.figures['sunscreen'], 'sunscreen', cleartarget=True)

    if 1:
        ax_out = layout.axes[('neutrogena_work', 'neutrogena_work_mono')]
        ax_out.set_ylim(-.001, .002)
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_sunscreen_160mg_2hands'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        ax_out.hlines(0, 0, 5, linewidth=0.5, color='black') 
        ax_out.set_xlim(0,3)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        
        ax_out = layout.axes[('neutrogena_work', 'neutrogena_work_control')]
        ax_out.set_ylim(-.001, .002)
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_sunscreen_160mg_2hands_control'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        ax_out.hlines(0, 0, 5, linewidth=0.5, color='black') 
        ax_out.set_xlim(0,3)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        
        ax_out = layout.axes[('neutrogena_work', 'neutrogena_work_diwater')]
        ax_out.set_ylim(-.001, .002)
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_sunscreen_160mg_2hands_diwater'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        ax_out.hlines(0, 0, 5, linewidth=0.5, color='black') 
        ax_out.set_xlim(0,3)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        
        ax_out = layout.axes[('neutrogena_work', 'neutrogena_work_mono15min')]
        ax_out.set_ylim(-.001, .002)
        path = mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_mono_sunscreen_160mg_2hands_15min'
        scatter_plot_summary_data_out(path, oddonly=False, evenonly=False, baselinesubtract=baselinesubtract, in_or_out='out', ax=ax_out, show_regression=False)
        ax_out.hlines(0, 0, 5, linewidth=0.5, color='black') 
        ax_out.set_xlim(0,3)
        figurefirst.mpl_functions.adjust_spines(ax_out, [])
        
        layout.append_figure_to_layer(layout.figures['neutrogena_work'], 'neutrogena_work', cleartarget=True)
        
        ax = layout.axes[('neutrogena_concentration', 'neutrogena_concentration')]
        plot_sunscreen_mortality(ax)
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(0,7)
        figurefirst.mpl_functions.adjust_spines(ax, [])
        layout.append_figure_to_layer(layout.figures['neutrogena_concentration'], 'neutrogena_concentration', cleartarget=True)
        
        ax = layout.axes[('sunscreen_comparison', 'sunscreen_comparison')]
        plot_sunscreen_type_mortality(ax)
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(0,8)
        figurefirst.mpl_functions.adjust_spines(ax, [])
        layout.append_figure_to_layer(layout.figures['sunscreen_comparison'], 'sunscreen_comparison', cleartarget=True)
        
        ax = layout.axes[('dimethicone', 'dimethicone')]
        plot_dimethicone_mortality(ax)
        ax.set_ylim(-0.1, 1.1)
        ax.set_xlim(0,4)
        figurefirst.mpl_functions.adjust_spines(ax, [])
        layout.append_figure_to_layer(layout.figures['dimethicone'], 'dimethicone', cleartarget=True)
        
        
        if redo_wordcloud:
            import sunscreen_word_cloud
            sunscreen_word_cloud.plot_wordcloud(layout.axes[('wordcloud', 'wordcloud')])
            layout.append_figure_to_layer(layout.figures['wordcloud'], 'wordcloud', cleartarget=True)
    
    layout.write_svg(mono_paper_locations.figure_template_locations.figure5 )
    
    
    
def make_calibration_check_figure():

    layout = get_paper_layout('figure_calibration_check.svg')
    
    axes = {}
    tags = ['before', 'in', 'underwater', 'out', 'after']
    for tag in tags:
        axes[tag] = layout.axes[('field_exps', tag)]
    
    summary_data_out, experiments_to_solutions, experiments_to_colors = make_figure(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_diwater_monowater_carb', axes=axes, solution_to_plot='check_2016', make_ax_spines_clean=False )
    
    ax_out = layout.axes[('field_exps', 'work_out')]
    scatter_plot_summary_data_out(mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_diwater_monowater_carb', oddonly=False, evenonly=False, baselinesubtract=False, in_or_out='out', ax=ax_out, show_regression=False)
    
    layout.append_figure_to_layer(layout.figures['field_exps'], 'field_exps', cleartarget=True)
    layout.write_svg('figure_calibration_check.svg' )
    
    # force trace plot is in mN
    # dip_speed*1e3 is in mm / sec
    
    
    
    
    
    
    

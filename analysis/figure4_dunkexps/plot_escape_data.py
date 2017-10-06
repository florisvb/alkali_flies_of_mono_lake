import csv
import pandas
import matplotlib.pyplot as plt
import numpy as np
import fly_plot_lib.plot as fpl
import fly_plot_lib.text as flytext
import figurefirst
import mono_paper_locations
import scipy.stats

colors = {'Mono Water': 'blue',
           'Na2CO3': 'red',
           'NaCl': 'green',
           'Na2SO4': 'orange',
           'K3PO4': 'cyan',
           'Na2SO3': 'brown',
           'NaOH': (0.001, 0.001, 0.001),
           'NaOH Na2SO4': 'brown',
           'Na2SiO3': 'orange',
           'K2CO3': 'magenta'}

hydration_sphere_delta_r = {'Na': 0.116,
                            'K': 0.074,
                            'CO3': 0.076,
                            'SO4': 0.043,
                            'PO4': 0.054,
                            'OH': 0.079,
                            'Cl': 0.043,
                            'SO3': 0.059,
                            'H2PO4': 0.033,
                            'NO3': 0.044,
                            'I': 0.026,
                            'F': 0.079,
                            'Br': 0.035,
                            'ClO4': 0.019}

ion_radius = {  'Na': 0.102,
                'K': 0.138,
                'CO3': 0.178,
                'SO4': 0.23,
                'PO4': 0.238,
                'OH': 0.133,
                'Cl': 0.181,
                'SO3': 0.2,
                'H2PO4': 0.2,
                'NO3': 0.179,
                'I': 0.22,
                'F': 0.133,
                'Br': 0.196,
                'ClO4': 0.25}

ion_charge = {  'Na': 1,
                'K': 1,
                'CO3': -2,
                'SO4': -2,
                'PO4': -3,
                'OH': -1,
                'Cl': -1,
                'SO3': -2,
                'H2PO4': -1,
                'NO3': -1,
                'I': -1,
                'F': -1,
                'Br': -1,
                'ClO4': -1}

ions_in_molecule = {    'NaCl': ['Na', 'Cl'],
                        'Na2CO3': ['Na', 'CO3'],
                        'Na2SO4': ['Na', 'SO4'],
                        'K2CO3': ['K', 'CO3'],
                        'K3PO4': ['K', 'PO4'],
                        'NaOH': ['Na', 'OH'],
                        'Na2SO3': ['Na', 'SO3']}

def get_delta_hydration_sphere_radius(molecule):
    ions = ions_in_molecule[molecule]
    return hydration_sphere_delta_r[ions[0]] - hydration_sphere_delta_r[ions[1]]

def get_charge_density_surface_area(ion):
    return ion_charge[ion] / (ion_radius[ion])

def get_volume_of_hydration_sphere(ion):
    total_r = ion_radius[ion] + hydration_sphere_delta_r[ion]
    total_vol = 4/3.*np.pi*total_r**3
    ion_vol = 4/3.*np.pi*ion_radius[ion]**3
    return total_vol - ion_vol

def get_sum_of_charge_density_surface_area(molecule):
    ions = ions_in_molecule[molecule]
    print molecule, ions
    print ions[0], ions[1]
    return get_charge_density_surface_area(ions[0]) + get_charge_density_surface_area(ions[1])

def plot_hofmeister_series():
    layout = get_paper_layout('hofmeister')
    ax = layout.axes[('hofmeister', 'hofmeister')]

    ions = ['CO3', 'SO4', 'PO4', 'OH', 'F', 'Cl', 'I', 'ClO4']
    ions = ions[::-1]

    Q_SA = []
    hyd_vol = []
    x = []
    for n, ion in enumerate(ions):
        x.append(n)
        Q_SA.append(np.abs(get_charge_density_surface_area(ion)))
        hyd_vol.append(get_volume_of_hydration_sphere(ion))

    ax.scatter(x, Q_SA, s=35, c='black')

    print '.....'
    print x, len(x)
    print Q_SA, len(Q_SA)
    print hyd_vol, len(hyd_vol)

    results = scipy.stats.linregress( np.array(x), np.array(Q_SA))
    print results
    slope = results[0]
    intercept = results[1]
    rval = results[2]
    pval = results[3]
    xx = np.linspace(-1,len(ions))
    yy = slope*xx + intercept
    ax.plot(xx,yy,zorder=-100,color='black')

    #x = np.linspace(0,30)
    #y = slope*x + intercept

    #ax.plot(x,y,zorder=-100,color='black')

    print ax.get_xlim()

    ax.set_xlim(-0.5,7.5)
    ax.set_ylim(2,14)
    
    figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, yticks=[2,14], xticks=x, linewidth=0.5)

    ax.set_xticklabels([])

    flytext.set_fontsize(ax.figure, 6)

    layout.append_figure_to_layer(layout.figures['hofmeister'], 'hofmeister', cleartarget=True)
    layout.write_svg(mono_paper_locations.figure_template_locations.figure_hofmeister )



def plot_correlation():
    compounds = ['NaCl', 'Na2CO3', 'NaOH', 'Na2SO4',  'K3PO4', 'K2CO3'] 
    dilutions = [16,8,4,2,1]
    data = extract_data_from_excel()

    layout = get_paper_layout('hofmeister')
    ax = layout.axes[('correlation', 'correlation')]
    #fig = plt.figure()
    #ax = fig.add_subplot(111)


    nstuck = []
    DrH_x_QSA = []
    c = []
    for compound in compounds:
        d = []
        for i, dilution in enumerate(dilutions):
            n_flies_stuck, N = get_nflies_trapped_for_compound_and_dilution(data, compound, dilution)
            print n_flies_stuck
            print N

            d.append( np.sum(n_flies_stuck)/float(np.sum(N)) )
        nstuck.append(np.mean(d))
        ions = ions_in_molecule[compound]

        if get_delta_hydration_sphere_radius(compound) > 0:
            DrH_x_QSA.append( np.abs(get_delta_hydration_sphere_radius(compound))*np.abs(get_charge_density_surface_area(ions[1]))**4 )
        else:
            DrH_x_QSA.append( np.abs(get_delta_hydration_sphere_radius(compound))*np.abs(get_charge_density_surface_area(ions[0]))**4 )

        if type(colors[compound]) == str:
            c.append(colors[compound])
        else:
            c.append('black')

    print nstuck
    print DrH_x_QSA
    print c

    ax.scatter( np.array(DrH_x_QSA), np.array(nstuck), c=np.array(c), s=35)

    results = scipy.stats.linregress(np.array(DrH_x_QSA), np.array(nstuck))
    print results
    slope = results[0]
    intercept = results[1]
    rval = results[2]
    pval = results[3]

    x = np.linspace(-100,700)
    y = slope*x + intercept

    ax.plot(x,y,zorder=-100,color='black')

    print ax.get_xlim()

    ax.set_xlim(-100,700)
    ax.set_ylim(0,1)
    
    figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], spine_locations={'left': 5, 'bottom': 5}, yticks=[0,1], xticks=[-100,700], linewidth=0.5)
    ax.set_xticklabels([-1,7])
    

    flytext.set_fontsize(ax.figure, 6)

    layout.append_figure_to_layer(layout.figures['correlation'], 'correlation', cleartarget=True)
    layout.write_svg(mono_paper_locations.figure_template_locations.figure_hofmeister )

def extract_data_from_excel(filename='virilis_escape_data.csv'):
    data = pandas.read_csv(filename)
    data = data.dropna()
    return data

def get_nflies_trapped_for_compound_and_dilution(data, compound, dilution, concentration='default'):
    if concentration == 'default':
        if compound == 'Mono Water':
            concentration = 1
        else:
            concentration = 0.5
    query = 'Compound == "' + compound + '" and Dilution == ' + str(dilution) + ' and Concentration == ' + str(concentration)
    return data.query(query).Captured.values, data.query(query).N.values

def plot_results_for_dilution(data, ax, dilution):

    compounds = ['NaCl', 'Mono Water', 'Na2CO3', 'NaOH', 'Na2SO4',  'Na2SO3', 'NaOH Na2SO4', 'K3PO4', 'Na2SiO3', 'K2CO3'] 
        
    for i, compound in enumerate(compounds):
        n_flies_stuck, N = get_nflies_trapped_for_compound_and_dilution(data, compound, dilution)

        a = None
        for j, n in enumerate(n_flies_stuck):
            tmp = np.array([1]*n + [0]*(N[j]-n)) # 20 flies per experiment default, omit water skiing flies
            if a is None:
                a = tmp     
            else:
                a = np.hstack((a, tmp))

        print compound, n_flies_stuck
        fpl.scatter_box(ax, i+1, a, xwidth=0.4, ywidth=0.1, color=colors[compound], edgecolor='none', flipxy=False, shading='95conf', markersize=1, linewidth=3, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)

    ax.set_ylim(-.02, 1.02)
    ax.set_yticks([0,1])

def plot_results_for_compound(data, ax, compound):
    dilutions = [16,8,4,2,1]
        
    for i, dilution in enumerate(dilutions):
        n_flies_stuck, N = get_nflies_trapped_for_compound_and_dilution(data, compound, dilution)

        a = None
        for j, n in enumerate(n_flies_stuck):
            tmp = np.array([1]*n + [0]*(N[j]-n)) # 20 flies per experiment default, omit water skiing flies
            if a is None:
                a = tmp     
            else:
                a = np.hstack((a, tmp))

        print compound, n_flies_stuck
        fpl.scatter_box(ax, i+1, a, xwidth=0.4, ywidth=0.1, color=colors[compound], edgecolor='none', flipxy=False, shading='95conf', markersize=1, linewidth=3, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)

    ax.set_ylim(-.02, 1.02)
    ax.set_yticks([0,1])

def plot_all_results_by_dilution(axes=None):
    data = extract_data_from_excel()

    dilutions = [16,8,4,2,1]
    if axes is None:
        axes = []
        fig = plt.figure()
        for i, d in enumerate(dilutions):
            ax = fig.add_subplot(1,len(dilutions),i+1)
            axes.append(ax)

    for i, d in enumerate(dilutions):
        plot_results_for_dilution(data, axes[i], d)

def get_paper_layout(fig='escape'):
    if fig == 'escape':
        svg = mono_paper_locations.figure_template_locations.figure4
    elif fig == 'hofmeister':
        svg = mono_paper_locations.figure_template_locations.figure_hofmeister
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

def plot_all_results_by_compound(axes=None):
    data = extract_data_from_excel()

    compounds = ['Mono Water', 'NaCl', 'Na2CO3', 'NaOH', 'Na2SO4', 'K3PO4', 'K2CO3']  
    if axes is None:
        if 0:
            axes = []
            fig = plt.figure()
            for i, compound in enumerate(compounds):
                ax = fig.add_subplot(1,len(compounds),i+1)
                axes.append(ax)
        else:
            layout = get_paper_layout()
            axes = [layout.axes[('virilis_escape', compound)] for compound in compounds]

    for i, compound in enumerate(compounds):
        plot_results_for_compound(data, axes[i], compound)

        if compound == 'Mono Water':
            figurefirst.mpl_functions.adjust_spines(axes[i], ['left'], spine_locations={'left': 5, 'bottom': 5}, linewidth=0.5)
        else:
            figurefirst.mpl_functions.adjust_spines(axes[i], [])


    flytext.set_fontsize(axes[i].figure, 6)

    layout.append_figure_to_layer(layout.figures['virilis_escape'], 'virilis_escape', cleartarget=True)
    layout.write_svg(mono_paper_locations.figure_template_locations.figure4 )



def plot_comparison_at_dilution(dilution=2):
    data = extract_data_from_excel()
    compounds = ['Mono Water', 'NaCl', 'Na2SO4', 'K3PO4', 'NaOH', 'Na2CO3', 'K2CO3']

    layout = get_paper_layout()
    ax = layout.axes[('comparison', 'comparison')]

    for i, compound in enumerate(compounds):
        n_flies_stuck, N = get_nflies_trapped_for_compound_and_dilution(data, compound, dilution)

        a = None
        for j, n in enumerate(n_flies_stuck):
            tmp = np.array([1]*n + [0]*(N[j]-n)) # 20 flies per experiment default, omit water skiing flies
            if a is None:
                a = tmp     
            else:
                a = np.hstack((a, tmp))

        print compound, n_flies_stuck
        fpl.scatter_box(ax, i+1, a, xwidth=0.4, ywidth=0.1, color=colors[compound], edgecolor='none', flipxy=False, shading='95conf', markersize=1, linewidth=3, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)

    ax.set_ylim(-.02, 1.02)
    ax.set_yticks([0,1])
    figurefirst.mpl_functions.adjust_spines(ax, [])

    flytext.set_fontsize(ax.figure, 6)

    layout.append_figure_to_layer(layout.figures['comparison'], 'comparison', cleartarget=True)
    layout.write_svg(mono_paper_locations.figure_template_locations.figure4 )




if __name__ == '__main__':

    plot_all_results_by_compound()
    plot_comparison_at_dilution()
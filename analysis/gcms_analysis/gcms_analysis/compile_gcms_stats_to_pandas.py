import pandas
import read_raw_gcms_data_to_pandas
import numpy as np
import sys

from mono_analysis import plot_raw_traces
import itertools
import matplotlib.pyplot as plt
import copy
import statsmodels.formula.api as sm
import fly_plot_lib.plot as fpl
import fly_plot_lib.colormaps
import figurefirst
import scipy.cluster
import fly_plot_lib.text as flytext

import mono_paper_locations

SPECIES = ['bluekelp', 'blackkelp', 'ephydra', 'santaana', 'oilfly', 'melanogaster', 'virilis']
ARBITRARY_HAIRINESS = {'bluekelp': 8,
                       'blackkelp': 5,
                       'ephydra': 10,
                       'santaana': 4,
                       'oilfly': 1,
                       'virilis': 6,
                       'melanogaster': 6}
                       
HAIR_TYPE_RATIO = {    'bluekelp': 8,                # large number means more small hairs per large hair | on lowerleg only
                       'blackkelp': 2,
                       'ephydra': 10,
                       'santaana': 1,
                       'oilfly': 3,
                       'virilis': 3,
                       'melanogaster': 2}

# calculate hairiness with: '/home/caveman/Documents/Projects/Alkali_Fly/SEM/20161118/Hairyness.ipynb'
hairiness_filename = mono_paper_locations.data_locations.hairiness_pickle
hairiness_pd = pandas.read_pickle(hairiness_filename)
def get_hairyness_for_species_and_bodypart_from_saved_data(pd, species, bodypart):
    if bodypart == 'all':
        return (pd[pd.species==species].mean(axis=1).values[0])
    elif bodypart == 'arbitrary':
        return ARBITRARY_HAIRINESS[species]     
    elif bodypart == 'hair_type_ratio':
        return HAIR_TYPE_RATIO[species]        
    else:
        return (pd[pd.species==species][bodypart].values[0])


hairyness = {}
for species in SPECIES:
    h = get_hairyness_for_species_and_bodypart_from_saved_data(hairiness_pd, species, 'all')
    hairyness[species] = h
    
pulvilli = { 'ephydra': 0,
             'bluekelp': 10,
             'virilis': 5,
             'santaana': 8,
             'melanogaster': 3,
             'blackkelp': 10,
             'oilfly': 10,
             }

scaling = {'ephydra': 6,      # approximate millimeters, length
           'bluekelp': 6.2,
           'blackkelp': 7,
           'oilfly': 2.3,
           'virilis': 3.5,
           'melanogaster': 2.5,
           'santaana': 1.5,
           }
           
dip_data_paths = {  'blackkelp': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/blackkelp_2016_diwater_monowater_carb',
                    'melanogaster': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/melanogaster_2016_diwater_monowater_carb',
                    'ephydra': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/ephydra_2016_diwater_monowater_carb',
                    'virilis': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/virilis_2016_diwater_monowater_carb',
                    'bluekelp': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/bluekelp_2016_diwater_monowater_carb',
                    'oilfly': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/oilflies_2016_diwater_monowater_carb',
                    'santaana': mono_paper_locations.data_locations.mono_lake_lab_data_basepath+'/santa_ana_2016_diwater_monowater_carb',
                    }

GCMS_path = mono_paper_locations.data_locations.mono_lake_gcms_raw_data


def rearrange_pd_to_single_output(pd):
    pd_diwater = copy.deepcopy(pd)
    pd_diwater['carbonate_concentration'] = 0.
    pd_diwater['water_type'] = 'diwater'
    pd = remove_size_correlation(pd_diwater, 'diwater')
    pd_diwater['work_out'] = pd['diwater_minus_size'].values
    
    pd_monowater = copy.deepcopy(pd)
    pd_monowater['carbonate_concentration'] = 2.
    pd_monowater['water_type'] = 'monowater'
    pd = remove_size_correlation(pd_monowater, 'monowater')
    pd_monowater['work_out'] = pd['monowater_minus_size'].values
    
    pd_carb = copy.deepcopy(pd)
    pd_carb['carbonate_concentration'] = 4.
    pd_carb['water_type'] = 'carb'
    pd = remove_size_correlation(pd_carb, 'carb')
    pd_carb['work_out'] = pd['carb_minus_size'].values
    
    pd_new = pandas.concat([pd_carb, pd_monowater, pd_diwater])
    
    return pd_new
    
def make_solution_size_summary_figure():
    layout = get_paper_layout()
    pd = compile_all_data(return_mean=False) 
    ''' False means each individual animal is brought in for the statistics, 
    instead of the species average. Note, however, that the same size is used for all individuals. 
    This could potentially be solved using the stored video data of the dipping experiments. 
    But it is not really necessary. Result is just not quite as strong as it could be. '''
    plot_size_correlations(pd, layout.axes[('solution_summary','solution_summary')])
    
    ax = layout.axes[('solution_summary','solution_summary')]
    #figurefirst.mpl_functions.adjust_spines(ax, ['left', 'bottom'], xticks=[0,1], yticks=[0,.001])
    figurefirst.mpl_functions.adjust_spines(ax, ['bottom'], xticks=[1,7.5], spine_locations={'left': 1, 'bottom': 1})
    ax.set_xticklabels(['1','7.5'])
    flytext.set_fontsize(ax.figure, 8)
    
    write_figure_to_svg(layout, 'solution_summary')

def print_rsq_for_solution_all_species(solution='diwater', remove_ephydra=False):
    pd = compile_all_data(return_mean=False)
    if remove_ephydra:
        pd = pd[pd.species!='ephydra']
    formula = solution + ' ~ size'
    result = sm.ols(formula=formula, data=pd).fit()
    print 'solution: ', solution
    print '     pval: ', result.pvalues['size']
    print '     rsq_adj: ', result.rsquared_adj
    print '     rsq: ', result.rsquared


def plot_size_correlations(pd, ax):
    # for best results use:
    # pd = compile_stats_to_pandas.compile_all_data(return_mean=False)
    pd_tmp = copy.deepcopy(pd)
    pd_tmp = pd_tmp[pd_tmp.species!='ephydra']

    def add_regression_line(pd, ax, solution):
        colors = {'diwater': 'green', 'monowater': 'blue', 'carb': 'red'}
        formula = solution + ' ~ size'
        result = sm.ols(formula=formula, data=pd).fit()
        xs = np.linspace( np.min(pd['size'])-0.5, np.max(pd['size'] )+0.5, 100)
        
        slope_lo = result.conf_int()[0]['size']
        slope_hi = result.conf_int()[1]['size']
        
        y1 = xs*slope_lo + result.params['Intercept']
        y2 = xs*slope_hi + result.params['Intercept']
        
        y = np.vstack([y1,y2])
        
        y_max = np.max(y, axis=0)
        y_min = np.min(y, axis=0)
        ax.fill_between(xs, y_max, y_min, edgecolor='none', facecolor=colors[solution], alpha=0.3)
        
        ax.plot(xs, xs*result.params['size']+result.params['Intercept'], colors[solution], zorder=-200)
        print 'solution: ', solution
        print '     pval: ', result.pvalues['size']
        print '     rsq_adj: ', result.rsquared_adj
        print '     conf_int: ', slope_lo, slope_hi


    for species in SPECIES:
        if species == 'ephydra':
            pd_tmp_species = pd[pd['species']==species]
            ax.plot(np.mean(pd_tmp_species['size'])-0.1, np.mean(pd_tmp_species['diwater']), '*', color='green', markeredgecolor='none', markersize=6, zorder=100)
            ax.plot(np.mean(pd_tmp_species['size'])+0.1, np.mean(pd_tmp_species['monowater']), '*', color='blue', markeredgecolor='none', markersize=6, zorder=100)
            ax.plot(np.mean(pd_tmp_species['size']), np.mean(pd_tmp_species['carb']), '*', color='red', markeredgecolor='none', markersize=6, zorder=100)
        else:
            pd_tmp_species = pd_tmp[pd_tmp['species']==species]
            ax.plot(np.mean(pd_tmp_species['size']), np.mean(pd_tmp_species['diwater']), '.', color='green', markeredgecolor='none', markersize=4, zorder=-100)
            ax.plot(np.mean(pd_tmp_species['size']), np.mean(pd_tmp_species['monowater']), '.', color='blue', markeredgecolor='none', markersize=4, zorder=-100)
            ax.plot(np.mean(pd_tmp_species['size']), np.mean(pd_tmp_species['carb']), '.', color='red', markeredgecolor='none', markersize=4, zorder=-100)
    
    for solution in ['diwater', 'monowater', 'carb']:
        add_regression_line(pd_tmp, ax, solution)
    
    if 0:
        pd_ephydra = pd[pd.species=='ephydra']
    
        fpl.scatter_box(ax, np.mean(pd_ephydra['size'].values)-0.2, pd_ephydra['diwater'].values, xwidth=0.1, ywidth=0.1, color='white', edgecolor='white', flipxy=False, shading='95conf', markersize=.5, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)    
        fpl.scatter_box(ax, np.mean(pd_ephydra['size'].values)-0.2, pd_ephydra['diwater'].values, xwidth=0.1, ywidth=0.1, color='green', edgecolor='green', flipxy=False, shading='95conf', markersize=.5, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)
        
        
        fpl.scatter_box(ax, np.mean(pd_ephydra['size'].values), pd_ephydra['monowater'].values, xwidth=0.1, ywidth=0.1, color='white', edgecolor='white', flipxy=False, shading='95conf', markersize=.5, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)
        fpl.scatter_box(ax, np.mean(pd_ephydra['size'].values), pd_ephydra['monowater'].values, xwidth=0.1, ywidth=0.1, color='blue', edgecolor='blue', flipxy=False, shading='95conf', markersize=.5, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)
        
        fpl.scatter_box(ax, np.mean(pd_ephydra['size'].values)+0.2, pd_ephydra['carb'].values, xwidth=0.1, ywidth=0.1, color='white', edgecolor='white', flipxy=False, shading='95conf', markersize=.5, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)
        fpl.scatter_box(ax, np.mean(pd_ephydra['size'].values)+0.2, pd_ephydra['carb'].values, xwidth=0.1, ywidth=0.1, color='red', edgecolor='red', flipxy=False, shading='95conf', markersize=.5, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=True)
        
    ax.hlines(0, 0, 8)
        
    ax.set_xlim([1,7.5])
    
def plot_correlation_with_size_removed(pd, output, inputs, ax=None):
    pd_new = rearrange_pd_to_single_output(pd)
    pd_new = pd_new[pd_new.water_type == output]
    pd_tmp = get_pd_tmp(pd_new, remove_species=[]) # whitened
    print pd_tmp.water_type
    
    ## stats
    
    maxnumcombinations = 3
    formulas = []
    output = 'work_out'
    for r in range(1,maxnumcombinations+1):
        combs = itertools.combinations(inputs, r)
        for comb in combs:
            f = output + ' ~ ' + ' + '.join(comb)
            formulas.append(f)
    
    formula_results = {}
            
    
    for formula in formulas:    
        # with ephydra
        result = sm.ols(formula=formula, data=pd_tmp).fit()
        formula_results[formula] = result.rsquared  
        
    
    
    formula_results_name = {1: [], 2: [], 3: []}
    formula_results_rsq = {1: [], 2: [], 3: []}
    for formula, rsq in formula_results.items():
        n = len(formula.split(' ~ ')[-1].split(' + '))
        formula_results_name[n].append(formula)
        formula_results_rsq[n].append(rsq)

    for n in range(1,4):
        order = np.argsort(formula_results_rsq[n])[::-1]
        formula_results_name[n] = [formula_results_name[n][o] for o in order]
        formula_results_rsq[n] = [formula_results_rsq[n][o] for o in order]
        
    if ax is None:
        layout = get_paper_layout()
        ax = layout.axes[('correlations', 'correlations')]
        write_svg = True
    else:
        write_svg = False
        
    print 'correlation r squareds'
    for n in range(1,4):
        for i, rsq in enumerate(formula_results_rsq[n]):
            if i==0:
                print n, rsq
                ax.plot(n, rsq, 'o', color='black', markersize=2, markeredgecolor='none')
            else:
                ax.plot(n, rsq, 'o', color='gray', markersize=2, markeredgecolor='none')
    
    # print best formulas
    for n in range(1,4):
        print formula_results_name[n][0], formula_results_rsq[n][0]
    
    ax.set_ylim(0,1)
    ax.set_xlim(0.9, 3.1)
    figurefirst.mpl_functions.adjust_spines(ax, [])
    
    if write_svg:
        layout.append_figure_to_layer(layout.figures['correlations'], 'correlations', cleartarget=True)
        layout.write_svg(mono_paper_locations.figure_template_locations.figure3)
    
    return formula_results

def get_paper_layout(subfig='species'):
    if subfig == 'chc':
        svg = mono_paper_locations.figure_template_locations.figure3_chc
    else:
        svg = mono_paper_locations.figure_template_locations.figure3
    layout = figurefirst.svg_to_axes.FigureLayout(svg)
    layout.make_mplfigures()
    return layout

def write_figure_to_svg(layout, figure_name, subfig='species'):
    layout.append_figure_to_layer(layout.figures[figure_name], figure_name, cleartarget=True)
    if subfig == 'chc':
        layout.write_svg(mono_paper_locations.figure_template_locations.figure3_chc)
    else:
        layout.write_svg(mono_paper_locations.figure_template_locations.figure3)

def simplify_gcms(gcms_data):
    
    Codd = 0
    Ceven = 0
    Codd_Me = 0
    Ceven_Me = 0
    Codd_alkene_diene = 0
    Ceven_alkene_diene = 0
    Codd_branched = 0
    Codd_straight = 0
    Ceven_branched = 0
    Ceven_straight = 0
    
    _Codd = ['C' + str(n) for n in [21,23,25,27,29,31,33]]
    _Ceven = ['C' + str(n) for n in [22,24,26,28,30,32]]
    
    for row in gcms_data.iterrows():
        row = row[1]
        if row.compound_name in _Codd:
            Codd += row.perc_total
        elif row.compound_name in _Ceven:
            Ceven += row.perc_total
        elif 'Me' in row.compound_name:
            Cidx = row.compound_name.index('C')
            C = row.compound_name[Cidx:Cidx+3]
            if C in _Codd:
                Codd_Me += row.perc_total
            elif C in _Ceven:
                Ceven_Me += row.perc_total
        elif 'alkene' in row.compound_name or 'diene' in row.compound_name:
            Cidx = row.compound_name.index('C')
            C = row.compound_name[Cidx:Cidx+3]
            if C in _Codd:
                Codd_alkene_diene += row.perc_total
            elif C in _Ceven:
                Ceven_alkene_diene += row.perc_total
        
        if 'Me' in row.compound_name or 'alkene' in row.compound_name or 'diene' in row.compound_name:
            Cidx = row.compound_name.index('C')
            C = row.compound_name[Cidx:Cidx+3]
            if C in _Codd:
                Codd_branched += row.perc_total
            elif C in _Ceven:
                Ceven_branched += row.perc_total
        #elif len(row.compound_name) == 3:
        #    if row.compound_name in _Codd or row.compound_name in _Ceven: 
        #        Cstraight += row.perc_total
    
    ret_time_avg_w = np.sum(gcms_data.ret_time_min*gcms_data.perc_total/(np.sum(gcms_data.perc_total)))
    corr_area = np.sum(gcms_data.corr_area)
    
    return {'RetTimeAvg': ret_time_avg_w,
            'corr_area': corr_area,
            'Codd': Codd,
            'Ceven': Ceven,
            'Codd_Me': Codd_Me,
            'Ceven_Me': Ceven_Me,
            'Codd_alkene_diene': Codd_alkene_diene,
            'Ceven_alkene_diene': Ceven_alkene_diene,
            #'Codd_straight': Codd_straight,
            'Codd_branched': Codd_branched,
            #'Ceven_straight': Ceven_straight,
            'Ceven_branched': Ceven_branched,
            }

def compile_data_for_species(species, return_mean=True):
    data_series = {}
    r, di_mono_carb = plot_raw_traces.scatter_plot_summary_data_out(dip_data_paths[species])
    
    if return_mean:
        data_series['diwater'] = np.mean(di_mono_carb['diwater_2016'])
        data_series['monowater'] = np.mean(di_mono_carb['monowater_2016'])
        data_series['carb'] = np.mean(di_mono_carb['carb'])
        data_series['species'] = species
        data_series['size'] = scaling[species]
        data_series['hairyness'] = hairyness[species]
        data_series['pulvilli'] = pulvilli[species]
        gcms_data = read_raw_gcms_data_to_pandas.parse_rows(GCMS_path, species, start=17)
        gcms_simple = simplify_gcms(gcms_data)
        
        data_series.update(gcms_simple)
        
        plt.close('all')
        return data_series
    
    else:
        n = len(di_mono_carb['carb'])
        data_series['diwater'] = di_mono_carb['diwater_2016']
        data_series['monowater'] = di_mono_carb['monowater_2016']
        data_series['carb'] = di_mono_carb['carb']
        data_series['species'] = [species]*n
        data_series['size'] = [scaling[species]]*n
        data_series['hairyness'] = [hairyness[species]]*n
        data_series['pulvilli'] = [pulvilli[species]]*n
        gcms_data = read_raw_gcms_data_to_pandas.parse_rows(GCMS_path, species, start=17)
        gcms_simple = simplify_gcms(gcms_data)
        
        data_series.update(gcms_simple)
        
        plt.close('all')
        return data_series
    
    
def compile_all_data(return_mean=True):

    if not return_mean:
        pds = []
        for species in SPECIES:
            species_data = compile_data_for_species(species, return_mean=return_mean)
            pd_species = pandas.DataFrame(species_data)
            pds.append(pd_species)
        pd = pandas.concat(pds)
        return pd
        
    else:
        pd = pandas.DataFrame()
        for species in SPECIES:
            species_data = compile_data_for_species(species, return_mean=return_mean)
            pd = pd.append(species_data, ignore_index=True)
        return pd
                
def analyze_pandas_dataframe(pd, formulas):
    
    
    pd_tmp = copy.copy(pd)
    pd_tmp = pd_tmp.drop('species', 1)

    # normalize / whiten            
    pd_tmp = pd_tmp - pd_tmp.mean() # zero mean
    pd_tmp = scipy.cluster.vq.whiten(pd_tmp) # variance of 1
    
    def get_mse_resid(formula, pd_tmp):
        result = sm.ols(formula=formula, data=pd_tmp).fit()
        return result.mse_resid
    
    results = {}
    for formula in formulas:
        mse_resid = get_mse_resid(formula, pd_tmp)
        results[formula] = mse_resid
        
    return results

def analyze_pandas_dataframe_detailed(pd, formulas):
    
    pd_tmp = copy.copy(pd)
    try:
        pd_tmp = pd_tmp.drop('species', 1)
    except:
        pass
        
    # normalize / whiten            
    pd_tmp = pd_tmp - pd_tmp.mean()
    pd_tmp = scipy.cluster.vq.whiten(pd_tmp)
    
    results = pandas.DataFrame()
    for formula in formulas:
        result = sm.ols(formula=formula, data=pd_tmp).fit()
        series = copy.copy(result.params)
        for key in series.keys():
            if result.pvalues[key] > 0.05:
                series[key] = np.nan
        series['formula'] = formula
        series['resid'] = result.mse_resid
        series['f_pvalue'] = result.f_pvalue
        series['df_model'] = result.df_model
        series['n_inputs'] = len(formula.split('~')[1].split('+'))
        results = results.append(series, ignore_index=True)
        
    return results

def get_model(pd, formula):
    pd_tmp = copy.copy(pd)
    try:
        pd_tmp = pd_tmp.drop('species', 1)
    except:
        pass
        
    # normalize / whiten            
    pd_tmp = pd_tmp - pd_tmp.mean()
    pd_tmp = scipy.cluster.vq.whiten(pd_tmp)
    
    result = sm.ols(formula=formula, data=pd_tmp).fit()
        
    return result
    
def get_OLS_formula_results(pd, output='carb', maxnumcombinations=5):
    res = []
    #inputs = ['size', 'hairyness', 'Ceven', 'Codd', 'Ceven_Me', 'Codd_Me', 'Ceven_alkene_diene', 'Codd_alkene_diene']
    inputs = ['pulvilli', 'size', 'hairyness', 'Codd', 'Ceven_Me', 'Codd_Me', 'Codd_alkene_diene']
    #inputs = ['RetTimeAvg', 'hairyness', 'size', 'pulvilli']
    #inputs = ['size', 'hairyness', 'Ceven', 'Codd', 'Cstraight', 'Cbranched']
    #inputs = ['size', 'hairyness', 'Cstraight', 'Cbranched']
    #inputs = ['size', 'hairyness', 'Ceven', 'Codd']
    formulas = []
    formula_combinations = {}
    
    for r in range(1,maxnumcombinations):
        combs = itertools.combinations(inputs, r)
        for comb in combs:
            f = output + ' ~ ' + ' + '.join(comb)
            formulas.append(f)
            formula_combinations[f] = comb
            
    results = analyze_pandas_dataframe_detailed(pd, formulas)
    
    return results
    

def get_pd_tmp(pd, remove_species=[]):

    pd_tmp = copy.copy(pd)
    for species in remove_species:
        pd_tmp = pd_tmp[pd_tmp.species!=species]
    

    # normalize / whiten       
    for key in pd_tmp.keys():     
        if key == 'species' or key == 'water_type' or key == 'carbonate_concentration':
            continue
        else:
            print key
            pd_tmp[key] = pd_tmp[key] - pd_tmp[key].mean() # zero mean
            pd_tmp[key] = scipy.cluster.vq.whiten(pd_tmp[key]) # variance of 1
    
    return pd_tmp
    
def remove_size_correlation(pd, solution):
    # leaving ephydra out
    pd_tmp = copy.deepcopy(pd)
    pd_noephydra = pd_tmp[pd_tmp.species!='ephydra']
    formula = solution + ' ~ size'
    result = sm.ols(formula=formula, data=pd_noephydra).fit()
    pd_tmp[solution+'_minus_size'] = pd_tmp[solution] - pd_tmp['size']*result.params['size']
    return pd_tmp
    
    
def plot_OLS_results(pd, solution, results=None, ax=None, skip_species=['ephydra']):

    pd_tmp = copy.copy(pd)
    for species in skip_species:
        pd_tmp = pd_tmp[pd_tmp.species!=species]
    

    if results is None:
        print 'getting results'
        results = get_OLS_formula_results(pd_tmp, solution)
        
    # solution = diwater, monowater, or carb
    layout = get_paper_layout()
    fname = solution + '_OLS'
    if ax is None:
        ax = layout.axes[(fname, fname)]
        svg = True
    else:
        svg = False
    
    # 
    keys_to_plot = ['Codd', 'Codd_Me', 'Codd_alkene_diene', 'Ceven', 'Ceven_Me', 'Ceven_alkene_diene', 'hairyness', 'size', 'pulvilli']
    #keys_to_plot = ['RetTimeAvg', 'hairyness', 'size', 'pulvilli']
    labels = []
    xticks = []
    
    n = -1
    for key in keys_to_plot:
        n += 1
        
        try:
            y_data = results[key].values
        except:
            continue
        indices = ~np.isnan(y_data)
        y_data = y_data[indices]#*(1-results['resid'].values[indices])
        color = results['resid'].values[indices]
        
        if len(y_data) >= 1:
        
            fpl.scatter_box(ax, n, y_data, xwidth=0.4, ywidth=0.1, color=(0.001, 0.001, 0.001), edgecolor='black', flipxy=False, shading='95conf', markersize=10, linewidth=1, use='mean', optimize_scatter_distance=False, optimize_scatter_distance_resolution=20, optimize_scatter_distance_y_scale=1, hide_markers=False, scatter_color=color, scatter_norm_minmax=[0,0.5])
            
            print key
            print y_data
            print n
            print
        
    
        labels.append(key)
        xticks.append(n)
    
    ax.hlines(0, -1, n+1) 
    ax.set_xlim(-0.5,n+0.5)
    ax.set_ylim(-2,2)
    #ax.set_ylim(-10,10)
    
    figurefirst.mpl_functions.adjust_spines(ax, [])
    #ax.set_xticks(xticks)
    #ax.set_xticklabels(labels, rotation='vertical')
    
    if svg:
        write_figure_to_svg(layout, fname)

def get_n_best_formulas(results, pd, n=5, df=4, sortby='resid'):
    tmp = results[results.df_model==df]
    sorted_resid = tmp[sortby].sort_values()
    best_formulas = [tmp.loc[i]['formula'] for i in sorted_resid.index[0:n]]
    
    OLS_results = []
    for formula in best_formulas:
        OLS_results.append( get_model(pd, formula) )
    
    print 'Formula, f_pvalue, rsq_adj, mse_resid'
    for i in range(0,5):
        try:
            print best_formulas[i], OLS_results[i].f_pvalue, OLS_results[i].rsquared_adj, OLS_results[i].mse_resid
        except:
            pass
    print
    
    return OLS_results
        
        

def plot_CHC_matrix(pd, layout):
    ax = layout.axes[('correlation_matrix', 'CHC_matrix')]
    
    HC_mixes = ['Codd', 'Codd_Me', 'Codd_alkene_diene', 'Ceven', 'Ceven_Me', 'Ceven_alkene_diene']

    CHC_matrix = np.zeros([len(SPECIES),len(HC_mixes)])

    for r, species in enumerate(SPECIES):
        for c, HC in enumerate(HC_mixes):
            CHC_matrix[r,c] = pd[pd.species==species][HC].values[0] / 100.
    
    ax.imshow(CHC_matrix, cmap=fly_plot_lib.colormaps.viridis, interpolation='nearest', origin='upper', extent=[0,len(HC_mixes),0,len(SPECIES)])
    ax.set_aspect('auto')
    figurefirst.mpl_functions.adjust_spines(ax, [])
    return layout 
    
def plot_hairyness(pd, layout):
    ax = layout.axes[('correlation_matrix', 'hairy_matrix')]
    
    hairyness = np.array([pd[pd.species==species]['hairyness'].values[0] for species in SPECIES])
    hairyness = hairyness.reshape([7,1])
    ax.imshow(hairyness, cmap=fly_plot_lib.colormaps.viridis, interpolation='nearest', origin='upper', extent=[0,1,0,len(SPECIES)])
    ax.set_aspect('auto')
    figurefirst.mpl_functions.adjust_spines(ax, [])
    return layout
    
def plot_AvgRetTime(pd, layout):
    ax = layout.axes[('correlation_matrix', 'AvgRetTime_matrix')]
    
    AvgRetTime = np.array([pd[pd.species==species]['RetTimeAvg'].values[0] for species in SPECIES])
    AvgRetTime = AvgRetTime.reshape([7,1])
    ax.imshow(AvgRetTime, cmap=fly_plot_lib.colormaps.viridis, interpolation='nearest', origin='upper', extent=[0,1,0,len(SPECIES)])
    ax.set_aspect('auto')
    figurefirst.mpl_functions.adjust_spines(ax, [])
    return layout
    
def plot_pulvilli(pd, layout):
    ax = layout.axes[('correlation_matrix', 'pulvilli_matrix')]
    
    pulvilli = np.array([pd[pd.species==species]['pulvilli'].values[0] for species in SPECIES])
    pulvilli = pulvilli.reshape([7,1])
    ax.imshow(pulvilli, cmap=fly_plot_lib.colormaps.viridis, interpolation='nearest', origin='upper', extent=[0,1,0,len(SPECIES)])
    ax.set_aspect('auto')
    figurefirst.mpl_functions.adjust_spines(ax, [])
    return layout
    
def plot_size(pd, layout):
    ax = layout.axes[('correlation_matrix', 'size_matrix')]
    
    size = np.array([pd[pd.species==species]['size'].values[0] for species in SPECIES])
    size = size.reshape([7,1])
    ax.imshow(size, cmap=fly_plot_lib.colormaps.viridis, interpolation='nearest', origin='upper', extent=[0,1,0,len(SPECIES)])
    ax.set_aspect('auto')
    figurefirst.mpl_functions.adjust_spines(ax, [])
    return layout
    
def make_correlation_matrix_figure(pd):
    layout = get_paper_layout(subfig='chc')
    
    layout = plot_CHC_matrix(pd, layout)
    layout = plot_hairyness(pd, layout)
    layout = plot_size(pd, layout)
    layout = plot_pulvilli(pd, layout)
    layout = plot_AvgRetTime(pd, layout)
    
    write_figure_to_svg(layout, 'correlation_matrix', subfig='chc')
    
def make_resid_colorbar():
    layout = get_paper_layout()
    
    fpl.colorbar(ax=layout.axes[('resid_colorbar','resid_colorbar')], ticks=None, ticklabels=None, colormap='jet', aspect='auto', orientation='vertical', filename=None, flipspine=False, show_spine=False)
    
    write_figure_to_svg(layout, 'resid_colorbar')
    
    
    
    
    
    
'''

pd = compile_stats_to_pandas.compile_all_data()
compile_stats_to_pandas.plot_OLS_results(pd, 'diwater')
compile_stats_to_pandas.plot_OLS_results(pd, 'monowater')
compile_stats_to_pandas.plot_OLS_results(pd, 'carb')
compile_stats_to_pandas.make_correlation_matrix_figure(pd)
compile_stats_to_pandas.make_resid_colorbar()

'''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

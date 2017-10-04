import pandas
import csv
import os
import numpy as np


def load_file_rows(filename):
    rows = []
    with open(filename, 'rb') as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            rows.append(row)
    return rows
    
def load_compound_names(filename):
    rows = load_file_rows(filename)
    names = {}
    for row in rows:
        ret_time, name = row
        names[float(ret_time)] = name
    return names

def get_filename(path, contains, does_not_contain=[]):
    cmd = 'ls ' + path
    ls = os.popen(cmd).read()
    all_filelist = ls.split('\n')
    try:
        all_filelist.remove('')
    except:
        pass
    filelist = []
    for i, filename in enumerate(all_filelist):
        if contains in filename:
            fileok = True
            for nc in does_not_contain:
                if nc in filename:
                    fileok = False
            if fileok:
                return os.path.join(path, filename)
                
def get_name_from_ret_time(compound_names, ret_time):
    err = np.min( np.abs( np.array(compound_names.keys()) - ret_time ) )
    if err < 0.02:
        idx = np.argmin( np.abs( np.array(compound_names.keys()) - ret_time ) )
        return compound_names[compound_names.keys()[idx]]
    else:
        return '?'
        
def parse_rows(path, species, start=17):
    quantitative_filename = get_filename(path, species+'.txt')
    rows = load_file_rows(quantitative_filename)

    compound_names_filename = get_filename(path, species+'_compound_names.csv')
    compound_names = load_compound_names(compound_names_filename)

    pd = pandas.DataFrame()
    columns = {'peak_number': 0, 'ret_time_min': 1, 'perc_max': 9, 'perc_total': 10, 'corr_area': 8}
    column_data = []
    column_names = ['peak_number', 'ret_time_min', 'perc_max', 'perc_total', 'corr_area', 'species', 'compound_name']
    
    for row in rows[17:]:
        print row
        print
        data = row[0].split(' ')
        data = [d for d in data if d != '']
        
        tmp = []
        for column_name in column_names:
            if column_name == 'species':
                tmp.append(species)
            elif column_name == 'compound_name':
                ret_time_index = columns['ret_time_min']
                ret_time = float(data[ret_time_index].rstrip('%'))
                compound_name = get_name_from_ret_time(compound_names, ret_time)
                tmp.append(compound_name)
            else:
                index = columns[column_name]
                print column_name, index
                num = float(data[index].rstrip('%'))
                tmp.append(num)        
                
        
        column_data.append(tmp)
    
    pd = pandas.DataFrame(np.array(column_data), columns=column_names)
    pd = pd.convert_objects(convert_numeric=True)
    return pd
    
    
    
    
    
    
    
    
    

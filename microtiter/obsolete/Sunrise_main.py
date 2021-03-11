# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 16:43:47 2018

@author: James
"""
import Concordia.Lab_Automation.microtiter.curve_analysis as curve_analysis
import Concordia.Lab_Automation.microtiter.plate_reader_tools as tools
import numpy as np
from collections import namedtuple

def make_plate(species_path=None,
               sunrise_path=None,
               treatment_path=None,
               strains_path=None,
               plate_name=None):

    paths = {'species path': species_path,
             'sunrise path': sunrise_path,
             'treatment path': treatment_path,
             'strains path': strains_path}

    for path_name, path in paths.items():
        if path:
            continue
        else:
            paths[path_name] = input(
                    'Paste in {}\n'.format(path_name)).strip('"')

    # species_path =
    # 'C:\Users\Owner\Dropbox\BioFoB\acid-tolerant_strains_02.xlsx'

    if plate_name == None:
        plate_name = input('enter name for plate\n')

    plate = curve_analysis.Plate(tools.read_sunrise(
            paths['sunrise path']), plate_name)
    plate.paths = paths

    plate.set_replicates('media', ['map', paths['treatment path']])
    plate.set_strains(paths['strains path'])
    plate.set_species(paths['species path'])

    plate.correct_for_blank()

    return plate


def get_integrals(plate, drop=False, timepoint=24, average=True):

    columns = plate({'cols': lambda x: x.integral[-1].index})[0]
    closest_timepoint = (np.abs(columns - timepoint)).argmin()
    timepoint = columns[closest_timepoint]

    integral = plate({'strain': lambda x: x.strain,
                      'species': lambda x: x.species,
                      'integral': lambda x: x.integral[-1].loc[timepoint],
                      'well': lambda x: x.name})

    if drop:
        integral = integral.drop(drop)
    integral.index.name = 'well'
    integral.set_index(['species', 'strain', 'well'], inplace=True)
    integral = integral.astype(float)
    integral.sort_index(level=['species', 'strain', 'well'], inplace=True)
    integral_groups = integral.groupby(
            level=['species', 'strain'])
    if average:
        average = integral_groups.mean()
        average.columns = ['integral mean']
        std = integral_groups.std()
        std.columns = ['integral stdev']
    
        results = average.join(std)
       
    else:
        results = integral
#    axis = plt(average)

    print('OMMITTED: {}'.format(drop))

    return integral, results

def draw_catplot(dataframe, **kwargs):
    import seaborn as sns
    import datetime
    
    x, y, hue = kwargs['x'], kwargs['y'], kwargs['hue']
    
    num_of_strains = len(dataframe[x].unique())
    #aspect = 10/num_of_strains
    aspect = 100/num_of_strains
    height = num_of_strains / 16
    catplot = sns.catplot(x=x, y=y, hue=hue, data=dataframe,
                          height=height, aspect=aspect)
    #catplot.set_xticklabels(rotation=90)

    time = str(datetime.datetime.now())[:-7].replace(':', '-')
    catplot.savefig(fname='catplot {}.svg'.format(time))
    catplot.savefig(fname='catplot {}.png'.format(time), dpi=200)

def draw_violin(x, y, data, savename):
    import seaborn as sns
    #plot = sns.catplot(data, kind='violin', cut=0)    
    plot = sns.catplot(x = x, y = y, data = data, kind='violin', cut=0, aspect=1.4, bw=0.25, palette='muted')
    plot.savefig(savename, dpi=200)
    
def sunrise_file_finder(root):
    def file_checker(files, check):
         for file in files:
             if '.xlsx' in file:
                 if str(check).lower() in str(file).lower():
                     return file
         else:
             return False
    
    import os
    kinds = ['aa20', 'ph3', 'ypd']
    file = namedtuple('file', ['path', 'parent', 'treatment'])
    files = []
    
    for entry in os.listdir(root):
        path = os.path.join(root,entry)
        if os.path.isdir(path):
            for kind in kinds:
                if file_checker(os.listdir(path), kind):
                    files.append(file(
                        path = os.path.join(root, entry,
                                file_checker(os.listdir(path), kind)),
                        parent = entry,
                        treatment = kind))
                    
    return(files)

def treatment_file_finder(root):
    def file_checker(files, check):
         for file in files:
             if '.xlsx' in file:
                 if str(check).lower() in str(file).lower():
                     return file
         else:
             return False
    
    import os
    kinds = ['Strain map']
    file = namedtuple('file', ['path', 'parent', 'treatment'])
    files = []
    
    for entry in os.listdir(root):
        path = os.path.join(root,entry)
        if os.path.isdir(path):
            for kind in kinds:
                if file_checker(os.listdir(path), kind):
                    files.append(file(
                        path = os.path.join(root, entry,
                                file_checker(os.listdir(path), kind)),
                        parent = entry,
                        treatment = kind))
                    
    return(files)

def congregator(files):
    dfs = []

    for file in files:
        dfs.append(tools.read_sunrise(file.path))
        dfs[-1].parent = file.parent
        dfs[-1].treatment = file.treatment
        
    return dfs
    
    
    
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 16:47:33 2018

@author: James
"""
import plate_reader_tools as tools
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def curve_maker(sunrise_path,
                strain_map_path=r"../microtiter/empty_strain_map.xlsx",
                raw_time=False):
    
    sunrise = tools.read_sunrise(sunrise_path)
    strain_map = tools.read_infinateM200_output(strain_map_path)
    sunrise.index.name = 'Time'
    strain_map = strain_map.dataframe.reset_index()
    strain_map.columns = ['well','name']
    melted = sunrise.reset_index().melt(id_vars='Time')
    melted.variable = tools.numbers_to_name(melted.variable)
    labeled = pd.merge(left=melted, right = strain_map,
                       left_on='variable', right_on='well')
    labeled=labeled.drop(columns=['variable'])
    labeled = labeled.rename(columns={'value':'OD595'})
    labeled.OD595 = labeled.OD595.astype(float)
    if raw_time:
        return labeled
    labeled.Time = labeled.Time.apply(lambda x: x.total_seconds()/3600)
    return(labeled)
    
def curve_viewer(sheet, context='talk', hue='name', order=False):
    sns.set_context(context)
    sheet = sheet.rename(columns={'name': hue})
    if order:
         sns.relplot(data = sheet, x='Time', y='OD595', kind = 'line',
                     hue = 'name', hue_order=order, aspect=2)
    else:
        sns.relplot(data = sheet, x='Time', y='OD595', kind = 'line',
                     hue = 'name', aspect=2)
def mu_max(curves, norm_eqs=None):
    if norm_eqs:
        for strain, norm_eq in norm_eqs.items():
            try:
                curves.update(curves.assign(OD595norm = norm_eq(curves[curves.name.str.contains(strain)].OD595)))
            except(Exception):
                curves = curves.assign(OD595norm = norm_eq(curves[curves.name.str.contains(strain)].OD595))
        curves.OD595 = curves.OD595norm
    #curves = curves.reset_index()
    #curves.rename(columns={'index': 'indexer'})
    window_size = 12
    curves.Time = pd.TimedeltaIndex(curves.Time, unit='h').round('T')
    data = (curves.set_index(['Time', 'name', 'well']).unstack([1,2]).resample('5T').mean()\
    .rolling(window_size).apply(lambda x: np.log(x[-1]/x[0])).OD595['2.5 hours':'15 hours']\
    .max()/(window_size/12)).drop('BLK').reset_index()
    data = data.assign(Strain = data.name.apply(lambda x: x.split(' ')[0]), 
                       Treatment = data.name.apply(lambda x: ' '.join(x.split(' ')[1:])))
    #data.pH = pd.to_numeric(data.pH, errors='ignore')
    data = data.rename(columns={0: 'Max growth rate'})
    sns.barplot(data=data, y='Treatment', x='Max growth rate', hue='Strain')
    return data

def mu_max_plot(growth_rates, strains, comparison='Cen'):
    #Strain = namedtuple('strain', ['key', 'name', 'marker', 'marker_size'])
    #strains = {'Cen': Strain('Cen','CEN.PK 113-7D',marker='s', marker_size=8),
    #       '03B': Strain('03B', 'Scheffersomyces amazonensis',marker='s', marker_size=8)}
    growth_rates = growth_rates[~growth_rates.Treatment.isnull()]
    fig, ax = plt.subplots()

    for key, strain in strains.items():
        straindata = growth_rates[growth_rates.Strain==key].groupby('Treatment')
        ax.errorbar(
            x=straindata.mean().index.values,
            y=straindata.mean().values,
            label=strain.name,
            marker=strain.marker,
            markersize=strain.marker_size,
            yerr=straindata.std().values)
    #plt.plot(df_pH.pH.values, df_pH.Cen.values, label='CEN.PK 113-7D', marker='s')
    #plt.plot(df_pH.pH.values, df_pH['03B'].values, label='Scheffersomyces amazonensis', marker=(8,1,45), markersize=13)
    ax.set(
        xlabel = 'Treatment',
        ylabel = '$\mu_{max}$',
        ylim = [0, 0.3],
        title = 'Growth rates across pH\n(YNB + 20 g/L AA)')
    ax.legend(bbox_to_anchor=(.5,1), loc='upper left', fontsize=9)
    plt.savefig(' '.join(growth_rates.Strain.unique())+'.png', bbox_inches='tight', dpi=300)

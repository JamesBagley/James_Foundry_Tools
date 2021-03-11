# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:58:48 2020


Version 1.1 
@author: James
"""
from microtiter.plate_reader_tools import read_infinateM200_output
import flowio
from pandas import DataFrame, concat
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, ceil

def getwell(well, folder):
    if len(well)<3:
        well = well[0]+'0'+well[1]
    return folder+'\\'+well+'.fcs'
    
def getstrain(strain, strain_map):
    S_map =  read_infinateM200_output(strain_map)
    S_map = dict(zip(S_map.wells.values(), S_map.wells.keys()))
    return(S_map[strain])
    
def getflowdata(strain):
    flowdata = flowio.flowdata.FlowData(getwell(getstrain(strain)))
    columns = list(map(lambda x: x['PnS'], flowdata.channels.values()))
    data = np.reshape(flowdata.events, (-1, flowdata.channel_count))
    return(DataFrame(data=data, columns=columns))
    
    
def multistrain(strains, folder, X='FSC-A', Y='FL2-A', xmax=3000000, ymax=30000, dim=None):
    if dim:
        fig, axes = plt.subplots(*dim)
    else:
        dim = [ceil(sqrt(len(strains))), ceil(sqrt(len(strains)))]
        fig, axes = plt.subplots(*dim)
    
    plot_loc = dict()
    for i in range(dim[0]*dim[1]):
        plot_loc[i] = [i//dim[1], i%dim[1]]
    for i, strain in enumerate(strains):
        data = getflowdata(strain, folder)

        if len(X.split('/'))>1:
            Xs = X.split('/')
            Xdata = data[Xs[0]]/data[Xs[1]]
            Xdata.name = X
        elif len(X.split('*'))>1:
            Xs = X.split('*')
            Xdata = data[Xs[0]]/data[Xs[1]]
            Xdata.name = X
        else:
            Xdata = data[X]
        
        if len(Y.split('/'))>1:
            Ys = Y.split('/')
            Ydata = data[Ys[0]]/data[Ys[1]]
            Ydata.name = Y
        
        elif len(Y.split('*'))>1:
            Ys = Y.split('/')
            Ydata = data[Ys[0]]/data[Ys[1]]
            Ydata.name = Y
            
        else:
            Ydata = data[Y]
        
        data = concat([Xdata, Ydata], axis='columns')

        data = data[(data[X]<xmax) & (data[Y]<ymax)]
        ax = axes[plot_loc[i][0], plot_loc[i][1]]
        ax.scatter(data[X], data[Y], s=0.05)
        ax.set_title(strain)
        
        #plt = sns.kdeplot(data[X], data[Y], shade=True, n_levels=    30)
        xlim=[0,xmax]
        ylim=[0,ymax]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)        
        ax.set_xlabel(X)
        ax.set_ylabel(Y)
    plt.show()

def onestrain(strain, folder, X='FSC-A', Y='FL2-A', xmax=3000000, ymax=30000, plot=True):
    data = getflowdata(strain, folder)

    if len(X.split('/'))>1:
        Xs = X.split('/')
        Xdata = data[Xs[0]]/data[Xs[1]]
        Xdata.name = X
    elif len(X.split('*'))>1:
        Xs = X.split('*')
        Xdata = data[Xs[0]]/data[Xs[1]]
        Xdata.name = X
    else:
        Xdata = data[X]
    
    if len(Y.split('/'))>1:
        Ys = Y.split('/')
        Ydata = data[Ys[0]]/data[Ys[1]]
        Ydata.name = Y
    
    elif len(Y.split('*'))>1:
        Ys = Y.split('/')
        Ydata = data[Ys[0]]/data[Ys[1]]
        Ydata.name = Y
        
    else:
        Ydata = data[Y]
    
    data = concat([Xdata, Ydata], axis='columns')
    data = data[(data[X]<xmax) & (data[Y]<ymax)]
    if plot:
        plt.scatter(data[X], data[Y], s=0.05)
        plt.title(strain)
        
        #plt = sns.kdeplot(data[X], data[Y], shade=True, n_levels=    30)
        xlim=[0,xmax]
        ylim=[0,ymax]
        plt.xlim(xlim)
        plt.ylim(ylim)    
        plt.xlabel(X)
        plt.ylabel(Y)
        plt.show()
    
    return data


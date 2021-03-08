# -*- coding: utf-8 -*-
"""
Created on Fri May 25 15:54:17 2018

@author: James
"""

import pandas as pd
from microtiter import plate_reader_tools
import sys
import random

class Plate:
    def __init__(self, name):
        self.name = name
        self.wells = {}
    
    def new_well(self, name, components):
        if not type(components) == dict:
            return(TypeError('components must be of type dict'))
        
        
class Media_Well:
    def __init__(self, name, number):
        self.name = name
        self.number = number

FILE = r'C:\Users\Owner\Documents\Concordia\1000 yeast strains\Biomek instructions\hydrolysate test 2 map.xlsx'

plate_maps_excel = pd.ExcelFile(FILE)

plates = plate_reader_tools.find_plates(plate_maps_excel)

for plate in plates:
    if 2 in plate.shape:
        sys.exit('plate wrong shape')

for plate in plates:
    plate.set_index(plate.iloc[0:, 0], inplace=True)
    plate.drop(['<>'], axis='rows', inplace=True)
    plate.drop([0], axis='columns', inplace=True)


concentrations = []
for i, plate in enumerate(plates):
    print(plate)
    concentrations.append(input(
            'plate {} looks like this, what is the source concentration in g/L?\
            \n'.format(i+1)))

melted_plates = []
for plate in plates:
    melted_plates.append(plate.T.melt())

melted_lists = []
for melted_plate in melted_plates:
    for well in melted_plate.index:
        value = float(melted_plate.iloc[well]['value'])

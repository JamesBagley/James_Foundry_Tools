# -*- coding: utf-8 -*-
"""
Created on Wed May  9 15:10:37 2018

@author: James
"""

import pandas as pd
import numpy as np
import datetime 

## Reads ODs from the 3rd floor plate reader, supports both matrix and list
## output formats

def read_infinateM200_output(target):
    file = pd.read_excel(target, 0, header=None)
    plate = find_plate(file)

    print(plate.shape)
    if plate.shape[1] > 2:
        MODE = 'MATRIX'
    else:
        MODE = 'LIST'

    print('Mode =', MODE)

    if MODE == 'MATRIX':
        plate_dict = {}
        counter = 0
        plate.set_index(plate[0], inplace=True)
        plate.drop(['<>'], axis='rows', inplace=True)
        plate.drop([0], axis='columns', inplace=True)
        for row in plate.index.values:
            for col in plate.columns.values:

                index = (str(row) + str(col))

                try:
                    value = float(plate[col][row])
                except ValueError:
                    value = plate[col][row]

                plate_dict[counter] = {'Value': value,
                                       'Well': index}
                counter += 1

        plate = pd.DataFrame.from_dict(plate_dict).T
        plate.set_index(plate['Well'], inplace=True)
        plate.drop(['Well'], axis='columns', inplace=True)

    plate.iloc[0, 0] = 'Well'
    if MODE == 'LIST':
        plate = plate.set_index(plate[0])
        plate.columns = plate.iloc[0]
        plate = plate.drop(['Well'])
        plate = plate.drop(['Well'], axis=1)
    return plate



##  Reads the xlsx outputs from tecan sunrise plate readers, see example file
##  "C:\Users\jbag2\OneDrive - Concordia University - Canada\Dropbox\Lab_Automation\example_file_structures\sunrise_growth_curve.xlsx"
##  if legth mismatch error increase number of rows skipped, e.g. need to skip
##  4 rows for example file

def read_sunrise(target, skipbottomrows=-1):
    
    if target == 'testing':
        target =\
        r"C:\Users\Owner\Concordia\Lab_Automation\example_files\sunrise_growth_curve.xlsx"
    curves = pd.read_excel(target, skiprows=1, header=None).T
    curves = curves.iloc[0:-2, 0:skipbottomrows]
    
    curves.set_index(curves[0].apply(lambda x: datetime.timedelta(seconds=int(x[:-1]))),inplace=True)
    curves.drop(0, axis='columns', inplace=True)
    curves.drop(1, axis='columns',inplace=True)
    curves.columns = np.arange(1, 97)
    return curves



##  takes output from read_sunrise and applies a rolling function to return 
##  growth rates across all wells at chosen interval, default is 60 minutes. 

def growth_rate(sunrise, window_size=60):
    _mininsec = 1/60
    _ws = str(window_size)+'min'
    def growth_func(x):
        try:
            rate = np.log2(x[-1]/x[0])/(((x.index[-1]-x.index[0])).seconds*(1/3600))
        except(ZeroDivisionError):
            rate = np.nan
        return rate
    return(sunrise.rolling(_ws).apply(growth_func))


##  mini function for finding plate data from infinite_M200 due to weird 
##  output file inconsistencies.

def find_plate(dataframe):
    origin = [np.argmax(dataframe[0].str.find('<>')), 0]
    final_col, final_row = find_dimensions(origin, dataframe)
    ydim = [origin[0], final_row]
    xdim = [origin[1], final_col]
    plate_df = dataframe.iloc[ydim[0]:ydim[1], xdim[0]:xdim[1]]

    return plate_df

##  same as find_plate but functions across multi sheet excell files
##
def find_plates(excel):
    plates = []
    for sheet_name in excel.sheet_names:
        try:
            sheet = excel.parse(sheet_name, header=None)
            plate = find_plate(sheet)
            plates.append(plate)
            print('success')
        except(AttributeError):
            print('failed')
            pass
    return plates

##  reads an excel file with the same structure as an infinite M200 plate
##  to get sample names 
def read_treatment_map(target):
    if target == 'testing':
        target = r'C:\Users\Owner\Concordia\Lab_Automation\example_files\treatment_map.xlsx'
    file = pd.read_excel(target, 0, header=None)
    treatment_map = find_plate(file)

    MODE = matrix_or_list(treatment_map)

    print('Mode =', MODE)

    if MODE == 'MATRIX':
        treatment_map_dict = {}
        treatment_map.set_index(treatment_map[0], inplace=True)
        treatment_map.drop(['<>'], axis='rows', inplace=True)
        treatment_map.drop([0], axis='columns', inplace=True)
        counter = 0
        for row in treatment_map.index.values:
            for col in treatment_map.columns.values:

                index = (str(row) + str(col))
                try:
                    value = float(treatment_map[col][row])
                except ValueError:
                    value = treatment_map[col][row]
                treatment_map_dict[counter] = {'Value': value,
                                               'Well': index}
                counter += 1
        treatment_map = pd.DataFrame.from_dict(treatment_map_dict).T
        treatment_map.set_index(treatment_map['Well'], inplace=True)
        treatment_map.drop(['Well'], axis='columns', inplace=True)

    if MODE == 'LIST':
        treatment_map = treatment_map.set_index(treatment_map[0])
        treatment_map.columns = treatment_map.iloc[0]
        treatment_map = treatment_map.drop(['Well'])
        treatment_map = treatment_map.drop(['Well'], axis=1)
    return treatment_map

##
def matrix_or_list(dataframe):
    if dataframe.shape[1] > 2:
        return 'MATRIX'
    else:
        return 'LIST'


def find_dimensions(origin, dataframe):
    final_column = find_dimension('cols', dataframe, origin)
    final_row = find_dimension('rows', dataframe, origin)

    return ([final_column, final_row])


def find_dimension(mode, dataframe, origin, previous=None):
    modes = ['rows', 'cols']
    axis = modes.index(mode)
    if not previous:
        previous = origin.copy()

    previous[axis] += 1
    current = previous
    current_value = str(dataframe.iloc[current[0], current[1]])

    if current_value == 'nan':
        return current[axis]
    else:
        try:
            return find_dimension(mode, dataframe, origin, current)
        except(IndexError):
            return current[axis]

def as_matrix(data, function=False, rows=list('ABCDEFGH'),
              columns=list(range(1, 13))):
    '''
    converts specifed rows and columns of the plate into a matrix, values
    are well objects by default but can accept a function arguement to get
    specific values such as the last OD reading
    '''
    matrix = []
    for row_index in rows:
        row = []
        for well in data.keys():
            well = data[well]
            if well.name[0] == row_index and int(well.name[1:]) in columns:
                # if a function is passed, run that function on the
                # well, else return the well object
                if function:
                    # if the function returns a not_a_number, return 0
                    if np.isnan(function(well)):
                        row.append(0)
                    # else return the result of the function
                    else:
                        row.append(function(well))
                else:
                    row.append(well)
        matrix.append(row)
    return matrix

## Converts well numbers to well position names
def numbers_to_name(numbers):
    if type(numbers) == int:
        row,col = numbers//12, (numbers-1)%12+1
        row = 'abcdefgh'[row].upper()
        return(row+str(col))
    else:
        names = []
        for num in numbers:
            num -= 1
            row,col = num//12, (num)%12+1
            row = 'abcdefgh'[row].upper()
            names.append(row+str(col))
        return names

## Converts well position names to numbers
def names_to_numbers(names):
    if type(names) == str:
        row = ord(names[0].lower())-97
        col = int(names[1:])
        return row*12+col
    else:
        numbers = []
        for name in names:
            row = ord(name[0].lower())-97
            col = int(name[1:])
            numbers.append(row*12+col)
        return numbers
            

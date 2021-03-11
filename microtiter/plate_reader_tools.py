# -*- coding: utf-8 -*-
"""
Created on Wed May  9 15:10:37 2018

@author: James
"""

import pandas as pd
import numpy as np
import datetime

class plate_96():
    def __init__(self):
        self.wells = {name: None for name in numbers_to_name(list(range(1,97)))}
    
    @property
    def matrix(self):
        return \
        pd.DataFrame(data=[self.wells.keys(), self.wells.values()]).T\
        .rename(columns={0:'well', 1:'value'})\
        .assign(row=lambda x: x.well.apply(lambda x: x[0])).\
        assign(column=lambda x: x.well.apply(lambda x: int(x[1:]))).\
        pivot(index='row', columns='column', values='value')
    
    @property
    def dataframe(self):    
        return\
        pd.DataFrame(index=self.wells.keys(), data=self.wells.values())\
        .rename(columns={0:'value'})
        
    def view(self):
        print(self.matrix.to_string())

def numbers_to_name(numbers):
    if type(numbers) == int:
        num = numbers - 1
        row, col = num // 12, num % 12 + 1
        row = 'abcdefgh'[row].upper()
        return(row+str(col))
    else:
        names = []
        for num in numbers:
            num -= 1
            row, col = num // 12, (num) % 12+1
            row = 'abcdefgh'[row].upper()
            names.append(row+str(col))
        return names

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
            
def read_infinateM200_output(target):
    file = pd.read_excel(io=target, header=None)
    plate = find_plate(file)
    plate = read_plate(plate)
    
    plate_obj = plate_96()
    for item in plate.iteritems():
        plate_obj.wells[item[0]] = item[1]
    
    
    return plate_obj


def read_sunrise(target):
    if target == 'testing':
        target =\
        r"C:\Users\Owner\Concordia\Lab_Automation\example_files\sunrise_growth_curve.xlsx"
    excel = pd.ExcelFile(target)
    raw_sheet = excel.parse()
    column_count = len(raw_sheet.columns)
    row_count = len(raw_sheet)

    curves = excel.parse().drop(columns=['Replicate Info', 'Well positions'],
                        index=[0,1])
    
    skip_rows = len(curves[pd.isna(curves.iloc[:,1])])
    curves = curves.iloc[:-skip_rows]
    
    headers = excel.parse(usecols=[column_count-2])
    
    #curves = pd.read_excel(target, skiprows=1, header=None).T
    times = excel.parse(skiprows=1, skipfooter=row_count-1, header=None).\
        drop(columns=[column_count-2, column_count-1]).values[0]
    times = list(map(lambda x: datetime.timedelta(seconds=int(x[:-1])),
                     times.tolist()))
    curves.columns = times

    plate_details = curves.iloc[-2:]
    run_details = curves.iloc[0, -3:]
    #date_details = curves.iloc[0, -4]
    #date_details = date_details.replace('/', ': ').split(': ')
    #date_details = [date_details[1], date_details[3]]
    #start_date = date_details[0].split('-')
    #start_date = datetime.date(*[int(val) for val in start_date])

    #start_time = date_details[1].split(':')
    #start_time = datetime.time(*[int(val) for val in start_time])

    #start_datetime = datetime.datetime.combine(start_date, start_time)
    #temperatures = curves.iloc[0:,0]

    #curves.columns = np.arange(1, 97)
    curves = curves.set_index(pd.Index(names_to_numbers(excel.parse()['Well positions'].dropna())))

    return curves.T


def read_plate(plate):
    if plate.shape[1] > 2:
        MODE = 'MATRIX'
    else:
        MODE = 'LIST'

    print('Mode =', MODE)

    if MODE == 'MATRIX':
        out_series = pd.Series()
        out_series.name = 'Well'

        plate.set_index(plate[0], inplace=True)
        plate.drop(['<>'], axis='rows', inplace=True)
        plate.drop([0], axis='columns', inplace=True)

        for row in plate.index.values:
            for col in plate.columns.values:
                index = (str(row) + str(col))
                value = plate[col][row]
                try:
                    value = float(value)
                except ValueError:
                    pass
                out_series = out_series.append(pd.Series(value, [index]))

    if MODE == 'LIST':
        plate.iloc[0, 0] = 'Well'
        plate = plate.set_index(plate[0])
        plate.columns = plate.iloc[0]
        plate = plate.drop(['Well'])
        plate = plate.drop(['Well'], axis=1)
        out_series = pd.Series(plate.iloc[:,0])

    return out_series


def find_plate(dataframe):
    origin = [np.argmax(dataframe[0].str.find('<>')), 0]
    final_col, final_row = find_dimensions(origin, dataframe)
    ydim = [origin[0], final_row]
    xdim = [origin[1], final_col]
    plate_df = dataframe.iloc[ydim[0]:ydim[1], xdim[0]:xdim[1]]

    return plate_df


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

def read_excel_sheets(file, names=None):
    excel = pd.ExcelFile(file)
    platelist = find_plates(excel)
    for i, plate in enumerate(platelist):
        if i == 0:
            output = read_plate(plate).rename(i)
        else:
            output = pd.merge(output, read_plate(plate).rename(i),
                              left_index=True, right_index=True)
    if names:
        names = dict(zip(range(len(names)),names))
        output = output.rename(columns=names)
    output.to_clipboard()
    return(output)

def structured_find_plates(excel):
    config = excel.parse('Config', header=None, index_col=0)
    sheets = config.loc['sheets']

    plates = {}
    for sheet in sheets:
        print(sheet)
        plates[sheet] = find_plate(excel.parse(sheet, header=None))
        
    return plates


def read_treatment_map(target):
    if target == 'testing':
        target = r'C:\Users\Owner\Concordia\Lab_Automation\example_files\treatment_map.xlsx'
    file = pd.read_excel(target, 0, header=None, dtype=str)
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
    are well objects by default but can accept a function argument to get
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
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:54:31 2018



@author: James
"""
import pandas as pd
import os.path

DEFAULT_OUTPUT_ORDER = ['culture_plate', 'culture_well', 'dilutant_plate',
                            'dilutant_well','destination', 'destination_well',
                            'culture_volume', 'dilutant_volume','row',
                            'column']

def media_plate(volume, media, wells=96):
    return None


def dilution(dilution_informtion,
             filename='test.csv',
             save_file=True):

    SOURCE_PLATE = 'CELL PLATE'
    DESTINATION_PLATE = '=DESTINATION'
    DILUTANT = 'MEDIA'
    DILUTANT_WELL = 4
    

    out_dict = {}
    for well in dilution_informtion.keys():
        out_dict[well] = {
                            'culture_plate': SOURCE_PLATE,
                            'culture_well': well,
                            'dilutant_plate': DILUTANT,
                            'dilutant_well': DILUTANT_WELL,
                            'destination': DESTINATION_PLATE,
                            'destination_well': well,
                            'culture_volume': dilution_informtion[well]['transfer_volume'],
                            'dilutant_volume': dilution_informtion[well]['dilutant_volume'],
                            'row': well[0],
                            'column': int(well[1:])
                          }    

    out_dataframe = pd.DataFrame.from_dict(out_dict, orient='index')
    out_dataframe.sort_values(by=['column', 'row'], inplace=True)
    
    
    out_dataframe = reorder_dilution_columns(out_dataframe)
    out_dataframe.drop(['row','column'], axis=1, inplace=True)
    if filename[-4:] != '.csv':
        filename += '.csv'
        
    filename = (os.path.join('.','output_files',filename))
    if save_file:
        out_dataframe.to_csv(filename)
    else:
        return(out_dataframe)


def matrix_to_row(matrix_id):

    LETTER_EQUIV =\
        {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}

    return(LETTER_EQUIV[matrix_id[0].upper()]*12+int(matrix_id[1:]))

def file_writer(dictionary):
    pass


def plate_dict_writer():
    pass

def reorder_dilution_columns(dataframe, desired_order = 'default'):
    if desired_order == 'default':
        desired_order = DEFAULT_OUTPUT_ORDER
        
    reordered_dataframe = dataframe[desired_order]
    return(reordered_dataframe)
    
    pass


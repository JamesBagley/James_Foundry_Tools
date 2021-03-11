# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:54:31 2018



@author: James Andrew Bagley
"""
import pandas as pd
import os.path


DEFAULT_OUTPUT_ORDER = ['culture_plate', 'culture_well', 'dilutant_plate',
                            'dilutant_well','destination', 'destination_well',
                            'culture_volume', 'dilutant_volume','row',
                            'column']

def dilution(dilution_information,
             filename='test.csv',
             save_file=True,
             randomize=False,
             inoculant_vol=10,
             intermediate_vol=200):

    SOURCE_PLATE = 'CELL PLATE'
    DESTINATION_PLATE = 'Intermediate'
    DILUTANT = 'MEDIA'
    DILUTANT_WELL = 4
    
    for key in dilution_information.keys():
        dilution_information[key]['destination'] = dilution_information[key]

    out_dict = {}
    for name, well in dilution_information.items():
        if max (well['transfer_volume']/inoculant_vol*intermediate_vol,
                well['dilutant_volume']/inoculant_vol*intermediate_vol) <= 200:
            out_dict[name] = {
                            'culture_plate': SOURCE_PLATE,
                            'culture_well': name,
                            'dilutant_plate': DILUTANT,
                            'dilutant_well': DILUTANT_WELL,
                            'destination': DESTINATION_PLATE,
                            'destination_well': name,
                            'culture_volume': well['transfer_volume']/inoculant_vol*intermediate_vol,
                            'dilutant_volume': well['dilutant_volume']/inoculant_vol*intermediate_vol,
                            'row': name[0],
                            'column': int(name[1:]),
                            'platerowcol':SOURCE_PLATE+DESTINATION_PLATE+name[1:]
                          }
        else:
            transfer_num = 0
            vals = [well['transfer_volume']/inoculant_vol*intermediate_vol,
                   well['dilutant_volume']/inoculant_vol*intermediate_vol]
                
            while max (vals) > 0:
                
                transfer_num += 1
                out_dict[name+' - '+str(transfer_num)] = {
                            'culture_plate': SOURCE_PLATE,
                            'culture_well': name,
                            'dilutant_plate': DILUTANT,
                            'dilutant_well': DILUTANT_WELL,
                            'destination': DESTINATION_PLATE,
                            'destination_well': name,
                            'culture_volume': min(vals[0], 200),
                            'dilutant_volume': min(vals[1], 200),
                            'row': name[0],
                            'column': int(name[1:]),
                            'platerowcol':SOURCE_PLATE+DESTINATION_PLATE+name[1:]}
                
                for i, val in enumerate(vals):
                    vals[i] = max(val-200, 0)
                
            
    out_dataframe = pd.DataFrame.from_dict(out_dict, orient='index')
    out_dataframe.sort_values(by=['column', 'row'], inplace=True)
    
    if filename[-4:] != '.csv':
        filename += '.csv'
        
    filename = (os.path.join('..', 'output_files', filename))
    if save_file:
        out_dataframe.to_csv(filename)

    else:
        return(out_dataframe)


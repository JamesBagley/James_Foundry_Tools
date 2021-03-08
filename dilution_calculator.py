# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:37:43 2018

@author: James
"""

import numpy as np
import sys

TOTAL_INOCULATION_VOLUME = 10
TARGET_OD = 0.1


def from_plate_reader(plate_df, inoculant_volume=TOTAL_INOCULATION_VOLUME,
                      target_OD=TARGET_OD):
    print(plate_df['Value'].values)
    smallest = np.min(plate_df['Value'].values)
    if smallest < target_OD:
        print('Lowest OD reading is lower than current target OD({})'.format(
                target_OD))
        target_equals_smallest = input('Would you like to set the target OD ' +
                                       'to equal the smallest OD? \n' +
                                       'y/n \n')
        target_equals_smallest = target_equals_smallest in \
            ('Y', 'y', 'yes', 'Yes')
        if target_equals_smallest:
            target_OD = smallest
        else:
            sys.exit('Target OD too high')

    dilutions = {}
    for well in plate_df.index.values:
        dilution_information = {}

        dilution_factor = target_OD/plate_df.loc[well]['Value']
        transfer_volume = dilution_factor * inoculant_volume
        dilutant_volume = inoculant_volume - transfer_volume

        dilution_information = {'dilution_factor': dilution_factor,
                                'transfer_volume': transfer_volume,
                                'dilutant_volume': dilutant_volume}

        dilutions[well] = dilution_information

    return dilutions

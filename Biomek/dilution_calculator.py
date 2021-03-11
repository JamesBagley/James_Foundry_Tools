# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:37:43 2018

@author: James Andrew Bagley
"""

from warnings import warn

DEFAULT_INOCULANT_VOLUME = 10
DEFAULT_TARGET_VOLUME = 180
DEFAULT_TARGET_OD = 0.1
DEFAULT_MEASURED_DIL = 1 # 1 = not diluted, 10 = diluted 1/10


def from_plate_reader(*, plate_df,
                      inoculant_volume=DEFAULT_INOCULANT_VOLUME,
                      target_volume=DEFAULT_TARGET_VOLUME,
                      target_OD=DEFAULT_TARGET_OD):
    
    if not inoculant_volume:
        inoculant_volume = DEFAULT_INOCULANT_VOLUME

    dilutions = {}
    for well, value in plate_df.value.items():
        dilution_information = {}
        dilution_factor = target_OD/(value)
        transfer_volume = dilution_factor * target_volume
        if transfer_volume > inoculant_volume:
            print(transfer_volume, inoculant_volume)
            warn(
                    '{} OD too low! Culture transfer set to max volume'.format(well),
                    category = Warning)
            transfer_volume = inoculant_volume
        dilutant_volume = inoculant_volume - transfer_volume

        dilution_information = {'dilution_factor': dilution_factor,
                                'transfer_volume': transfer_volume,
                                'dilutant_volume': dilutant_volume}

        dilutions[well] = dilution_information
    return dilutions

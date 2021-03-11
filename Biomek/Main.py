# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:49:03 2018

@author: James Andrew Bagley

PRINCIPLE:
    To start a growth curve experiment with a number of different strains
    the starting ODs of each strain have to be brought in line as the 
    precultures will naturally deviate. As the number of strains increases,
    this becomes a non-trivial problem. Automating this step makes it feasible
    to conduct experiments on much larger scales, testing >96 strains 
    simultaneously.
    
    To do this you need:
        A) accurate OD readings for each strain
        B) sufficient cell counts to reach the target OD in each strain
    
    The solution provided assumes that conditions are on different plates and
    as a result it is designed to create an "intermediate plate" in which all
    wells have the same OD, the intermediate plate is then used to inoculate 
    the experiment plates.
    
    1 preculture plate -> 1 intermediate plate -> n experimental plates
    where n = number of conditions 
    or if replicates are on different plates:
        n = number of conditions * number of replicates 
        
    
    
    
    
HOW TO OPERATE:
    
    1) Take a sample from all 96 wells of preculture and dilute til OD is 
       approximately 0.4 (when measuring in a microtiter plate with a total
       volume of 180uL)
       
       NOTE: 0.3-0.5 are the ideal OD reading range as they are in the most
       linear range of the OD standard curve (in a plate with 180uL of liquid)
              
    2) Measure ODs
       
    3) Using the "equal_OD_dilutions" function input all relevent paramaters:
           measurement dilution,
           target starting OD in experiment,
           Volume of experiment culture (including inoculum),
           volume of inoculant into experiment culture,
           volume of intermediate plate,
           location of wells that contain blanks (used for OD calculations),
           location of wells that should be ignored 
           (ignored wells are made into blanks in the intermediate plate),
           path for the OD measurement file generated in 2,
           desired filename for the generated biomek instructions file.
           
    4) Take biomek instructions file and read using a biomek transfer from file
       step.
           Diluent and culture transfers are stored in seperate columns and
           need to be handled by different transfer from file steps as they use
           cell culture transfers and media transfers should use different 
           liquid types and pipetting techniques
        
    
"""

from microtiter.plate_reader_tools import read_infinateM200_output
import dilution_calculator
import biomek_file_writer


def equal_OD_dilutions(measurement_dilution,
                       target_OD=0.1,
                       target_volume=180,
                       inoculant_volume=False,
                       file_in=False,
                       file_out=False,
                       test_mode=False,
                       ignore_wells=set(),
                       blank_wells=set(),
                       intermediate_volume=200):

    TESTING_FILE = \
       r'C:\Users\Owner\Concordia\Lab_Automation\example_files\plate_reader_file_format.xlsx'
    DEFAULT_INOCULANT_VOLUME = 10
    ignore_wells.update(blank_wells)

    if test_mode:
        if not file_in:
            file_in = TESTING_FILE
        if not file_out:
            file_out = 'test'
        if not inoculant_volume:
            inoculant_volume = DEFAULT_INOCULANT_VOLUME
    else:
        if not file_in:
            file_in = input("what's the address for your tecan file?\n")
        if not file_out:
            file_out = input("what should the name of the output be?\n")

    plate = read_infinateM200_output(file_in).dataframe
    if blank_wells:
        blank_val = plate.loc[blank_wells].mean()
        print('Blank average = {}'.format(blank_val))
        plate = plate - blank_val
    plate = plate * measurement_dilution
    print(plate)
    too_low = False #set(plate.where(lambda x: x<0.1).index.values).difference(ignore_wells)
    
    if too_low:
        print('Well ODs less than 0.1!', too_low, 'Wells dropped!', sep='\n')
        ignore_wells.update(too_low)
    
    dilutions_df = plate.drop(ignore_wells)
    dilution_information = dilution_calculator.from_plate_reader(
            plate_df=dilutions_df,
            inoculant_volume=inoculant_volume,
            target_OD=target_OD,
            target_volume = target_volume)
    for well in ignore_wells:
        dilution_information[well] = {'dilution_factor': 0,
                                      'transfer_volume': 0,
                                      'dilutant_volume': inoculant_volume}
    biomek_file_writer.dilution(dilution_information,
                                filename=file_out,
                                inoculant_vol=inoculant_volume,
                                intermediate_vol=intermediate_volume)

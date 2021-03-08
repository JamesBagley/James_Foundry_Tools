# -*- coding: utf-8 -*-
"""
Created on Thu May 10 18:49:03 2018

@author: James
"""

from microtiter import plate_reader_tools
import dilution_calculator
import biomek_file_writer


def equal_OD_dilutions(file_in=False,
                       file_out=False,
                       innoculation_vol=False,
                       target_OD=0.1,
                       test_mode=False):

    TESTING_FILE = \
       r'C:\Users\Owner\Concordia\Lab_Automation\example_files\plate_reader_file_format.xlsx'
    DEFAULT_INOCULATION_VOLUME = 10
    if test_mode:
        if not file_in:
            file_in = TESTING_FILE
        if not file_out:
            file_out = 'test'
        if not innoculation_vol:
            innoculation_vol = DEFAULT_INOCULATION_VOLUME

    else:
        if not file_in:
            file_in = input("what's the address for your tecan file?\n")
        if not file_out:
            file_out = input("what should the name of the output be?\n")

    plate = plate_reader_tools.read_infinateM200_output(file_in)

    dilution_information = dilution_calculator.from_plate_reader(
            plate, innoculation_vol, target_OD)
    biomek_file_writer.dilution(dilution_information, filename=file_out)

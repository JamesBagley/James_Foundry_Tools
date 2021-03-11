# -*- coding: utf-8 -*-
"""
Created on Sat May 26 20:44:08 2018

TODO:
Put in actual math for find_phases

@author: James
"""

import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
import seaborn as sns
import Concordia.Lab_Automation.microtiter.plate_reader_tools as tools
import pandas as pd


class Plate():
    '''
    Generates a list of well objects, each storing their growth curve, a spline
    representation along with the corresponding derivative and some features
    for easier statistical analysis.

    Has methods to generate heatmaps and a matrix representation of the plate
    TODO: implement better replicate handling
    '''

    def save(self, filename=None):
        import pickle
        
        if not filename:
            filename = '{}.pickle'.format(self.name)
        file = open(filename, 'wb')
        pickle.dump(self, file)
        file.close()

    def __init__(self, curves, name,
                 dims =[list('ABCDEFGH'),list(range(1, 13))]):
        '''
        initiating a new plate object, takes the pandas dataframe generated
        from the read sunrise module in plate_reader_tools and a string as a
        name for the plate.

        Generates a dictionary of every well in the plate using the well class
        and stores the original dataframe
        '''

        self.well_list = {}
        self.sample_list = {}
        self.blank_list = {}
        self.sunrise_dataframe = curves
        self.name = name
        self.dimensions = dims
        for curve in curves:
            self.well_list[curve] = self.Well(curve, curves[curve])

    def __getitem__(self, key):
        try:
            return [self.well_list[k] for k in key]
        except(TypeError):
            return self.well_list[key]
        except(KeyError):
            return(self.well_list[tools.names_to_numbers(key)])

    class Well():
        '''
        Class representing wells in a plate, calculates common statistics
        automatically
        '''
        def __init__(self, number, curve, blank=False, blanked_curve=False,
                     strain=False, species=False):
            self._curve = curve
            # self.integral, self.spline, self.derivative = \
            #     fit_spline(self.curve, der=[-1, 0, 1])
            # try:
            #     self.phases = find_phases(self.derivative)
            # except(IndexError):
            #     self.phases = [np.NaN]*4
            # try:
            #     self.lag_phase_length, self.mu_max, self.saturation_time =\
            #         read_curve(curve)
            # except(IndexError):
            #     self.lag_phase_length, self.mu_max, self.saturation_time =\
            #         [np.NaN]*3
            # if self.curve[-1]>self.curve[0]*1.25:
            #    try:
            #        index = curve.index.values/np.timedelta64(1,'m')
            #        self.rkslc = growth_curve_regression.get_growth_char(
            #                index, curve.values.astype(np.float64))
            #        self.mu = self.rkslc[0]
            #    except(RuntimeError):
            #        self.mu = 0
            # else:
            #    self.mu = 0
            # self.mu_spline = self.derivative/self.spline

            self._number = number
            self._treatments = {}
            self._blank = blank
            self._blanked_curve = blanked_curve
            self._strain = strain
            self._species = species

        @property
        def integral(self, blanked=True):
            return self.spline(blanked=blanked, der=-1)

        def spline(self, **kwargs):
            defaults = {'log': False, 'visualize': False, 'save': False,
                        'filename': False, 'der': 0, 'mu': False,
                        'blanked': False, 'resample': False}
            context = {}
            context.update(defaults)
            context.update(kwargs)
            resample = context.pop('resample')
            if context.pop('blanked'):
                curve = self.blank_curve
            else:
                curve = self.raw_curve
            if resample:
                curve = curve.astype(np.float16).resample(str(resample)+'T').mean()
                

            return fit_spline(curve, **context)

        @property
        def mu(self):
            try:
                if not self.blank_curve:
                    return Exception('plate must be blanked before taking mu')
            except(ValueError):
                pass
            mu_curve = self.spline(blanked=False, mu=True, resample=60)
            return mu_curve

        @property
        def blank(self, verbose=False):
            if verbose:
                try:
                    print('retrieving blank information for {}({})'
                          .format(self.name, self.number))
                except(NameError):
                    print('retrieving blank infromation')
            return self._blank

        @blank.setter
        def blank(self, blank):
            self._blank = blank

        @property
        def treatments(self):
            return self._treatments

        @treatments.setter
        def treatment(self, treatment):
            self._treatments[treatment[0]] = treatment[1]

        def color_code(self):
            phases = find_phases(self.derivative)
            color_phases(self.spline, phases)

        @property
        def raw_curve(self):
            return self._curve

        @raw_curve.setter
        def curve(*, self, curve):
            self._curve = curve

        @property
        def blank_curve(self):
            return self._blanked_curve

        @blank_curve.setter
        def blank_curve(self, blank_curve):
            self._blanked_curve = blank_curve

        @property
        def name(self):
            return tools.numbers_to_name(self._number)

        @property
        def number(self):
            return self._number

        @property
        def strain(self):
            return self._strain

        @strain.setter
        def strain(self, strain):
            self._strain = strain

        @property
        def species(self):
            return self._species

        @species.setter
        def species(self, species):
            self._species = species

    def get_replicates(self, rep_type, func):
        averaged_replicates = []
        stdev_replicates = []
        try:
            if rep_type.lower() == 'quadrant':
                cols_in_quad = len(self.dimensions[1])/2
                rows_in_quad = len(self.dimensions[0])/2
                treatments = cols_in_quad*rows_in_quad
                for treatment in range(1, int(treatments)+1):
                    wells = [treatment,
                             treatment+cols_in_quad,
                             treatment+treatments,
                             treatment+treatments+cols_in_quad]
                    wells = self[wells]
                    wells = np.array([func(wells[well]) for well in wells])
                    averaged_replicates.append(wells.mean())
                    stdev_replicates.append(wells.std())
        except(AttributeError):
            pass
        try:
            if rep_type[0].lower() == 'sequential':
                reps = rep_type[1]
                treatments =\
                    len(self.dimensions[0])*len(self.dimensions[1])/reps
                for treatment in range(int(treatments)):
                    wells = list(range(1+treatment*reps, 5+treatment*reps))
                    wells = self[wells]
                    wells = np.array([func(wells[well]) for well in wells])
                    averaged_replicates.append(wells.mean())
                    stdev_replicates.append(wells.std())
        except(AttributeError):
            pass

        replicate_info = list(zip(averaged_replicates, stdev_replicates))
        return(replicate_info)

    def set_replicates(self, treatment_name, rep_type, treatment_values=False):
        try:
            if rep_type.lower() == 'quadrant':
                cols_in_quad = len(self.dimensions[1])/2
                rows_in_quad = len(self.dimensions[0])/2
                treatments = int(cols_in_quad*rows_in_quad)
                if not(treatment_values):
                    treatment_values = range(1, treatments+1)
                treatments_values = list(zip(range(1, treatments+1), treatment_values))
                print(treatments)
                for treatment, value in treatments_values:
                    wells = [treatment,
                             treatment+cols_in_quad,
                             treatment+treatments,
                             treatment+treatments+cols_in_quad]
                    wells = self[wells]
                    for well in wells:
                        well.treatment = [treatment_name, value]
        except(AttributeError):
            pass
        try:
            if rep_type[0].lower() == 'sequential':
                reps = rep_type[1]
                treatments =\
                    int(len(self.dimensions[0])*len(self.dimensions[1])/reps)
                if not(treatment_values):
                    treatment_values = range(1, treatments+1)
                treatments_values = list(zip(range(treatments), treatment_values))
                for treatment, value in treatments_values:
                    wells = list(range(1+treatment*reps, reps+1+treatment*reps))
                    wells = self[wells]
                    for well in wells:
                        wells[well].treatment = [treatment_name, value]
        except(AttributeError):
            pass

        if rep_type[0].lower() == 'map':
            filename = rep_type[1]
            treatment_map = tools.read_treatment_map(filename)
            #treatment_name = treatment_map.columns[0]
            x =  zip(treatment_map.index.values, treatment_map['Value'].values)
            for well, treatment in x:
                try:
                    if '(B)' in treatment:
                        treatment = treatment[:-3]
                        self[well].blank = True
                        try:
                            treatment = float(treatment)
                        except(ValueError):
                            pass
                except(TypeError):
                    pass
                self[well].treatment = [treatment_name, treatment]

    def set_strains(self, filename):
        strain_map = tools.read_treatment_map(filename)
        x = zip(strain_map.index.values, strain_map['Value'].values)
        for well, strain in x:
            try:
                if strain < 99:
                    strain = '0' + str(int(strain))
                else:
                    strain = str(int(strain))
            except(TypeError):
                pass
            self[well].strain = strain

    def set_species(self, filename):
        species_series = pd.read_excel(
                filename,
                sheet_name='overview', header=0, usecols=[0, 2],
                index_col=0, skipfooter=3, squeeze=True)

        def species_func(well):
            try:
                well.species = species_series[well.strain]
            except(KeyError):
                pass
                #well.species = 'well.strain'

        wells = pd.DataFrame.from_dict(self.well_list, orient='index')
        wells.applymap(species_func)

    def __call__(self, funcs=False):
        out_dict = dict()
        for key, well in self.well_list.items():
            vals = list()
            for name, func in funcs.items():
                vals.append(func(well))
            if len(funcs) > 1:
                out_dict[well.name] = vals
            else:
                out_dict[well.name] = vals[0]

        if len(funcs) > 1:
            df = pd.DataFrame.from_dict(out_dict)
            df = df.T
            df.columns = funcs.keys()
            return df
        else:
            return pd.Series(out_dict)

    def correct_for_blank(self, blank_value=False):
        args = {'curve': lambda x: x.raw_curve,
                'treatments': lambda x: x.treatments,
                'blank': lambda x: x.blank}
        df = self(args)
        newcurve = pd.DataFrame([*df.pop('curve')])
        newcurve.columns = pd.MultiIndex.from_product([['curve'], newcurve.columns])
        newcurve.index = df.index
        df.columns = pd.MultiIndex.from_product([['well_info'], df.columns])
        df = df.join(newcurve)
        treatment_types = list(set(df['well_info']['treatments'].map(lambda x: tuple(x.keys()))))
        if len(treatment_types) > 1:
            print('NOT SUPPORTED')
            raise(Exception('mismatching treatments not supported'))
        # TODO: actually support mismatching treatments
        treatment_types = treatment_types[0]
        treatment_value_set = {}

        for treatment in treatment_types:
            treatment_value_set[treatment] = set(map(lambda x: x[treatment], df['well_info']['treatments']))    
            treatment_values = df['well_info']['treatments'].apply(lambda x: x[treatment])
            #treatment_values.name = ['well_info', treatment]
            treatment_values = pd.DataFrame(treatment_values)
            treatment_values.columns = pd.MultiIndex.from_product(
                    [['well_info'], [treatment]])
            df = df.join(treatment_values)

        if blank_value:
            return df.curve - blank_value

        blank_df = df[df['well_info']['blank']]
        # return(blank_df)
        blank_ODs = {}
        for treatment in treatment_value_set.keys():
            treatment_blanks = dict()
            for value in treatment_value_set[treatment]:
                treatment_blank = blank_df[blank_df['well_info']['treatments'].apply(lambda x: x[treatment])==value]['curve']
                #treatment_blank = pd.DataFrame([*treatment_blank])
                treatment_blank = treatment_blank.agg('mean')
                treatment_blank.name = value
                treatment_blanks[value] = treatment_blank
            blank_ODs[treatment] = treatment_blanks
        #blank_ODs = pd.DataFrame.from_dict(blank_ODs)
        def blank_finder_wrapper(treatment):
            def blank_finder(row):
                treatment_value = row['well_info'][treatment]
                blank_value = blank_ODs[treatment][treatment_value]
                return blank_value
            return blank_finder

        #return(blank_ODs)
        well_blank_values = df.apply(blank_finder_wrapper(treatment_types[0]),\
                                                          axis = 'columns')
        well_blank_values.columns = pd.MultiIndex.from_product([['media_blank'], well_blank_values.columns])
        #for treatment in treatment_types[1:]:
        #    well_blank_value = df.apply(blank_finder_wrapper(treatment),\
        #                                axis = 'columns')
        #    well_blank_value.name = treatment
        #    well_blank_values.join(well_blank_values)
        
        df = df.join(well_blank_values)
        def set_blank_value(well):
            self[well.name].blank_curve = well.curve-well.media_blank
        
        df.T.apply(set_blank_value)
        return blank_ODs, df, well_blank_values

    def make_sns_heatmap(self, characteristic,
                         cmap=sns.light_palette("green"),
                         ret=False,
                         context='paper',
                         dim=[list('ABCDEFGH'), list(range(1, 13))],
                         replicates=False):
        '''
        creates a seaborn heatmap of any well characteristic (specified via
        function), also accepts a specified colormap (cmap), can return the
        figure (ret), accepts context which changes relative sizes of elements,
        and can restrict the heatmap to certain wells
        '''

        data = tools.as_matrix(self.well_list,
                                            characteristic,
                                            *dim)
        with sns.plotting_context(context):
            heatmap = sns.heatmap(data,
                                  cmap=cmap,
                                  annot=True,
                                  linewidths=.5,
                                  fmt='.2f')
            heatmap.set_xticklabels(dim[1])
            heatmap.set_yticklabels(dim[0], rotation=0)
            fig = heatmap.get_figure()
            fig.set_size_inches(13, 8)
            fig.savefig("output.svg")
            if ret:
                return heatmap


'''
    def make_plotly_heatmap(self, characteristic, cmap=False, ret=False):
        from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
        import plotly.plotly as py
        import plotly.graph_objs as go
        #py.init_notebook_mode(connected=True)
        data = plate_reader_tools.as_matrix(self.well_list,
                                            characteristic,
                                            *dim)
        values = []
        for row in data:
            this_row = []
            for well in row:
                if np.isnan(characteristic(well)):
                    this_row.append(0)
                    #print(True)
                else:
                    this_row.append(characteristic(well))
                    #print(False)
            values.append(this_row)
        heatmap = go.Heatmap(z = values,
                             y = list('ABCDEFGH'),
                             x = list(range(1,13)),
                             autorange='reversed')
        data = [heatmap]
        plot(data, filename='labelled-heatmap')
        
        if ret:
            return heatmap
'''


            

def read_curve(curve):
    #check_processing(curve)
    spline, derivative = fit_spline(curve, der=[0, 1])
    phases = find_phases(derivative)
    
    lag_phase_length = phases['lag'][1] - phases['lag'][0]
    mu_max = find_max_mu(curve)
    saturation_time = phases['saturation'][0]
    
    #print('lag_phase_length', lag_phase_length, '\n'+\
    #      'mu_max', mu_max, '\n'+\
    #      'saturation_time', saturation_time)
    
    return lag_phase_length, mu_max, saturation_time

def check_processing(curve):
    #necessary?
    pass


def find_phases(derivative):
    start = np.argmin(derivative.iloc[0:int((len(derivative)/5))])

    lag_phase = [derivative.index[0]]
    lag_phase.append(start+derivative.index[np.argwhere(derivative[start:]>0.002)[0][0]])
    
    growth_phase = [lag_phase[1]]
    growth_phase.append(growth_phase[0]+derivative.index[np.argwhere(derivative[growth_phase[0]:]<0.002)[0][0]])

    stagnation_phase = [growth_phase[1]]
    stagnation_phase.append(stagnation_phase[0]+derivative.index[np.argwhere(derivative[stagnation_phase[0]:]<0.0004)[0][0]])

    saturation_phase = [stagnation_phase[1]]
    saturation_phase.append(derivative.index[-1])

    phases = {'lag': lag_phase,
              'growth': growth_phase,
              'stagnation': stagnation_phase,
              'saturation': saturation_phase}
    return phases


def color_phases(spline, phases):
    if not phases:
        spline, derivative = fit_spline(spline, der=[0, 1])
        phases = find_phases(derivative)
    fig, ax = plt.subplots()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('OD595')
    colors = ['y', 'g', 'b', 'k']
    print(spline.index)
    for i, phase in enumerate(phases.keys()):
        index = spline[phases[phase][0]:phases[phase][1]].index
        data = spline[phases[phase][0]:phases[phase][1]].values
        ax.plot(index, data, color=colors[i])

    plt.show()


def fit_spline(curve, log=False,
               visualize=False, save=False,
               filename=None, der=0, mu=False):
    '''
    Creates a spline representation of a growth curve, can also be used to
    produce the integral or any derivative of the growth curve itself
    '''
    if mu == False:
        try:
            mu = 'mu' in der
        except(TypeError):
            mu = der == 'mu'
        if mu:
            der.remove('mu')
    if save and not filename:
        exit
    if save and not visualize:
        visualize = True
    if type(der) != list:
        der = [der]
    try:
        index = curve.index/np.timedelta64(1, 'h')
    except(TypeError):
        try:
            index = list(map(lambda x: x/pd.Timedelta(hours=1), curve.index))
        except(TypeError):
            index = curve.index

    if log:
        curve = curve.apply(lambda x: np.log(x))

    # creating t c and k values of the spline, s value represents degree of
    # smoothing.
    tck = interpolate.splrep(index, curve.values, s=1.0*10**-3)
    splines = []
    for derivative in der:
        if derivative >= 0:
            spline = interpolate.splev(index, tck, der=derivative)
            spline = pd.Series(spline, index)
            splines.append(spline)
        if derivative < 0:
            Bspline = (interpolate.BSpline(*tck, extrapolate=False))
            integral = interpolate.splantider(Bspline, -derivative)
            spline = (interpolate.splev(index, integral, der=0))
            spline = pd.Series(spline, index)
            splines.append(spline)
    if mu:
        mu = interpolate.splev(index, tck, 1) / \
                       interpolate.splev(index, tck, 0)
        mu = pd.Series(mu, index)
        splines.append(mu)
        print('mu')

    if visualize:
        fig, axis1 = plt.subplots()
        axis1.set_xlabel('Time (hours)')
        axis1.set_ylabel('OD595', color='r')
        if 0 in der:
            axis1.plot(index, splines[0], '-',
                       index, curve.values, '.', color='r')
        else:
            axis1.plot(index, splines[0], '-', color='r')

        axis1.tick_params(axis='y', labelcolor='r')

        secondary_axes = {}
        for i, spline in enumerate(splines[1:]):
            secondary_axes['axis'+str(i)] = axis1.twinx()
            secondary_axes['axis'+str(i)].axhline(color='b')
            secondary_axes['axis'+str(i)].set_ylabel('der ='+str(i), color='b')
            secondary_axes['axis'+str(i)].plot(index, spline, color='b')
            secondary_axes['axis'+str(i)].tick_params(axis='y', labelcolor='b')

        fig.tight_layout()

        plt.show()

        if save:
            directory = r'C:\Users\Owner\Concordia\Lab_Automation\output_files\\'
            fig.savefig(directory+filename)

    return splines
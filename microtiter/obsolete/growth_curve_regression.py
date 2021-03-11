# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 14:58:05 2018

@author: James
"""
from scipy.optimize import curve_fit as curve_fit
import matplotlib.pyplot as plt


def growth_curve(time, rate, carrying_capacity,
                 initial_population, lag_phase, c):
    e = 2.71828
    return carrying_capacity/(1+initial_population*e**(-rate*time)-lag_phase)+c


def get_growth_char(index, curve):
    # these bounds aren't really based on anything in particular, just trial
    # and error
    bounds = ((0, 0., 0., -5, -10), (.2, 3., 2., 5, 6000))
    characters, covariances = curve_fit(growth_curve, index, curve,
                                        bounds=bounds)
    rate, carrying_capacity, initial_population, lag_phase, c = characters
    return rate, carrying_capacity, initial_population, lag_phase, c


def logistic_regression_graph(index, curve):
    popt = tuple(get_growth_char(index, curve))
    actual_results = (index, curve)
    regression_plot = (index, growth_curve(index, *popt))
    plt.plot(*actual_results, 'r', *regression_plot, 'b')

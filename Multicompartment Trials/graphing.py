# -*- coding: utf-8 -*-
"""

Graphing class design to make graph objects.
Input can be dataframe or arrays.

- set colours and themes for ions and voltages
- set style of no major gridlines etc.
- integration with bqplot for user interface


@author: E Shorer
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd



class graph(object):

    def __init__(self, Title = ''):
        self.title = Title
        self.fig, self.axes = plt.subplots()
        self.x_arr = []
        self.y_arr = []
        self.graph_form = {}

    def final_vals(self, dataframe):
        self.df = dataframe

    def set_x_axis(self, x_arr, x_title="x axis"):
        self.x_arr = x_arr
        self.x_title = x_title

    def set_y_axis(self, y_arr, y_title="y axis"):
        self.y_arr = y_arr
        self.y_title = y_title

    def plot_it(self):
        sns.despine()
        plt.plot(self.x_arr, self.y_arr)

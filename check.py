# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 12:46:40 2018

@author: pnolan86
"""

from netCDF4 import Dataset
root = Dataset('windagedata.nc','r')
vars = root.variables
t = vars['time'][:]
root.close()
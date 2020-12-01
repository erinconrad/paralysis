#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 08:47:06 2020

@author: erinconrad
"""

# Parameters
xls_path = "/Users/erinconrad/Desktop/residency stuff/R25/paralysis/data/test.xlsx"

# import numpy, scipy, xlrd, others
import xlrd
import scipy
import scipy.stats


# Load data file
book = xlrd.open_workbook(xls_path)
sheet = book.sheet_by_index(0)
data = [[sheet.cell_value(r, c) for c in range(sheet.ncols)] for r in range(sheet.nrows)]


# Statistical testing
oddsratio, pvalue = scipy.stats.fisher_exact(data)
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 08:42:48 2023

@author: quemeb
"""

import pandas as pd
import numpy as np
import pyreadstat

# Importing dataset
data, meta = pyreadstat.read_sav(r"C:\Users\bryan\Desktop\USC_research\Withers\EXPLORING MENTAL HEALTH AND RESILIENCY.sav")

# Dropping made up dicho variables
cols_to_drop = [col for col in data.columns if "_01" in col]
data.drop(columns=cols_to_drop, inplace=True)
data = data[data.Ageyears.between(18, 75)]

# Making all variables lower case
data.columns = map(str.lower, data.columns)

# Labeling insurance
data['insu'] = data['insu'].map({1: "None", 2: "Partial", 4: "Full Private", 10: "Full Public"})

# Dropping more stuff
data = data[data.race01 != 4]
data = data[~data.religion01.isin([3, 4, 5, 6])]

# Creating MENTALHEALTH score from scratch
vars_to_replace = [f'a0{num}' for num in range(35, 48)]
for var in vars_to_replace:
    data[var].replace(2, 0, inplace=True)
data['mentalhealth'] = data[vars_to_replace].sum(axis=1)
data.dropna(subset=['mentalhealth'], inplace=True)

# Creating binary Mental Health outcome 
data['mental_cat'] = (data['mentalhealth'] >= 1).astype(int)
data.drop(columns=vars_to_replace + ['mentalhealth'], inplace=True)

# Back to data cleaning
cols_to_drop = [
    'heightincm', 'weightinkg', 'ver', 'surveyversion', 'n', 'respondentid',
    'a187', 'a188', 'a189', 'a190', 'a191', 'a174', 'a175', 'a176', 'a180', 'pleaseindicateyourworkorstudystatus01'
]
data.drop(columns=cols_to_drop, inplace=True)

# RESILIENCE
data['res'] = data[['a156', 'a157', 'a158', 'a159', 'a160', 'a161']].mean(axis=1)
data.dropna(subset=['res'], inplace=True)
data.drop(columns=['resilience', 'a156', 'a157', 'a158', 'a159', 'a160', 'a161'], inplace=True)

# Generating Ordinal Resilience 
data['res_cat'] = 0
data.loc[data['res'] >= 3, 'res_cat'] = 1
data.loc[data['res'] >= 4.31, 'res_cat'] = 2
data.drop(columns=['r3', 'res'], inplace=True)

# COVID_EFFECT
data['cov_effect'] = data[['a162', 'a163', 'a164', 'a165', 'a166', 'a167']].sum(axis=1)
data.drop(columns=['covid_effect', 'a162', 'a163', 'a164', 'a165', 'a166', 'a167'], inplace=True)

# COVID_FEAR
data['cov_fear'] = data[['a168', 'a169', 'a170', 'a171', 'a172']].sum(axis=1)
data.drop(columns=['covid_fear', 'a168', 'a169', 'a170', 'a171', 'a172'], inplace=True)

# COVID_PHYCH
data['cov_psych'] = data[['a173', 'a177', 'a178', 'a179', 'a181', 'a182', 'a183', 'a184', 'a185', 'a186']].sum(axis=1)
data.drop(columns=['covid_psych', 'a173', 'a177', 'a178', 'a179', 'a181', 'a182', 'a183', 'a184', 'a185', 'a186'], inplace=True)

# WELLBEING
data['cov_wellbeing'] = data[['a211', 'a212', 'a213', 'a214', 'a215']].sum(axis=1)
data['well_cat'] = (data['cov_wellbeing'] >= 17.5).astype(int)
data.drop(columns=['a211', 'a212', 'a213', 'a214', 'a215', 'wellbeing', 'wellbeing01', 'cov_wellbeing'], inplace=True)

# COVID_BURDEN
burden_vars = [f'a{num}' for num in range(192, 211)]
for var in burden_vars:
    data[var].replace(77, 0, inplace=True)
data['cov_burden'] = data[burden_vars].sum(axis=1)
data.drop(columns=burden_vars, inplace=True)

# COVID_CONTACT
contact_vars = ['a002', 'a003', 'a005']
for var in contact_vars:
    data[var] = np.where(data[var] == 2, 1, np.where(data[var] == 1, 2, 0))
data['cov_contact'] = data[contact_vars].sum(axis=1)
data.drop(columns=contact_vars, inplace=True)

# Define y and x based on the provided globals
y = 'mental_cat'
x_continuous = ['ageyears', 'cov_effect', 'cov_fear', 'cov_psych', 'cov_burden', 'cov_contact']
x_categorical = [col for col in data.columns if col not in x_continuous + [y, 'res_cat']]
x = x_continuous + x_categorical



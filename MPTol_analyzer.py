#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""@Author: Vivek Srinivas - Baliga Lab, ISB

This is a program to calculate drug tolerance from multiwell plate assay
Program requires two files - 'Labels.csv' and 'Readings.txt' for example
    ##refer to examples files for format## 
    Step 1: Plot growth
            plot_growth_in_wells(results_text_file, labels, instrument = "BioTek"/"BioAnalyzer")
    Step 2: Identify OD threshold for 'Start of growth point determination'
    Step 3: Calculate drug tolerance with command -
            calculate_plot_tolerance(results_text_file, labels, ODT = __,Drug = __, instrument = "BioTek"/"BioAnalyzer" )
"""

import csv, re
import pandas as pd
from pandas import DataFrame as DF
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt
import seaborn as sns
from datetime import datetime as DT
import warnings

#- - - - - - - - - - - - ########### - - - - - - - - - - - - - - #
def get_mins(i):
    h,m,s = (float(j) for j in  i.split(":"))
    return (h*60.0)+m

def brkstr(string):
    return [float(i) for i in string.split(";")]

def read_results_biotek(results_text_file):
    with open(results_text_file,"r") as f:
        lines = [ i[:-2] for i in f.readlines()]
    data1 = lines[lines.index("OD600:600"):]
    data = data1[:data1.index("Results")]
    data_split = np.array([i.split("\t") for i in data[2:-1]]).T
    Time_points = [str(get_mins(i)) for i in data_split[0][1:]]
    return data_split[2:], Time_points

def read_results_bioanalyzer(results_text_file):
    data = pd.read_csv(results_text_file,header=0).transpose()
    time = data[0,:].values
    array = [[i]+list(l.values) for i,l in data.iloc[1:].iterrows()]
    return time, array

def read_results_text_file(results_text_file, labels, instrument_type):
    labels = pd.read_csv(labels,index_col=0,header = 0)
    ret_DF = DF()
    if instrument_type == "BioTek":
        data_split, Time_points = read_results_biotek(results_text_file)
    elif instrument_type=="BioAnalyzer":
        data_split, Time_points = read_results_bioanalyzer(results_text_file)
    for n,zz in enumerate(data_split):
        if zz[0] in labels.index:
            ret_DF.at[n,"Well"] = zz[0]
            ret_DF.at[n,"Drug"] = labels.loc[zz[0]].Drug
            ret_DF.at[n,"Concentration"] = labels.loc[zz[0]].Concentration
            ret_DF.at[n,"Incubation_time"] = labels.loc[zz[0]].Incubation_time
            ret_DF.at[n,"Sample_type"] = labels.loc[zz[0]].Type
            ret_DF.at[n,"Time_points"]= ";".join(Time_points)
            ret_DF.at[n,"OD"] = ";".join(zz[1:])
        else:
            print "%s well is not labelled correctly and will be skipped"%zz[0]
    return ret_DF

def plot_growth_in_wells(results_text_file,labels,**kwargs):
    """Requires 'results text file' and 'labels csv file'.
    -- Use this to determine OD threshold""" 
    data = read_results_text_file(results_text_file,labels, kwargs["instrument"])
    for n,i in data.iterrows():
        if i.Sample_type == "Control":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "red", linewidth=1, alpha = 1, label = "Controls - Well %s"%i.Well)
        else:
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "grey", linewidth=0.5, alpha = 0.5)
    plt.xlabel("Time (mins)",fontsize=16)
    plt.ylabel("OD600",fontsize=16)
    plt.legend(loc = "upper left")
    plt.show()

def get_sgt(od,time, odt):
    warnings.simplefilter('ignore', np.RankWarning)
    model = np.polyfit(od,time,2)
    det_model = np.poly1d(model)
    sgt = det_model(odt)
    if sgt > min(time) and sgt < max(time)+np.mean(time):
        ret_sgt = sgt
    else:
        ret_sgt = np.nan
    return ret_sgt

def add_sgt(results_text_file, labels, odt, instrument):
    data = read_results_text_file(results_text_file,labels, instrument)
    for n,i in data.iterrows():
        data.at[n,"SGT"] = get_sgt(brkstr(i.OD),brkstr(i.Time_points),odt)
    return data

def plot_sgt_in_wells(results_text_file,labels,ODT, **kwargs):
    """Requires 'results text/csv file' and 'labels csv file'.
    -- Use this to determine OD threshold""" 
    data = add_sgt(results_text_file,labels,ODT, kwargs["instrument"])
    for n,i in data.iterrows():
        if i.Sample_type == "Control":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "red", linewidth=1, alpha = 1, label = "Controls - Well %s"%i.Well)
            plt.plot(i.SGT,ODT,"o",color = "red",markersize=10)
        else:
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "grey", linewidth=0.5, alpha = 0.5)
            plt.plot(i.SGT,ODT,"o",color = "grey",markersize=7,alpha=0.7)
    plt.xlabel("Time (mins)",fontsize=16)
    plt.ylabel("OD600",fontsize=16)
    plt.legend(loc = "upper left")
    plt.show()



def calculate_tolerance(results_text_file,labels,ODT, instrument):
    data = add_sgt(results_text_file,labels,ODT, instrument)
    data = data.drop(["Time_points","OD"], axis =1)
    control = np.nanmean(data[data.Sample_type.isin(["Control"])].SGT.values)
    if control < 0:
        control = 10**5
    else:
        control = control
    for n,i in data.iterrows():
        data.at[n,"log2_surviving_cells"] = np.log2(2**(-(i.SGT-control)))
    return data
        
def calculate_plot_tolerance(results_text_file,labels,**kwargs):
    """Requires 'results text/csv file', 'labels csv file', ODT[optional], Drug and instrument type.
    -- Use this to calculate and plot % cell survival"""
    if "ODT" not in kwargs.keys():
        kwargs["ODT"] = 0.15
    else:
        pass
    if "drug" not in kwargs.keys():
        print "Provide drug name"
    else:
        data = calculate_tolerance(results_text_file,labels,kwargs["ODT"], kwargs["instrument"])
        if kwargs["drug"] not in data.Drug.values:
            print "Enter correct name of the drug"
        elif kwargs["instrument"] not in data.Drug.values:
            print "Enter instrument used for reading OD ('BioTek'/'BioAnalyzer')"
        else:
            data = data[data.Drug.isin([kwargs["drug"]])]
            tolerance = data.pivot(index="Incubation_time",columns="Concentration",values="log2_surviving_cells")
            sns.heatmap(tolerance, cmap="RdYlGn", linewidths=0.30,cbar_kws={'label': 'log2 (Fraction of surviving cells)'})
            plt.title("%s treatment profile"%kwargs["drug"])
            plt.savefig("%s_%s.pdf"%(results_text_file.split(".")[0],kwargs['drug']))
            plt.clf()
            data.to_csv("%s_%s_analyzed.csv"%(results_text_file.split(".")[0],kwargs['drug']))

    

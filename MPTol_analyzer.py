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
            calculate_plot_tolerance(results_text_file, labels, ODT = __, seperateby = "drug"/"type", drug/type = __, instrument = "BioTek"/"BioAnalyzer" )
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

def brkstr2(string):
    return string.split

def convert_to_strings(listx):
	new_list = []
	for i in listx:
		new_list.append(str(i))
	return new_list
   

def read_results_biotek(results_text_file):
    with open(results_text_file,"r") as f:
        lines = [ i[:-2] for i in f.readlines()]
    p1, p2 = re.compile("Time"), re.compile("Result")
    idp1, idp2 = list(filter(p1.match,lines)),list(filter(p2.match,lines))    
    idx1, idx2 = lines.index(idp1[-1]), lines.index(idp2[-1])
    data = lines[idx1:idx2-1]
    data_split = np.array([i.split("\t") for i in data]).T
    Time_points = [str(get_mins(i)) for i in data_split[0][1:]]
    od_readings = data_split[2:]
    ork = len([i for i in od_readings[0] if i != ""])
    return od_readings[:,:ork], Time_points[:ork-1]

def read_results_bioanalyzer(results_text_file):
    data = pd.read_csv(results_text_file,header=0).transpose()
    time = data.iloc[0,:].values
    Time_points = [str(get_mins(i)) for i in time]
    od_readings = [convert_to_strings([i]+list(l.values)) for i,l in data.iloc[1:,:].iterrows()]
    return od_readings, Time_points

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
    if sgt > min(time) and sgt < 2.0*(max(time)):#:+np.mean(time):
        ret_sgt = sgt
    elif sgt > 2.0*(max(time)):#:+np.mean(time):
        ret_sgt = 2.0*(max(time))#:+np.mean(time)
    else:
        ret_sgt = 0.5
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
    control_sgts = data[data.Sample_type.isin(["Control"])].SGT.values
    control = np.nanmean(control_sgts)
    if control > 0: 
        control = control
    else:
        control = max(control_sgts)
    print "control SGT = %s"%control
    for n,i in data.iterrows():
        if i.SGT < control:
            data.at[n,"log2_surviving_cells"] = np.log2(control-i.SGT)
        else:
            data.at[n,"log2_surviving_cells"] = 0.0
    return data
        
def calculate_and_plot_tolerance(results_text_file,labels,**kwargs):
    """Requires 'results text/csv file', 'labels csv file', ODT[optional], Drug and instrument type.
    -- Use this to calculate and plot % cell survival"""
    if "ODT" not in kwargs.keys():
        kwargs["ODT"] = 0.15
    else:
        pass
    if "seperateby" not in kwargs.keys():
        print "Provide an attribute to segregate data"
    else:
        if kwargs["seperateby"] in ["drugs", "drug","Drug","Drugs"]:
            if "drug" not in kwargs.keys():
                print "Provide drug name"
            else:
                data = calculate_tolerance(results_text_file,labels,kwargs["ODT"], kwargs["instrument"])
                if kwargs["drug"] not in data.Drug.values:
                    print "Enter correct name of the drug"
                else:
                    data = data[data.Drug.isin([kwargs["drug"]])]
                    tolerance = data.pivot(index="Incubation_time",columns="Concentration",values="log2_surviving_cells")
                    sns.heatmap(tolerance, cmap="RdBu", vmin = 0,vmax =12, linewidths=0.30,cbar_kws={'label': 'log2 (Fraction of surviving cells)'})
                    plt.title("%s treatment profile"%kwargs["drug"])
                    plt.savefig("%s_%s.pdf"%(results_text_file.split(".")[0],kwargs['drug']))
                    plt.clf()
                    data.to_csv("%s_%s_analyzed.csv"%(results_text_file.split(".")[0],kwargs['drug']))
        elif kwargs["seperateby"] in ["sample", "type","sampletype","sample_type"]:
            if "type" not in kwargs.keys():
                print "Provide the sample type"
            else:
                data = calculate_tolerance(results_text_file,labels,kwargs["ODT"], kwargs["instrument"])
                if kwargs["type"] not in data.Sample_type.values:
                    print "Enter the correct sample type"
                else:
                    data = data[data.Sample_type.isin([kwargs["type"]])]
                    tolerance = data.pivot(index="Incubation_time",columns="Concentration",values="log2_surviving_cells")
                    sns.heatmap(tolerance, cmap="RdBu",vmin=0,vmax=12, linewidths=0.30,cbar_kws={'label': 'log2 (Fraction of surviving cells)'})
                    plt.title("%s treatment profile of %s"%(data.Drug.values[0],kwargs['type']))
                    plt.savefig("%s_%s.pdf"%(results_text_file.split(".")[0],kwargs['type']))
                    plt.clf()
                    data.to_csv("%s_%s_analyzed.csv"%(results_text_file.split(".")[0],kwargs['type']))
                

    

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""@Author: Vivek Srinivas - Baliga Lab, ISB

This is a program to calculate drug tolerance from multiwell plate assay
Program requires two files - 'Labels.csv' and 'Readings.txt' for example
    ##refer to examples files for format## 
    Step 1: Plot growth
            plot_growth_in_wells(results_text_file, labels, instrument = "BioTek1"/"BioTek2"/"BioAnalyzer")
    Step 2: Identify OD threshold for 'Start of growth point determination'
    Step 3: Calculate drug tolerance with command -
            calculate_plot_tolerance(results_text_file, labels, ODT = __, instrument = BioTek1"/"BioTek2"/"BioAnalyzer")
            ----OR----
            collate_calculate_plot_estimated_survival_for_strain_type(results_text_file, labels, ODT = __, instrument = BioTek1"/"BioTek2"/"BioAnalyzer",Cohort_name = "__",experiment = ["__"])
"""

import csv, re,os
import pandas as pd
from pandas import DataFrame as DF
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt
import seaborn as sns
from datetime import datetime as DT
import warnings, math
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm

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
    
def blank_zero_hour(odl):
    ret_odl = []
    for i in odl:
        if i <= odl[0]:
            ret_odl.append(str(0.0))
        else:
            ret_odl.append(str(float(i)-float(odl[0])))
    return ret_odl

def read_results_biotek1(results_text_file):
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

def read_results_biotek2(results_text_file):
    data = pd.read_csv(results_text_file,skiprows = range(0,12), header=0).transpose()
    time = data.iloc[0,:].values
    Time_points = [str(get_mins(i)) for i in time]
    od_readings = [convert_to_strings([i.strip()]+list(l.values)) for i,l in data.iloc[2:,:].iterrows()]
    return od_readings, Time_points

def read_results_bioanalyzer(results_text_file):
    data = pd.read_csv(results_text_file,header=0).transpose()
    time = data.iloc[0,:].values
    Time_points = [str(get_mins(i)) for i in time]
    od_readings = [convert_to_strings([i]+list(l.values)) for i,l in data.iloc[1:,:].iterrows()]
    return od_readings, Time_points

def read_results_text_file(results_text_file, labels, instrument_type,**kwargs):
    labels = pd.read_csv(labels,index_col=0,header = 0)
    ret_DF = DF()
    if instrument_type == "BioTek1":
        data_split, Time_points = read_results_biotek1(results_text_file)
    elif instrument_type == "BioTek2":
        data_split, Time_points = read_results_biotek2(results_text_file)
    elif instrument_type=="BioAnalyzer":
        data_split, Time_points = read_results_bioanalyzer(results_text_file)
    for n,zz in enumerate(data_split):
        if zz[0] in labels.index:
            ret_DF.at[n,"Well"] = zz[0]
            ret_DF.at[n,"Drug"] = labels.loc[zz[0]].Drug
            ret_DF.at[n,"Concentration"] = labels.loc[zz[0]].Concentration
            ret_DF.at[n,"Incubation_time"] = labels.loc[zz[0]].Incubation_time
            ret_DF.at[n,"Sample_type"] = labels.loc[zz[0]].Type
            ret_DF.at[n,"Experiment"] = labels.loc[zz[0]].Experiment
            ret_DF.at[n,"Time_points"]= ";".join(Time_points) 
            ret_DF.at[n,"OD"] = ";".join(blank_zero_hour(zz[1:]))
        else:
            print "%s well is not labelled correctly and will be skipped"%zz[0]
    if "seg" in kwargs.keys():
        if kwargs["seg"] is True:
            ret_DF = ret_DF.loc[ret_DF["Experiment"].isin([kwargs["exp"]])]
        else:
            pass
    else:
        pass
    return ret_DF

def plot_growth_in_wells(results_text_file,labels,**kwargs):
    """Requires 'results text file' and 'labels csv file'.
    -- Use this to determine OD threshold"""
    if "seg" in kwargs.keys():
        data = read_results_text_file(results_text_file,labels, kwargs["instrument"],seg = kwargs["seg"],exp = kwargs["Experiment"])
        print data.shape
    else:
        data = read_results_text_file(results_text_file,labels, kwargs["instrument"])
    for n,i in data.iterrows():
        if i.Sample_type == "ND-control":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "blue", linewidth=1, alpha = 1)
        elif i.Sample_type == "MDK-control":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "red", linewidth=1, alpha = 1)
        elif i.Sample_type == "Sample":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "grey", linewidth=0.5, alpha = 0.5)
        else:
            pass
    plt.xlabel("Time (mins)",fontsize=16)
    plt.ylabel("OD600",fontsize=16)
    #plt.legend(loc = "upper left")
    plt.show()

def get_sgt(od,time, odt):
    warnings.simplefilter('ignore', np.RankWarning)
    try:
        model = np.polyfit(od,time,3)
        det_model = np.poly1d(model)
        sgt = det_model(odt)
        if sgt < min(time):
            ret_sgt = min(time)
        elif sgt > 2.0*(max(time)):
            ret_sgt = 2.0*(max(time))
        else:
            ret_sgt = sgt
    except:
        ret_sgt = 2.0*(max(time))
    return ret_sgt

def add_sgt(results_text_file, labels, odt, instrument,**kwargs):
    if "seg" in kwargs.keys():
        data = read_results_text_file(results_text_file,labels, instrument,seg = kwargs["seg"],exp = kwargs["Experiment"])
        print data.shape
    else:
        data = read_results_text_file(results_text_file,labels, instrument)
    for n,i in data.iterrows():
        data.at[n,"SGT"] = get_sgt(brkstr(i.OD),brkstr(i.Time_points),odt)
    return data

def plot_sgt_in_wells(results_text_file,labels,ODT, **kwargs):
    """Requires 'results text/csv file' and 'labels csv file'.
    -- Use this to determine OD threshold"""
    if "seg" in kwargs.keys():
        data = add_sgt(results_text_file,labels,ODT, kwargs["instrument"],seg = kwargs["seg"],exp = kwargs["Experiment"])
    else:
        data = add_sgt(results_text_file,labels,ODT, kwargs["instrument"])
    for n,i in data.iterrows():
        if i.Sample_type == "ND-control":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "blue", linewidth=1, alpha = 1)
        elif i.Sample_type == "MDK-control":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "red", linewidth=1, alpha = 1)
        elif i.Sample_type == "Sample":
            plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--",color = "grey", linewidth=0.5, alpha = 0.5)
        else:
            pass
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


def get_control_values(data):
    ret_dict = {}
    time_points = brkstr(data.Time_points.values[0])
    for i in list(set(data.Experiment.values)):
        sub_data = data.loc[data["Experiment"].isin([i])]
        sgtvals = sub_data.SGT.values
        if "ND-control" in sub_data.Sample_type.values:
            mean_ND_control = np.nanmean(sub_data[sub_data.Sample_type.isin(["ND-control"])].SGT.values)
        else:
            print "ND-control is absent in %s - experiment, value of fastest growing culture among the experimental set will be used (This may result in inaccurate results)"%i
            mean_ND_control = min(sgtvals)
        if "MDK-control" in sub_data.Sample_type.values:
            mean_MDK_control = np.nanmean(sub_data[sub_data.Sample_type.isin(["MDK-control"])].SGT.values)
        else:
            print "MDK-control is absent in %s - experiment, average value of 48 hours will be used (This may result in inaccurate results)"%i
            mean_ND_control = 2880.0
        if mean_ND_control > 0 and mean_ND_control < max(time_points)\
           and mean_MDK_control > 0: 
            pass
        else:
            mean_ND_control = min(sgtvals)
            mean_MDK_control = 2880.0
        ret_dict[i] = {"mean_ND_control":mean_ND_control,"mean_MDK_control":mean_MDK_control}
        print "For experiment - %s : Mean SGTs of ND-control = %0.3f and MDK-control = %0.3f"%(i,mean_ND_control, mean_MDK_control)
    return ret_dict

def normalize_to_t0(data):
    t0dict = {}
    for n1,rows1 in data.iterrows():
        if float(rows1.Incubation_time) == 0.0:
            if rows1.Fraction_of_surviving_cells > 0.0:
                t0dict["%s-%s"%(rows1.Experiment,rows1.Concentration)] = rows1.Fraction_of_surviving_cells
            else:
                t0dict["%s-%s"%(rows1.Experiment,rows1.Concentration)] = 1.0
        else:
            pass
    for n2,rows2 in data.iterrows():
        if rows2.Fraction_of_surviving_cells > 0.0:
            ss = "%s-%s"%(rows2.Experiment,rows2.Concentration)
            data.at[n2,"Normed_fraction_of_surviving_cells"] = rows2.Fraction_of_surviving_cells/t0dict[ss]
        else:
            data.at[n2,"Normed_fraction_of_surviving_cells"] = 0.0
    return data

def calculate_tolerance(results_text_file,labels,ODT, instrument):
    data = add_sgt(results_text_file,labels,ODT, instrument)
    control_vals = get_control_values(data)
    time_points = brkstr(data.Time_points.values[0])
    data = data.drop(["Time_points","OD"], axis =1)
    data = data[data.Sample_type.isin(["Sample"])]
    for n,i in data.iterrows():
        exp = i.Experiment
        if i.SGT < control_vals[exp]["mean_MDK_control"]:
            if i.SGT < min(time_points):
                SGT_for_calc = min(time_points)
            else:
                SGT_for_calc = i.SGT
            data.at[n,"Fraction_of_surviving_cells"] = 1.01-(2.0**((control_vals[exp]["mean_MDK_control"]-SGT_for_calc)/-control_vals[exp]["mean_ND_control"]))
        else:
            data.at[n,"Fraction_of_surviving_cells"] = 0.0
    return normalize_to_t0(data)
        
       
def calculate_and_plot_tolerance(results_text_file,labels,**kwargs):
    savepath = results_text_file.split("/")[0]+"/"
    if "Results" in os.listdir(savepath):
        pass
    else:
        os.mkdir(savepath+"Results")
    """Requires 'results text/csv file', 'labels csv file', ODT[optional] and instrument type.
    -- Use this to calculate and plot % cell survival"""
    if "ODT" not in kwargs.keys():
        kwargs["ODT"] = 0.15
    else:
        pass
    data = calculate_tolerance(results_text_file,labels,kwargs["ODT"],kwargs["instrument"])
    for experiment in list(set(data.Experiment.values)):
        sub_data = data[data.Experiment.isin([experiment])]
        tolerance = sub_data.pivot(index="Incubation_time",columns="Concentration",values="Normed_fraction_of_surviving_cells")
        ax = sns.heatmap(tolerance,cmap = "RdBu",vmin=0.01,center = 1.0, vmax=1.0001,linewidths=0.30,cbar_kws={'label': 'Fraction of surviving cells'})
        ax.invert_yaxis()
        plt.title("%s treatment profile in %s"%(sub_data.Drug.values[0],experiment))
        plt.savefig("%s/%s.pdf"%(savepath+"Results",experiment))
        plt.clf()
        sub_data.to_csv("%s/%s_analyzed.csv"%(savepath+"Results",experiment))

def kill_curve_surface(param_list,k):
    Conc, inc_time = param_list
    return 1-(((k*Conc)/(1+Conc))*inc_time)

def calculate_expected_asymptote(data):
    Conc = np.array([rows.Concentration for n, rows in data.iterrows()])
    inc_time = np.array([rows.Incubation_time for n, rows in data.iterrows()])
    ydata = np.array([rows.Normed_fraction_of_surviving_cells for n, rows in data.iterrows()])
    popt,pcov = curve_fit(kill_curve_surface,(Conc,inc_time),ydata)
    estimated_asymptote = pd.DataFrame()
    print popt
    for c in np.arange(data.Concentration.min(),1000.0,50):#,data.Concentration.max(),1.0):
        for t in np.arange(data.Incubation_time.min(),94.0,8.0):#data.Incubation_time.max(),1.0):
            n = "%s-%s"%(c,t)
            estimated_asymptote.at[n,"xMIC"]=c
            estimated_asymptote.at[n,"Time"]=t
            estimated_asymptote.at[n,"Survival"]=kill_curve_surface([c,t],*popt)
    return estimated_asymptote
        
    




def collate_calculate_plot_estimated_survival_for_strain_type(results_text_file,labels,**kwargs):
    savepath = results_text_file.split("/")[0]+"/"
    if "Results" in os.listdir(savepath):
        pass
    else:
        os.mkdir(savepath+"Results")
    """Requires 'results text/csv file', 'labels csv file', ODT[optional] and instrument type.
    -- Use this to calculate and plot % cell survival"""
    if "ODT" not in kwargs.keys():
        kwargs["ODT"] = 0.15
    else:
        pass
    data = calculate_tolerance(results_text_file,labels,kwargs["ODT"],kwargs["instrument"])
    sub_data = data[data.Experiment.isin(kwargs["experiment"])]
    exp_tol = calculate_expected_asymptote(sub_data)
    tolerance = exp_tol.pivot(index="Time",columns="xMIC",values="Survival")
    ax = sns.heatmap(tolerance,cmap = "Greens",vmin=0.01, vmax=1,linewidths=0.30,\
                     cbar_kws={'label': 'Estimated fraction of cells surviving'})#,\
                     #norm=LogNorm(vmin=exp_tol.Survival.min(), vmax=exp_tol.Survival.max()))
    ax.invert_yaxis()
    plt.title("Probable kill curve for %s treatment profile in %s"%(sub_data.Drug.values[0],kwargs["Cohort_name"]))
    plt.savefig("%s/%s.pdf"%(savepath+"Results",kwargs["Cohort_name"]))
    plt.clf()
    sub_data.to_csv("%s/%s.csv"%(savepath+"Results",kwargs["Cohort_name"]))



def plot_growth_in_wells_of_experiment(results_text_file,labels,**kwargs):
    """Requires 'results text file' and 'labels csv file'.
    -- Use this to determine OD threshold"""
    data = read_results_text_file(results_text_file,labels, kwargs["instrument"])
    for experiment in list(set(data.Experiment.values)):
        if experiment in kwargs["explist"]:
            sub_data = data[data.Experiment.isin([experiment])]
            for n,i in sub_data.iterrows():
                if i.Concentration == kwargs["conc"] and i.Incubation_time == kwargs["inctime"] and i.Sample_type not in ["MDK-control","No cell"]:
                    plt.plot(brkstr(i.Time_points),brkstr(i.OD),"--", linewidth=3, alpha = 1,label = experiment,color = kwargs["col_dict"][experiment])
    plt.xlabel("Time (mins)",fontsize=16)
    plt.ylabel("OD600",fontsize=16)
    #plt.legend(loc = "upper left")
    #plt.show()

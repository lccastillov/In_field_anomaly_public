#####################################################################

## ---------------------------
##
## Script name: Match_satellite_plot_info
##
## Purpose of script:
#                   1. Match the field date per plot with the outputs of a
##                  histogram analysis that classifies pixels into anomalous and non-anomalous.
##                  The stats per field are stored into a hdf5 file
##                 2. Fit a linear regression to predict yield using only statistics derived from different
##                 vegetation indices rasters

##
## Author: Liliana Castillo Villamor. Affiliation: Aberystwyth University

##
## Date Created: 01/03/2021
##
## Copyright (c) Liliana Castillo Villamor, 2021
## Email: lic42@aber.ac.uk
##
## ---------------------------
##
## Notes:
## the hdf5 file is created by using the script in_field_anomaly_detection.py
##
## ---------------------------

#Import required packages

from Functions_Multitemp import Retrieve_data
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

## set working directory when using Docker
m_root = '/data'

#Path to the database with field data
Database=m_root+'/Ibague/Ground_Data_Ibague/Database/IbagueBD.db'

#path to the hdf5 file with the info f anomalies per plot
hdf5_field_anomalies=m_root+'/Ibague/Anomalies_interpretation/hdf5_fieldanomalies.hdf5'

#Folder that contains the anomaly rasters (for different VIs)
anomalies_folder=m_root+'/Ibague/Anomalies'

#csv files with te  list of all the images available
list_s2_csv=m_root+'/Ibague/Anomalies/S2bands.csv'

## List of all the VIs to be analised
s2_vi_list=['NDVI8A','NDVI8','EVI8A','EVI8', 'GCI8A', 'GCI8','GNDVI8A', 'GNDVI8','RECI8A_RE_B5', 'RECI8A_RE_B6',  'RECI8A_RE_B7', 'RECI8_RE_B5', 'RECI8_RE_B6',  'RECI8_RE_B7', \
            'RENDVI8A_RE_B5', 'RENDVI8A_RE_B6','RENDVI8A_RE_B7', 'RENDVI8_RE_B5', 'RENDVI8_RE_B6','RENDVI8_RE_B7', 'SAVI8A','SAVI8',\
            'NDWI8_SWIR1', 'NDWI8_SWIR2', 'NDWI8A_SWIR1', 'NDWI8A_SWIR2']



s2_dates_np=pd.to_datetime(s2_dates_np, format='%Y%m%d')

#file with the age of the plots at each of the available scenes
df_age, = Retrieve_data.plots_age(Database, s2_dates_np)
for vi in s2_vi_list:
    prefix ='SEN2'

    #imports hdf5 file with info onf in-field anomalies per plot
    # and matches them with the following field data:
    # Emergence and harvest dates, variety, plant density.
    df_std, df=Retrieve_data.import_hdf5_anomalies(m_root, vi, prefix, df_age, hdf5_field_anomalies, anomalies_folder)
    df, = Retrieve_data.conditions_stages(df)
    #Convert Yield column (Object) into float type
    df["Yield"] = df.Yield.astype(float)

    #print yield histogram

    """
    sns_plot=sns.boxplot(x='stage', y='mean_vi', hue='anomalous', data=df)
    plt.xticks(rotation=90)
    plt.show()
    plt.close()
    """

    # Remove "null" values
    df = df[df['Yield'].notnull()]
    df = df[df['Variety'].notnull()]
    df = df[df['stage'].notnull()]
    df = df[df['Emergence month'].notnull()]

    #trim blank spaces
    df.Variety = df.Variety.str.strip()


    df=df[(df.mean_vi >=-1)]
    print("List unique varieties", df.Variety.unique())
    #Run linear regression to relate plot satege, emerence month, variety an mean VI
    #and prints summary

    mod3 = smf.ols (formula = 'Yield ~ stage+Emergence_month +mean_vi + Variety', data = df).fit()
    print("mod ",vi,mod3.summary())




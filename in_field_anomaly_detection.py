## ---------------------------
##
## Script name: Match_satellite_plot_info
##
## Purpose of script:
#                   1. Retrieve in-field anomalies using only Vegetation indices derived from optical imagery.
#                   2. It stores the statistics per field into a hdf5 file. The stats include percecntage of aomalies,
#                   mean VI value per plot, mean VI values for the anomalous areas.

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

import rsgislib
import glob
import rsgislib.imagecalc
import rsgislib.imagecalc.calcindices
import rsgislib.imageutils
import rsgislib.rastergis
import os
from rsgislib import imageutils as imageutils
import datetime
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
from osgeo import gdal
from osgeo import gdalnumeric



###############################################################
####     Define common variables     #####

#Location of the root folder
m_root='/data'

## set working directory when using Docker
datatype = rsgislib.TYPE_32FLOAT
gdalformat = 'KEA'
Sensor = 'sen2'
prefix='SEN2'

StackVIFolder=m_root+'/StacksVIs'
AnomaliesPath=m_root+'/Anomalies/'

## List of all the VIs to be analised
VIList=['NDVI8A','NDVI8','EVI8A','EVI8', 'GCI8A', 'GCI8','GNDVI8A', 'GNDVI8','RECI8A_RE_B5', 'RECI8A_RE_B6',  'RECI8A_RE_B7', 'RECI8_RE_B5', 'RECI8_RE_B6',  'RECI8_RE_B7', \
        'RENDVI8A_RE_B5', 'RENDVI8A_RE_B6','RENDVI8A_RE_B7', 'RENDVI8_RE_B5', 'RENDVI8_RE_B6','RENDVI8_RE_B7', 'SAVI8A','SAVI8',\
        'NDWI8_SWIR1', 'NDWI8_SWIR2', 'NDWI8A_SWIR1', 'NDWI8A_SWIR2']




entrada_plot_histograms=input(" Do you want to plot the histograms Y/N: ")
print("You said ", entrada_plot_histograms)


def histAnalysis():
    print("\n HISTOGRAM ANALYSIS started at", datetime.datetime.now()," \n")
    import FunctionsHistoAnalysis
    driver = gdal.GetDriverByName(gdalformat)
    

    #Raster with rasterised crop plots
    clumps_image=m_root+'/Extents/Plots_Escobal_20191023_S2.kea'

    for VI in VIList:

        print("Running histogram analysis over ", VI, "images")

        #Output Anomalies VI_Image (For current Vegetation Index)
        Anom_class_image=m_root+'/Anomalies/'+prefix+'_'+VI+'stack_Anom.kea'

        #Stack VI image
        VI_Image =m_root+'/StacksVIs/'+prefix+'_'+VI+'_stack.kea'
        bandnames=imageutils.getBandNames(VI_Image)


        #HDF5 file that will store plots statistics
        histometrics=m_root+"/Anomalies/"+prefix+"_"+"HistoMetrics"+VI+".hdf5"

        #####################################
        # ££££ Run the function ££££ #        
        FunctionsHistoAnalysis.multiplehistothresholds (VI_Image, clumps_image, Anom_class_image, driver, \
                                                        histometrics, entrada_plot_histograms, m_root)
        print("Finished histogram analysis over ", VI, "images")
    print("\n ENDED HISTOGRAM ANALYSIS AT ", datetime.datetime.now())

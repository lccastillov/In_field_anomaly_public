########
### THis script gets thresholdsto distinguish between anomalous and  non-anomalous values  ####
####                            using histogram analysis                                   ####
###                It finally stores the statistics into a hdf5 file
#### In also turns into 99 the value of those plots that are partially or completely cover by clouds


from osgeo import gdal
from osgeo import gdalnumeric
import numpy as np
from rsgislib import imageutils
import scipy
from scipy import stats
from enum import Enum
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import rsgislib
import os


# Function to establish thresholds based on the histograms

#The outputs (thresholds) depend on the bin size

def PlotPlixels(SegID, ClumpsArray, ImageArray):

    class RSGISRATThresMeasure(Enum):
        kurtosis = 1
        skewness = 2
        combined = 3
        auto = 4
    
    # Get the raster coordinates of the pixels corresponding to the segment ID:
    y_idx, x_idx = np.where(ClumpsArray== SegID)

    # Get the individual pixel values:
    Pixel_Values = []
    for i in range(len(y_idx)):
        Row, Column = y_idx[i], x_idx[i]
        Pixel_Values.append(ImageArray[Row][Column])
        del Row, Column
    del y_idx, x_idx

    # Stores VI pixel values in an numpy array
    Pixel_Values = np.array(Pixel_Values)
    
    #print("Pixel_Values \n", Pixel_Values)
    
    #Retrieves the maximum, minimum, and mean values of the array of pixels whitin the current segment
    
    minPixelVal=(np.min(Pixel_Values))
    #minpixelval is a value that is returned as a list (the minimum value found in the array.

    sizeArrayPixels=Pixel_Values.size
    print("sizeArrayPixels", sizeArrayPixels, "type", type(sizeArrayPixels))
    return [minPixelVal, sizeArrayPixels]
    

def HistoThresholds(SegID, ClumpsArray, ImageArray, m_root, band_name, entrada_plot_histograms):
    
    class RSGISRATThresMeasure(Enum):
        kurtosis = 1
        skewness = 2
        combined = 3
        auto = 4
    
    # Get the raster coordinates of the pixels corresponding to the segment ID:
    y_idx, x_idx = np.where(ClumpsArray== SegID)

    # Get the individual pixel values:
    Pixel_Values = []
    for i in range(len(y_idx)):
        Row, Column = y_idx[i], x_idx[i]
        Pixel_Values.append(ImageArray[Row][Column])
        del Row, Column
    del y_idx, x_idx

    # Stores VI pixel values in an numpy array
    Pixel_Values = np.array(Pixel_Values)
    
    #print("Pixel_Values \n", Pixel_Values)
    
    #Retrieves the maximum, minimum, and mean values of the array of pixels whitin the current segment
    
    maxPixelVal=(np.max(Pixel_Values))
    minPixelVal=(np.min(Pixel_Values))
    meanPixelVal=(np.mean(Pixel_Values))
    stdPixelVal=(np.std(Pixel_Values))
    numpixels=np.size(Pixel_Values)
    

    #Retrieves the first and third quartiles of the array of pixels
    lq = np.percentile(Pixel_Values, 25)
    uq = np.percentile(Pixel_Values, 75)
    
    # Get shape of the Pixel_Values array
    n = Pixel_Values.shape[0]

    iqr = uq - lq
     
    #Optimal bin size as proposed by Freedman-Diaconis rule
    binSize = 2 * iqr * n**(-1/3)
    
    if binSize==0:
        
        binSize=1
        numBins=10

    #print("Bin Size = ", binSize)
    
    Pixel_Values = Pixel_Values[np.isfinite(Pixel_Values)]
    
    #Estimates the number of bins of the histogram (Total)
    numBins =  int((np.max(Pixel_Values) - np.min(Pixel_Values))/binSize)+2


    #print("The whole histogram has ",numBins, " bins")
    
    #Creates the histogram and the bin edges arrays
    hist, bin_edges = np.histogram(Pixel_Values, bins=numBins)
    
    #Calculates kurtosis and skewness for all the pixels within the plot
    PixelValsKurt=scipy.stats.kurtosis(Pixel_Values)
    
    #print("Pixel_Values kurtosis", PixelValsKurt, "\n")
    
    PixelValsSkew=scipy.stats.skew(Pixel_Values)
    #print("Pixel_Values skewness", PixelValsSkew, "\n")
    
    lqNumBins = int((lq - bin_edges[0])/binSize)+1
    uqNumBins = int((bin_edges[-1]-uq)/binSize)+1

    #Create numpy arrays to store metrics of partial histograms

    kurtosisVals = np.zeros((lqNumBins,uqNumBins), dtype=np.float)
    skewnessVals = np.zeros((lqNumBins,uqNumBins), dtype=np.float)
    meanVals= np.zeros((lqNumBins,uqNumBins), dtype=np.float)
    
    #Arrays to store the indices of each combination of bins
    lowBins = np.zeros((lqNumBins,uqNumBins), dtype=np.int)
    upBins = np.zeros((lqNumBins,uqNumBins), dtype=np.int)

    for lowBin in range(lqNumBins):
        for upBin in range(uqNumBins):
            #print("Bin [" + str(lowBin) + ", " + str(numBins-upBin) + "]")
            
            #Creates temporal histogram. A dynamic one that starts in the current low bien and ends in the upbin 
            histTmp = hist[lowBin:(numBins-upBin)]
            
            #print("Current temporary histogram \n",histTmp )
                        
            #Stores the current low bin and upbin in an array
            lowBins[lowBin,upBin] = lowBin
            upBins[lowBin,upBin] = numBins-upBin
            
            #Choose the thresholds to cut off value from the original array of pixels
            minLimit=bin_edges[lowBin+1]
            maxLimit=bin_edges[(numBins-upBin)]
            
            #Create an array with the pixels that are between the current thresholds
            subPixelsArr=np.extract((Pixel_Values>=minLimit)  & (Pixel_Values<=maxLimit),  Pixel_Values)
            
            
            #Stores the kurtosis and skewness of all the potential histograms in arrays
            #By default, Kurtosis is calculated using Fisher's definition (normal ==> 0.0). If False, Pearson’s definition is used (normal ==> 3.0)

            
            kurtosisVals[lowBin,upBin] = scipy.stats.kurtosis(subPixelsArr)
            skewnessVals[lowBin,upBin] = scipy.stats.skew(subPixelsArr)
            
            
            """
            This is wrong, because it estimates the kurtosis of the histogram, and not the kurtosis of the values within the thresholds
            kurtosisVals[lowBin,upBin] = scipy.stats.kurtosis(histTmp)
            skewnessVals[lowBin,upBin] = scipy.stats.skew(histTmp)
            """
                                 
            #print("kurtosisVals")
            #print(kurtosisVals[lowBin,upBin])
            
    ##### Estimates absolute values of Kurtosis and Skweness
    kurtosisValsAbs = np.absolute(kurtosisVals)
    skewnessValsAbs = np.absolute(skewnessVals)
    #print("Kurtosis Range: [" + str(np.min(kurtosisValsAbs)) + ", " + str(np.max(kurtosisValsAbs)) + "]") 
    #print("Skewness Range: [" + str(np.min(skewnessValsAbs)) + ", " + str(np.max(skewnessValsAbs)) + "]") 
    kurtosisValsNorm = (kurtosisValsAbs-np.min(kurtosisValsAbs)) / (np.max(kurtosisValsAbs)-np.min(kurtosisValsAbs))
    skewnessValsNorm = (skewnessValsAbs-np.min(skewnessValsAbs)) / (np.max(skewnessValsAbs)-np.min(skewnessValsAbs))

    #print("Kurtosis Norm Range: [" + str(np.min(kurtosisValsNorm)) + ", " + str(np.max(kurtosisValsNorm)) + "]") 
    #print("Skewness Norm Range: [" + str(np.min(skewnessValsNorm)) + ", " + str(np.max(skewnessValsNorm)) + "]") 

    #Combined is the sum of absolute values of skewness and kurtosis
    combined = kurtosisValsNorm + skewnessValsNorm
    #combined = kurtosisValsAbs + skewnessValsAbs
    #print(combined)

    #Retrieves the Indices of the minimum elements of the arrays that contain the min kurtosis, skewness and combined.
    #These indices (tuple) correspond to combination of upper and lower bins


    minKurt = np.unravel_index(np.argmin(kurtosisValsAbs, axis=None), kurtosisValsAbs.shape)
    minSkew = np.unravel_index(np.argmin(skewnessValsAbs, axis=None), skewnessValsAbs.shape)
    minComb = np.unravel_index(np.argmin(combined, axis=None), combined.shape)

    #print("Kurtosis bin indexes: ", minKurt)
    #print("Skewness bin indexes: ", minSkew)
    #print("Combined bin indexes: ", minComb)


    #minKurt[0]=row index ; minKurt[1]=column index
    #minkurt [0] muestra el indice de la bin en la que se deberia cortar el histograma para producir la minima kurtosis

    lowBinKurt = minKurt[0]
    #I look in the bin edges the low bin of the lower kurtosis and add half. 
    
    lowerThresKurt = bin_edges[lowBinKurt] + (binSize/2)
    #print("lowerThresKurt",lowerThresKurt)
    #upBinKurt=number of bins above the value. that has the lowest kurtosis ( in the right side)
    #upBinKurt [0] muestra la cantidad de bins que se deben quitar a la derecha para producir la minima kurtosis

    upBinKurt = numBins-minKurt[1]
    upperThresKurt = bin_edges[upBinKurt] + (binSize/2)
    print("No Change Data Range (Kurtosis): [" + str(lowerThresKurt) + "," + str(upperThresKurt) + "]")
    #print("upperThresKurt",upperThresKurt)
    
    
    lowBinSkew = minSkew[0]
    lowerThresSkew = bin_edges[lowBinSkew] + (binSize/2)
    upBinSkew = numBins-minSkew[1]
    upperThresSkew = bin_edges[upBinSkew] + (binSize/2)
    #print("No Change Data Range (Skewness): [" + str(lowerThresSkew) + "," + str(upperThresSkew) + "]")

    lowBinComb = minComb[0]
    lowerThresComb = bin_edges[lowBinComb] + (binSize/2)
    upBinComb = numBins-minComb[1]
    upperThresComb = bin_edges[upBinComb] + (binSize/2)
    print("No Change Data Range (Combined): [" + str(lowerThresComb) + "," + str(upperThresComb) + "]")


    lowerThres = 0.0
    upperThres = 0.0

    thresMeasure=RSGISRATThresMeasure.auto
    if thresMeasure == RSGISRATThresMeasure.kurtosis:
        lowerThres = lowerThresKurt
        upperThres = upperThresKurt
    elif thresMeasure == RSGISRATThresMeasure.skewness:
        lowerThres = lowerThresSkew
        upperThres = upperThresSkew
    elif thresMeasure == RSGISRATThresMeasure.combined:
        lowerThres = lowerThresComb
        upperThres = upperThresComb
    elif thresMeasure == RSGISRATThresMeasure.auto:
        """Quiere decir que si en alguno de los extremos, la diferencia los thresholds obtenidos from
        skewness y kurtosis es muy grande (>interquartile range -iqr) se da prioridad a los 
        thresholds asociados a la skewness y son los que se definen como thresholds
        Estos thresholds representan dónde se deberia cortar el histograma original para minimizar ya sea la kurtosis o la skewness
        Hay uno superior y uno inferior, ya que pueden haber varias combinaciones en ambas colas
        En caso que esta diferencia en los extremos no sea tan grande, entonces se utilizan los valores combinados
        Es decir those combinations that minimize both kurtosis and skewness
        """
        if (abs(lowerThresKurt-lowerThresSkew) > (uq-lq)) or (abs(upperThresKurt-upperThresSkew) > (uq-lq)):
            lowerThres = lowerThresSkew
            upperThres = upperThresSkew
        else:
            lowerThres = lowerThresComb
            upperThres = upperThresComb
        print("No Change Data Range (auto): [" + str(lowerThres) + "," + str(upperThres) + "]")
    else:
        raise Exception("Don't understand metric for threshold provided must be of type ThresMeasure")

    OptimalArray=np.extract((Pixel_Values>=lowerThres)  & (Pixel_Values<=upperThres),  Pixel_Values)
    print("There are", OptimalArray.size, " pixels that are considered not anomalous")

    #Get the number o pixels below the lower threshold
    LowAnomalies=Pixel_Values < lowerThres
    numLowAnomalies=len(Pixel_Values[LowAnomalies])
    print("There are", numLowAnomalies, " anomalous pixels below the lower threshold ")
    
    #Get the number o pixels above the upper threshold
    UpperAnomalies=Pixel_Values > upperThres
    numUpperAnomalies=len(Pixel_Values[UpperAnomalies])
    print("There are", numUpperAnomalies, " anomalous pixels above the upper threshold ")
    percLowAnomalies=numLowAnomalies/Pixel_Values.size
    print("percLowAnomalies",percLowAnomalies)
    percUpperAnomalies=numUpperAnomalies/Pixel_Values.size

    # In case the user choses to plot the histograms
    if entrada_plot_histograms == "Y":
        # print("You chose to print the histograms")
        plotHistoStatistics(lowerThres, upperThres, lowerThresKurt, upperThresKurt, lowerThresSkew, upperThresSkew,
                            lowerThresComb, upperThresComb, bin_edges, hist, binSize, m_root, band_name, SegID)

    #Find th statisticas for the pixels that generate the most normal distribution
    maxOptHist=(np.max(OptimalArray))
    minOptHist=(np.min(OptimalArray))
    meanOptHist=(np.mean(OptimalArray))
    stdOptHist=(np.std(OptimalArray))
    modeOptHist=stats.mode(OptimalArray)
    kurtosisOptHist=scipy.stats.kurtosis(OptimalArray)
    skewnessOptHist=scipy.stats.skew(OptimalArray)
    
    
    maxPixelVal=(np.max(Pixel_Values))
    minPixelVal=(np.min(Pixel_Values))
    meanPixelVal=(np.mean(Pixel_Values))
    stdPixelVal=(np.std(Pixel_Values)) 
    
    
    return [lowerThres, upperThres,lowerThresKurt , upperThresKurt, lowerThresSkew , upperThresSkew, lowerThresComb, upperThresComb, bin_edges, hist, binSize, maxPixelVal, minPixelVal, meanPixelVal, stdPixelVal, \
          maxOptHist, minOptHist, meanOptHist,  stdOptHist, modeOptHist, kurtosisOptHist,  skewnessOptHist, PixelValsKurt, PixelValsSkew, numpixels, numLowAnomalies ,numUpperAnomalies,percLowAnomalies ,percUpperAnomalies]

def plotHistoStatistics(lowerThres, upperThres,lowerThresKurt , \
                        upperThresKurt, lowerThresSkew , upperThresSkew, \
                        lowerThresComb, upperThresComb, bin_edges, hist, binSize, m_root, band_name, clump ):

    exportplot = m_root + '/Ibague/Histograms/PS_B_' + str(band_name) + 'Plot_' + str(clump) + '.jpg'
    print("Se va a imprimir el histogramas en ", exportplot)

    center = (bin_edges[:-1] + bin_edges[1:]) / 2
    plt.bar(center, hist, align='center', width=binSize)
    showAllThreshPlot=True
    if showAllThreshPlot:
        plt.vlines(lowerThresKurt, 0, np.max(hist), color='y', linewidth=1, label='Kurtosis Lower')
        plt.vlines(upperThresKurt, 0, np.max(hist), color='y', linewidth=1, label='Kurtosis Upper')
        plt.vlines(lowerThresSkew, 0, np.max(hist), color='r', linewidth=1, label='Skewness Lower')
        plt.vlines(upperThresSkew, 0, np.max(hist), color='r', linewidth=1, label='Skewness Upper')
        plt.vlines(lowerThresComb, 0, np.max(hist), color='g', linewidth=1, label='Combined Lower')
        plt.vlines(upperThresComb, 0, np.max(hist), color='g', linewidth=1, label='Combined Upper')
    else:
        plt.vlines(lowerThres, 0, np.max(hist), color='r', linewidth=1, label='Lower Threshold')
        plt.vlines(upperThres, 0, np.max(hist), color='r', linewidth=1, label='Upper Threshold')
    plt.grid(True)
    plt.legend(loc=0)
    plt.savefig(exportplot)
    plt.close()


def multiplehistothresholds (VI_Image, clumps_image, Anom_class_image, driver, histometrics, entrada_plot_histograms, m_root):


    #Bring dates
    bandnames=imageutils.getBandNames(VI_Image)
    print("Band names", bandnames)
    
    
    #First I open the image in gdal 
    RasterImage = gdal.Open(VI_Image)
    
    #I create a copy of the image in which I will store the classified image
    PixelClumps=driver.CreateCopy(Anom_class_image,RasterImage , strict=0)
    
    
    #I open the clumps image in gdal 
    RasterClumps = gdal.Open(clumps_image)
    
    # Read the clumps image as an array
    ClumpsArray=RasterClumps.GetRasterBand(1).ReadAsArray()
    
    
    #Identify the total number of plots (clumps) available
    N_Segs = np.max(ClumpsArray)
    
    
    #Creating panda dataframes for lower and upper thresholds
    LowerArray = pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    UpperArray = pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    
    lowerThresKurtArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    upperThresKurtArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    
    lowerThresSkewArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    upperThresSkewArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    
    lowerThresCombArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    upperThresCombArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    
    
    #Creating panda dataframes for mean values of all the pixels within each crop plot
    meanValsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    maxValsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    minValsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    stdValsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    skewValsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    kurtValsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    NumPixelsArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    numLowAnomaliesArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    numUpperAnomaliesArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    percLowAnomaliesArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    percUpperAnomaliesArray=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    
    #Creating panda dataframes for statistics of values within the normal distribution pf the plot
    
    meanArray_ND=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    maxArray_ND=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    minArray_ND=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    stdArray_ND=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    skewArray_ND=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    kurtArray_ND=pd.DataFrame(0, index=np.arange(1,N_Segs+1,1), columns=bandnames, dtype=np.float64)
    
      
    print("datatypes",LowerArray.dtypes)
    for band in range(RasterImage.RasterCount) :
    #for band in range(3,4) :
        arrayposit=band
        band += 1
        ImageArray=RasterImage.GetRasterBand(band).ReadAsArray()
        
        band_name=bandnames[arrayposit]
            
        #Create a zero array that will overwrite PixelClumps-- IT might be moved outside the for when working with more plots
            
        OutputArray = np.zeros_like(RasterClumps, dtype='uint8')
        np.max(OutputArray)
    
           
        ImageArray=RasterImage.GetRasterBand(band).ReadAsArray()
        
        #for clump in range(200, 220):
        for clump in range(1, N_Segs+1):
    
            Lowerlist = []
            Upperlist = []
    
            print("GETTING BAND: ",band_name, "band_number",band, "polygon",clump)
    
    
            #Get raster band to be written
            PixelClumpsBand=PixelClumps.GetRasterBand(band)
    
            minPixelVal, sizeArrayPixels=PlotPlixels(SegID=clump, ClumpsArray=ClumpsArray, ImageArray=ImageArray)
            #print("minoixelval",minPixelVal[0],  type(minPixelVal))
            print("sizeArrayPixels", sizeArrayPixels, "type", type(sizeArrayPixels))

            
            if (minPixelVal == -999 or sizeArrayPixels<=5):


                print("the clump", clump, "has a cloud")
                # Apply classification thresholds to the output image:
                OutputArray = np.where((OutputArray == 0) & (ClumpsArray == clump), -999, OutputArray)

                #print("Hey guera do not forget to remove the #####")
                PixelClumpsBand.WriteArray(OutputArray)
        
                #print("final output array /n", OutputArray)
                
                
                
                #Store the thresholds for each band and each date in a panda dataframe
        
                LowerArray.loc[clump][str(band_name)] = -999
                UpperArray.loc[clump][str(band_name)] = -999
        
                lowerThresKurtArray.loc[clump][str(band_name)] = -999
                upperThresKurtArray.loc[clump][str(band_name)] = -999
        
                lowerThresSkewArray.loc[clump][str(band_name)] = -999
                upperThresSkewArray.loc[clump][str(band_name)] =-999
        
                lowerThresCombArray.loc[clump][str(band_name)] = -999
                upperThresCombArray.loc[clump][str(band_name)] = -999
                meanValsArray.loc[clump][str(band_name)] = -999
                minValsArray.loc[clump][str(band_name)] =-999
                maxValsArray.loc[clump][str(band_name)] = -999
                #modeValsArray.loc[clump][str(band_name)] = modePixelVal
                stdValsArray.loc[clump][str(band_name)] = -999
                kurtValsArray.loc[clump][str(band_name)] = -999
                skewValsArray.loc[clump][str(band_name)] = -999
                NumPixelsArray.loc[clump][str(band_name)]=-999
                numLowAnomaliesArray.loc[clump][str(band_name)]=-999
                numUpperAnomaliesArray.loc[clump][str(band_name)]=-999
                percLowAnomaliesArray.loc[clump][str(band_name)]=-999
                percUpperAnomaliesArray.loc[clump][str(band_name)]=-999

      
                maxArray_ND.loc[clump][str(band_name)]= -999
                meanArray_ND.loc[clump][str(band_name)]= -999
                minArray_ND.loc[clump][str(band_name)]= -999
                stdArray_ND.loc[clump][str(band_name)]= -999
                skewArray_ND.loc[clump][str(band_name)]= -999
                kurtArray_ND.loc[clump][str(band_name)]= -999
            else:
                
            

                lowerThres, upperThres, lowerThresKurt, upperThresKurt, lowerThresSkew, upperThresSkew, lowerThresComb, upperThresComb, \
                bin_edges, hist, binSize, maxPixelVal, minPixelVal, meanPixelVal, stdPixelVal, maxOptHist, minOptHist, meanOptHist,  \
                stdOptHist, modeOptHist, kurtosisOptHist,  skewnessOptHist,PixelValsKurt, PixelValsSkew, numpixels, numLowAnomalies ,\
                numUpperAnomalies,percLowAnomalies ,percUpperAnomalies=\
                    HistoThresholds(SegID=clump, ClumpsArray=ClumpsArray, ImageArray=ImageArray, m_root=m_root, band_name=band_name, entrada_plot_histograms=entrada_plot_histograms)
                
                
        
                # Apply classification thresholds to the output image:
                OutputArray = np.where((OutputArray == 0) & (ClumpsArray == clump) & (ImageArray < lowerThres), 1, OutputArray) # Classify pixels below lower threshold.
                OutputArray = np.where((OutputArray == 0) & (ClumpsArray == clump) & (ImageArray > upperThres), 3, OutputArray) # Classify pixels above upper threshold.
                OutputArray = np.where((OutputArray == 0) & (ClumpsArray == clump) & (ImageArray >= lowerThres) & (ImageArray <= upperThres), 2, OutputArray) # Classify pixels near the centre of the distribution.
        
        
                # write the output array:
                
                #print("Hey guera do not forget to remove the #####")
                PixelClumpsBand.WriteArray(OutputArray)
        
                print("final output array /n", OutputArray)
                             
                
                #Store the thresholds for each band and each date in a panda dataframe
        
                LowerArray.loc[clump][str(band_name)] = lowerThres
                UpperArray.loc[clump][str(band_name)] = upperThres
        
                lowerThresKurtArray.loc[clump][str(band_name)] = lowerThresKurt
                upperThresKurtArray.loc[clump][str(band_name)] = upperThresKurt
        
                lowerThresSkewArray.loc[clump][str(band_name)] = lowerThresSkew
                upperThresSkewArray.loc[clump][str(band_name)] = upperThresSkew
        
                lowerThresCombArray.loc[clump][str(band_name)] = lowerThresComb
                upperThresCombArray.loc[clump][str(band_name)] = upperThresComb
                meanValsArray.loc[clump][str(band_name)] = meanPixelVal
                minValsArray.loc[clump][str(band_name)] = minPixelVal
                maxValsArray.loc[clump][str(band_name)] = maxPixelVal
                #modeValsArray.loc[clump][str(band_name)] = modePixelVal
                stdValsArray.loc[clump][str(band_name)] = stdPixelVal
                kurtValsArray.loc[clump][str(band_name)] = PixelValsKurt
                skewValsArray.loc[clump][str(band_name)] = PixelValsSkew
                NumPixelsArray.loc[clump][str(band_name)]=numpixels
                numLowAnomaliesArray.loc[clump][str(band_name)]=numLowAnomalies
                numUpperAnomaliesArray.loc[clump][str(band_name)]=numUpperAnomalies
                percLowAnomaliesArray.loc[clump][str(band_name)]=percLowAnomalies
                percUpperAnomaliesArray.loc[clump][str(band_name)]=percUpperAnomalies
      
                maxArray_ND.loc[clump][str(band_name)]= maxOptHist
                meanArray_ND.loc[clump][str(band_name)]= meanOptHist
                minArray_ND.loc[clump][str(band_name)]= minOptHist
                stdArray_ND.loc[clump][str(band_name)]= stdOptHist
                skewArray_ND.loc[clump][str(band_name)]= skewnessOptHist
                kurtArray_ND.loc[clump][str(band_name)]= kurtosisOptHist
        

            
        del lowerThres, upperThres, PixelClumpsBand
        print("Poniendo los nombres de las bandas")
    
        imageutils.setBandNames(Anom_class_image, bandnames)
        rsgislib.imageutils.popImageStats(Anom_class_image, True, 0., True)
       
    
        #minimosarray=minimosarray.append(pd.Series(Lowerlist, index=bandnames),ignore_index=True)
    #print("LowerArray",LowerArray)
    #print("UpperArray",UpperArray)
    
    del  RasterClumps, RasterImage, PixelClumps, VI_Image

    #############################################################################################
    ####################     Store in hdf5 file #################################################
    
    
    """
    mode : {‘a’, ‘w’, ‘r+’}, default ‘a’
    Mode to open file:
    ‘w’: write, a new file is created (an existing file with the same name would be deleted).
    ‘a’: append, an existing file is opened for reading and writing, and if the file does not exist it is created.
    ‘r+’: similar to ‘a’, but the file must already exist
    """
    
    #Storing the data in hdf5 files. In HDFView chesk the axis 0 and axis1 folders. 
    #it seems that the index starts in 0, but it starts in 1
    #To store the lower and upper thresholds finally chosen
    LowerArray.to_hdf(histometrics,key='Thresholds/LowerThreshold',mode='a',table=True, data_columns=True)
    UpperArray.to_hdf(histometrics,key='Thresholds/UpperThreshold',mode='a',table=True, data_columns=True)
    
    #To store the lower and upper thresholds only based on the minimum kurtosis
    lowerThresKurtArray.to_hdf(histometrics,key='Thresholds/LowerThresKurtosis',mode='a',table=True, data_columns=True)
    upperThresKurtArray.to_hdf(histometrics,key='Thresholds/UpperThresKurtosis',mode='a',table=True, data_columns=True)
    
    
    #To store the lower and upper thresholds only based on the minimum skewness
    lowerThresSkewArray.to_hdf(histometrics,key='Thresholds/LowerThresSkewness',mode='a',table=True, data_columns=True)
    upperThresSkewArray.to_hdf(histometrics,key='Thresholds/UpperThresSkewness',mode='a',table=True, data_columns=True)
    
    lowerThresCombArray.to_hdf(histometrics,key='Thresholds/LowerThresCombined',mode='a',table=True, data_columns=True)
    upperThresCombArray.to_hdf(histometrics,key='Thresholds/UpperThresCombined',mode='a',table=True, data_columns=True)
    
    
    #To store the lower and upper thresholds only based on the mean value per plot per date
    meanValsArray.to_hdf(histometrics,key='PlotStats/PlotMeanVal',mode='a',table=True, data_columns=True)
    minValsArray.to_hdf(histometrics,key='PlotStats/PlotMinVal',mode='a',table=True, data_columns=True)
    maxValsArray.to_hdf(histometrics,key='PlotStats/PlotMaxVal',mode='a',table=True, data_columns=True)
    stdValsArray.to_hdf(histometrics,key='PlotStats/PlotStdVal',mode='a',table=True, data_columns=True)
    skewValsArray.to_hdf(histometrics,key='PlotStats/PlotSkewVal',mode='a',table=True, data_columns=True)
    kurtValsArray.to_hdf(histometrics,key='PlotStats/PlotKurtVal',mode='a',table=True, data_columns=True)
    NumPixelsArray.to_hdf(histometrics,key='PlotStats/NumPixels',mode='a',table=True, data_columns=True)
    numLowAnomaliesArray.to_hdf(histometrics,key='PlotStats/numLowAnomalies',mode='a',table=True, data_columns=True)
    numUpperAnomaliesArray.to_hdf(histometrics,key='PlotStats/numUpperAnomalies',mode='a',table=True, data_columns=True)
    percLowAnomaliesArray.to_hdf(histometrics,key='PlotStats/percLowAnomalies',mode='a',table=True, data_columns=True)
    percUpperAnomaliesArray.to_hdf(histometrics,key='PlotStats/percUpperAnomalies',mode='a',table=True, data_columns=True)  

    
    meanArray_ND.to_hdf(histometrics,key='NormalPlotStats/NormalMean',mode='a',table=True, data_columns=True)
    maxArray_ND.to_hdf(histometrics,key='NormalPlotStats/NormalMax',mode='a',table=True, data_columns=True)
    minArray_ND.to_hdf(histometrics,key='NormalPlotStats/NormalMin',mode='a',table=True, data_columns=True)
    stdArray_ND.to_hdf(histometrics,key='NormalPlotStats/NormalStd',mode='a',table=True, data_columns=True)
    skewArray_ND.to_hdf(histometrics,key='NormalPlotStats/NormalSkew',mode='a',table=True, data_columns=True)
    kurtArray_ND.to_hdf(histometrics,key='NormalPlotStats/NormalKurt',mode='a',table=True, data_columns=True)
    
    
    #Store statistics for new histograms
    
    h = h5py.File(histometrics)
    
    h.close()




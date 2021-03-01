def import_hdf5_anomalies(m_root, VI, prefix, df_age, hdf5_field_anomalies, anomalies_folder):
    #imports hdf5 file with info onf in-field anomalies per plot
    # and matches them with the following field data:
    # Emergence and harvest dates, variety, plant density.
    import pandas as pd
    import numpy as np
    hdf5_file = anomalies_folder + '/'+ prefix + '_HistoMetrics' + VI + '.hdf5'
    csv_image_anomalies_vs_plots = m_root + '/Ibague/Anomalies_interpretation_grid/' + prefix +'_'+VI+ '_csv_image_anomalies_vs_plots.csv'
    csv_image_anomalies_vs_plots_v2=m_root + '/Ibague/Anomalies_interpretation_grid/' + prefix +'_'+VI+ '_csv_image_anomalies_vs_plots_v2.csv'
    df = df_age.copy()
    # Retrieves the number of anomalies per plot age

    # Read hdf file into pandas dataframe
    hdf = pd.HDFStore(hdf5_file, mode='r')

    # Shows all the keys within the hdf5 file
    # print("hdfkeys",hdf.keys())

    # Retrieve percentage of Lower and upper anomalies per plot per date (from hdf5 file)
    df_mean_vi = hdf.get('/PlotStats/PlotMeanVal')
    df_norm_mean_vi = hdf.get('/NormalPlotStats/NormalMean')
    df_low_mean_vi = hdf.get('/LowPlotStats/low_Mean')
    df_up_mean_vi = hdf.get('/UpPlotStats/up_Mean')
    df_std_vi = hdf.get('/PlotStats/PlotStdVal')

    DF_percLowerAnomalies = hdf.get('/PlotStats/percLowAnomalies')
    DF_numUpperAnomalies = hdf.get('/PlotStats/percUpperAnomalies')

    # Create empty lists to store temporary lower and upper anomalies
    mean_vi_Array = []
    df_std_vi_Array=[]
    mean_normal_vi_array = []
    mean_low_vi_array = []
    mean_up_vi_array = []
    perc_lowerAnomalies = []
    perc_upperAnomalies = []

    # Loop to go through all the indexes of the query and retrieve the respective percentage of anomalies
    for index, row in df.iterrows():
        plot = row['Id_plot']
        date = row['Image_date'].replace("-", "")

        meanNDVI = df_mean_vi.loc[int(float(plot))][str(date)]
        mean_normal_vi = df_norm_mean_vi.loc[int(float(plot))][str(date)]
        low_mean_VI = df_low_mean_vi.loc[int(float(plot))][str(date)]
        up_mean_VI = df_up_mean_vi.loc[int(float(plot))][str(date)]
        percLowerAnomalies = DF_percLowerAnomalies.loc[int(float(plot))][str(date)]
        numUpperAnomalies = DF_numUpperAnomalies.loc[int(float(plot))][str(date)]
        std_vi=df_std_vi.loc[int(float(plot))][str(date)]

        #print("mean",meanNDVI, "low",low_mean_VI, "up_mean_VI",up_mean_VI,'percLowerAnomalies',percLowerAnomalies,"numUpperAnomalies", numUpperAnomalies)

        # Show which of the plots are potentially clouded for specific dates
        # if meanNDVI == -999:
        # print("The plot", row['Codigo'], " is covered by clouds on ", row['Image_date'])

        # Append the current values to the lists of lower and upper anomalies respectively
        mean_vi_Array.append(meanNDVI)
        df_std_vi_Array.append(std_vi)
        perc_lowerAnomalies.append(percLowerAnomalies)
        perc_upperAnomalies.append(numUpperAnomalies)
        mean_normal_vi_array.append(mean_normal_vi)
        mean_low_vi_array.append(low_mean_VI)
        mean_up_vi_array.append(up_mean_VI)

    # Create  new columns in the dataframe that will store the percentage of lower and uper anomalies
    mean_vi_colum_title = 'mean_' + VI
    df['mean_vi'] = mean_vi_Array
    df['std_vi']=df_std_vi_Array
    df['perc_LowerAnomalies'] = perc_lowerAnomalies
    df['perc_UpperAnomalies'] = perc_upperAnomalies
    df['normalmean_vi'] = mean_normal_vi_array
    df['lowmean_vi'] = mean_low_vi_array
    df['upmean_vi'] = mean_up_vi_array
    #print("df_preliminar", df.head(10))

    # Remove the plots potentially clouded
    df = df[df.mean_vi != -999]
    #I commmented this to avoid problems in vis with band7
    #df = np.where((df.mean_vi > 0.000000001) & (df.mean_vi < -0.00001) )
    df = df[(df.mean_vi > 0.000000001) | (df.mean_vi < -0.00001)]
    #df[(df['Salary_in_1000'] >= 100) & (df['Age'] < 60) & df['FT_Team'].str.startswith('S')]

    # We create a new pandas dataframe that will onlye keep the records with more than 0 anomalies (we will call it Standardized)
    df_std = df[df.perc_LowerAnomalies > 0.000001]
    df_std = df[df.perc_UpperAnomalies > 0.000001]
    df.to_hdf(hdf5_field_anomalies, key='Anom_images/' + prefix + '/Plot_ages_vs_anom_image', mode='a', table=True,
              data_columns=True)
    df.to_csv(csv_image_anomalies_vs_plots_v2, index=False)

    return [df_std, df]
    hdf5_field_anomalies.close()


def plots_age(Database, imageList):
    import sqlite3
    import pandas as pd
    import numpy as np
    print("imageListc", imageList)
    # Query that presents info per plot and cycles dates, including plot age when the images were taken
    # AnomImage=Stack Image with anomalies
    # Database=DAtabase
    # phen_stage=Desired phenological stage

    conn = sqlite3.connect(Database)
    c = conn.cursor()
    conn.commit()

    # Create an array to store the selected records
    plot_info = np.empty([0, 10])
    names = ["Codigo", "Id_plot", "Image_date", "age", "Emergence date", "Emergence month", "Harvest_date", "Variety",
             "Seed_Density_kg_Ha", "Yield"]
    print("imageListb", imageList)
    for id in imageList:
        # Watch out the date format must be with the symbol "2012-03-14"

        id = id.date()

        c.execute("""
                             SELECT 
                              Ibague_Plots.Codigo,
                              Ibague_Plots.Id_plot,  
                              ? AS "Image_date", 
                              julianday(?) - julianday(Ibague_Cycles.Emergence_Date) AS "age",
                              Ibague_Cycles.Emergence_Date AS "Emergence date",
                              strftime('%m', Ibague_Cycles.Emergence_Date) AS "Emergence month",
                              Ibague_Cycles.Harvest_Date AS "Harvest_date",
                              Ibague_Cycles.Variety,
                               Ibague_Cycles.Seed_Density_kg_Ha, 
                               Ibague_Cycles.Yield


                          FROM Ibague_Plots
                          INNER JOIN Ibague_Cycles ON  Ibague_Cycles.Id_plot = Ibague_Plots.Id_plot


                          WHERE
                            julianday(?) >= julianday(Ibague_Cycles.Emergence_Date)
                            AND
                             (julianday(Ibague_Cycles.Harvest_Date)>=julianday(?)
                             OR 
                            julianday(Ibague_Cycles.Harvest_Date)=0)

                            """
                  , (id, id, id, id))
        alist = c.fetchall()
        len(alist)
        a = np.array(alist)

        if a.size != 0:
            plot_info = np.concatenate((plot_info, a), axis=0)
    df_age = pd.DataFrame(plot_info, columns=names)
    ## remove duplicates
    df_age = df_age.drop_duplicates(subset=['Id_plot', 'Emergence date', 'Image_date'])

    df_age=df_age[['Codigo','Id_plot','Image_date','Emergence date','Emergence month','Harvest_date','age',\
        'Variety','Seed_Density_kg_Ha','Yield']]

    return [df_age]

def conditions_stages (df):
    #It classifies the phenological stages based on the plot age
    #Check https://www.uaex.edu/publications/pdf/mp192/chapter-2.pdf
    def f(row):
        if row['age'] <= 12:
            val = 'Seedling'
        elif row['age'] <= 40 :
            val = 'Tillering'
        elif row['age'] <= 45:
            val = 'Vegetative lag'
        elif row['age'] <= 65:
            val = 'Panicle growth'
        elif row['age'] <= 95:
            val = 'Booting'
        else:
            val = 'Ripening'
        return val

    df['stage']=df.apply(f, axis=1)
    return[df]

def retrieve_pixelvalues_in_array(SegID, ClumpsArray, ImageArray):
    #It retrieves all the VI pixel values in an array to be further analise
    import numpy as np
    # Get the raster coordinates of the pixels corresponding to the segment ID:
    y_idx, x_idx = np.where(ClumpsArray == SegID)

    # Get the individual pixel values:
    Pixel_Values = []
    for i in range(len(y_idx)):
        Row, Column = y_idx[i], x_idx[i]
        Pixel_Values.append(ImageArray[Row][Column])
        del Row, Column
    del y_idx, x_idx

    # Stores VI pixel values in an numpy array
    Pixel_Values = np.array(Pixel_Values)
    sizeArrayPixels = Pixel_Values.size
    if sizeArrayPixels <= 5 :
        validplot = 'no'
    else:
        validplot = 'si'
    # print("Pixel_Values \n", Pixel_Values)


    return [validplot, Pixel_Values]

def retrieve_mean_yields_plot(yield_pr_images_list,plot_yields_folder,df_age,m_root, vi_temp, prefix,\
                                                                   hdf5_field_anomalies,
                                                                   anomalies_folder,ClumpsArray):
    #Retrieves the mean yield per plot using the data obtained from the sensour mounted on the rice harvester

    import pandas as pd
    import numpy as np
    import gdal
    import os
    #Funtion that creates a dataframe with the yield per plot
    #calculated as the mean off all yield pixels within the plot

    df_yield_mean = pd.DataFrame(columns=['Id_plot', 'Harvest_date', 'av_yield'])

    for yield_pr_image in yield_pr_images_list:

        name_file = str(yield_pr_image).replace(plot_yields_folder, "").replace("/", "").replace("_", "")

        # Identify current harvest date and plot name
        current_harvest_date = name_file[-14:-6]
        print("starting analysis current_harvest_date", current_harvest_date)
        current_plot_name = str(name_file).replace(current_harvest_date, "").replace("yield", "").replace("pr.tif", "")
        print("starting analysis plot name= ", current_plot_name)

        # first I choose the rows that have this id in the df_age
        id_lote = df_age[df_age.Codigo == current_plot_name].copy()
        # print("id_lote",id_lote)
        id_lote = int(id_lote['Id_plot'].iloc[0])
        folder_yield_plot = plot_yields_folder + '/plot_' + current_plot_name + '_' + current_harvest_date
        if not os.path.exists(folder_yield_plot):
            os.makedirs(folder_yield_plot)

        print("current_plot_name", current_plot_name, " current_harvest_date", current_harvest_date)

        # I open yield raster for the current plot
        yield_open_raster = gdal.Open(yield_pr_image)

        # Retrieve yield array for current image
        yield_array = yield_open_raster.GetRasterBand(1).ReadAsArray()
        validplot, yield_pixel_values = retrieve_pixelvalues_in_array(id_lote, ClumpsArray,
                                                                                    yield_array)
        mean_yield_plot = np.mean(yield_pixel_values)
        print("yield_pixel_values shape", yield_pixel_values.shape)

        ##Change date data to remove dash
        df_age['Image_date_no_dash'] = df_age['Image_date'].str.replace("-", "")
        df_age['Harvest_date_no_dash'] = df_age['Harvest_date'].str.replace("-", "")


        # Get pandas dataframe with info mean values of VIS and other stats for the age associated to each image
        df_std_temp, df_temp = import_hdf5_anomalies(m_root, vi_temp, prefix, df_age,
                                                                   hdf5_field_anomalies,
                                                                   anomalies_folder)
        df_current_plot = df_temp[
            (df_temp['Id_plot'] == id_lote) & (df_temp['Harvest_date_no_dash'] == current_harvest_date)]

        if df_current_plot.empty == False:
            ##Retreive emergence date for each plot for which there is available yield
            emergence_date = df_current_plot['Emergence date'].values[0]
            print("nombre_lote", current_plot_name, " id ", id_lote, " emergence_date ", emergence_date)
            new_row = {'Id_plot': id_lote, 'Harvest_date': current_harvest_date, 'av_yield': mean_yield_plot}
            df_yield_mean = df_yield_mean.append(new_row, ignore_index=True)

        print("df_yield_mean", df_yield_mean.head(2))



    return[df_yield_mean]

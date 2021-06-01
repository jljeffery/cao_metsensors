
# =========================================================================
# MET SENSORS PYTHON MODULE
# Tools to produces NetCDF data from a day of time series met data files
#
# Author:        Judith Jeffery, RAL 
# History:       Based on Python code written for Leosphere lidar by Chris Walden 
# Version:	 1.0
# Last modified: 09/08/17
# =========================================================================
module_version = 1.0

# -------------------------------------------------------------------------
# Import required tools
# -------------------------------------------------------------------------
import numpy as np
import os, re, sys, getopt, shutil, zipfile, string, pwd
import netCDF4 as nc4

import datetime
import time
import calendar
import scipy.signal
from pylab import *
import module_data_object_python3
import module_distrometer_format5
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
# ----------------------
# Define useful function
# ----------------------
def in_interval(seq,xmin,xmax):
    for i, x in enumerate(seq):
        if x>=xmin and x<xmax:
            yield i

# ------------------------------------------------------------------------
# Pluvio raingauge netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_pluvio(nday):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon

    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['pldcrg_ch','pldcrg_ch']	#Need a file for each dataset read, but can be same one. May want them separate for different corrections 
    nvar = np.size(chids)
    n = 0 	#Number of data points

    # -------------------------
    # Define various parameters
    # -------------------------

    missing_value = -1.0E+20


    # ----------------------
    # Initialize data arrays
    # ----------------------
    vals        = missing_value * np.ones((10000,nvar))
    timesecs    = np.zeros((10000))
    input_missing =  np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_moe_home2/range/pluvio/"
    #path_out = "/data/netCDF/files/amof-pluvio_chilbolton/"
    #graph_path_out = "/data/netCDF/graphs/amof-pluvio_chilbolton/" + str(year) + '/' + strmon + '/'        #Still need to add directory for name of product
    #path_out = "/data/netCDF/files/ncas-rain-gauge-4/ncas-rain-gauge-4_cao_"
    path_out = "/data/amof-netCDF/ncas-rain-gauge-4/ncas-rain-gauge-4_cao_"
    data_version = "_v1.0"	#Changed temporarily so as to show this was generated with v3 python
    graph_path_out = "/data/amof-netCDF/graphs/ncas-rain-gauge-4" + data_version + "/" + str(year) + '/' + strmon + '/'
    data_product = "_precipitation"
    #path_out = "/home/jla/netCDF/files/amof-pluvio_chilbolton/"
    #text_path_out = "/mnt/wilma_home/jla/halo/text_profile/"
    print(path_in, graph_path_out)

    template_file_path = "/home/jla/python/metsensors/ncas-rain-gauge-4_metadata-template.yaml"      #YAML template

    #------------Setting up YAML data object----------------------
    #Need export PYTHONPATH=/home/jla/python/global_modules
    handler = module_data_object_python3.Handler()
    data_object_id = handler.load_template(template_file_path)
    print(data_object_id)
    #handler.show_requirements_for_template(data_object_id)

    # ---------------------------------------------------------------
    # Loop to process raw data from current and next day's data files
    # Since first data point in file is at 00:00:00, first data point
    # contains data entirely from the previous day.
    # Hence, need to open file from next day too.
    # ---------------------------------------------------------------

    for day_incr in range(2):

        nday_file = nday + day_incr


        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        print(nday_file, datevals)


        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "pldc_" + datestring_now + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print(infiles)
            nfiles = len(infiles)
            if nfiles == 0:
                print('No files on this day')

        else:
            nfiles = 0
            print('No directory found')

        # ----------------------
        # Process data files
        # ----------------------
        print('Starting to process files')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            for z in range(9):      #9 header lines, may eventually want some of them
                line = f.readline()

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #print fdata
                fdatatmp = line.split()
                #print len(line), line
                startnum = date2num(datetime.datetime(int(year_now),int(month_now),int(day_now),int(line[6:8]),int(line[9:11]),int(line[12:14])))
                if startnum > nday and startnum <= (nday + 1):
                    timesecs[n] = 3600.0*float(line[6:8]) + 60.0*float(line[9:11]) + float(line[12:14])
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0
                    if len(line) > 45:	#There's data in the line beyond the quantities we're writing to a file
                        fdata=fdatatmp[1].split(';')
                        tempstr = fdata[0]
                        char_flag = 0
                        for m in range(len(tempstr)):
                            char_test = tempstr[m].isalpha()
                            if char_test:	#Alphabet character found
                                char_flag = char_flag + 1
                        if tempstr[0] == '+' and char_flag == 0:	#Lines with errors don't tend to start records with + 
                            vals[n,0] = float(fdata[0])
                        else:
                            vals[n,0] = missing_value
                            input_missing[n,0] = 1
                        tempstr = fdata[1]
                        char_flag = 0
                        for m in range(len(tempstr)):
                            char_test = tempstr[m].isalpha()
                            if char_test:	#Alphabet character found
                                char_flag = char_flag + 1
                        if fdata[1][0] == '+' and char_flag == 0: 
                            vals[n,1] = float(fdata[1])
                        else:
                            vals[n,1] = missing_value
                            input_missing[n,1] = 1
                    else:	#No relevant data in line or data missing
                        input_missing[n,0] = 1
                        input_missing[n,1] = 1

                    n += 1

            f.close()

        print('no. values  = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        input_missing = input_missing[0:n,:]
        timesecs   = timesecs[0:n]

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]	#Assumes evenly spaced, may want to improve this

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        #Make any qc_flag values which were missing at input == 2. These won't have been read from the correction file
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)
        #qc_flag = qc_flag + input_missing	#Add value of 1 from any values that were set as missing when reading data
        #Can't use valid_min_max from function as won't have the input_missing values
        #valid_min_max = corr_details[1]


        print('qc_flag shape = ',qc_flag.shape)
        print('valid_min_max = ',valid_min_max)

        #vals_qc = np.where(qc_flag <> 2, vals, 0)
        #vals_qc = np.where(qc_flag <> 3, vals_qc, missing_value)
        #vals = np.where(qc_flag <> 3, vals, missing_value)
        vals = np.where(qc_flag != 2, vals, missing_value)


        #print(qc_flag[0:49,0])


    #--------------------------------------------
    # Sorting out day/time information
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        #lengths_of_dimensions = time_details[8]
        #substitutions = time_details[9]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])
        print(first_last_datetime)

        print('vals shape = ',vals.shape)

        # -------------------------------
        # Open new netCDF file for output
        # -------------------------------
        #cfarr_head = 'amof-pluvio_chilbolton_'
        cfarr_head = 'ncas-rain-gauge-4_cao_'
        #out_file   = os.path.join(path_out,cfarr_head+datestring+'.nc')
        out_file   = path_out+datestring+data_product+data_version+'.nc'
        print('Opening new NetCDF file ' + out_file)

        lengths_of_dimensions = {"time": n}
        print(epoch_timesecs[0])
        print(epoch_timesecs[n-1])
        print(datetime.datetime.now())
        print(datetime.datetime(year,month,day,0,0,0))
        print(datetime.datetime(year,month,day,23,59,59))
        print(datetime.date(year,month,day))

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr), "end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec": np.amax(sc_arr), "min_rr": valid_min_max[0,0], "max_rr": valid_min_max[1,0], "min_thick": valid_min_max[0,1], "max_thick": valid_min_max[1,1], "uname": user_id, "codename": __file__, "machinename": computer_id}

        data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
        data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
        data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]      #Need to include fraction of day
        data_object["variables"]["year"]["values"][:] = yr_arr[:]
        data_object["variables"]["month"]["values"][:] = mn_arr[:]
        data_object["variables"]["day"]["values"][:] = dy_arr[:]
        data_object["variables"]["hour"]["values"][:] = hr_arr[:]
        data_object["variables"]["minute"]["values"][:] = mi_arr[:]
        data_object["variables"]["second"]["values"][:] = sc_arr[:]
        data_object["variables"]["rainfall_rate"]["values"][:] = vals[:,0]
        data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,1]
        data_object["variables"]["qc_flag_rainfall_rate"]["values"][:] = qc_flag[:,0]
        data_object["variables"]["qc_flag_thickness_of_rainfall_amount"]["values"][:] = qc_flag[:,1]
        exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
        print(exit_code)



        variables_to_plot = ["rainfall_rate", "thickness_of_rainfall_amount"]
        generate_netcdf_graphs(out_file, graph_path_out, variables_to_plot, datestring)



# ------------------------------------------------------------------------
# Chilbolton Campbell datalogger netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_met(nday,num_today):


    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['oatnew_ch', 'rhnew_ch', 'QFE_ch','ws_ch', 'wd_ch', 'rg001dc_ch', 'rg006dc_ch', 'rg008dc_ch', 'rg004tb_ch', 'rg009dc_ch']
    nvar = np.size(chids)
    print('nvar = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    #wscal = 0.01
    #wdcal = 0.045
    #wdoff = 180 
    #tcal = 0.02
    #toff = -40.0
    #rhcal = 0.02
    #stdrgcal = 0.0033	#Nominal cal. for standard rate raingauge
    daily_file_read = 0	#Flag showing whether we are reading data from an ongoing logger file (= 0) or an already-written daily file (= 1)
    daily_file_exist = 0	#Flag showing whether a daily file already exists(= 1)
    format_threshold1 = 737635   #Date when added raingauges to CR1000X = 20200729
    format_threshold2 = 737867   #Date when real data from gauge near flux compound (009) = 20210318
    if nday < format_threshold1:
        nn_max = 3	#Just met sensors, no raingauges, before 20200729
    elif nday >= format_threshold1 and nday < format_threshold2:
        nn_max = 8      #Met sensors, raingauge and multiple raingauges at Chilbolton
    else:
        nn_max = 9

 
    missing_value = -1.0E+20

    #---------------------------------------
    # Read purge times from hmp155_purge.txt
    #---------------------------------------

    #purge_limit = load_hmp155_purgetime(nday)	


    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    #timesecs    = np.zeros((10000))
    #vals        = missing_value*np.ones((10000,nvar))
    input_missing =  np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,(nvar+1)))	#Making it nvar+1 is a fudge to make substitutions line work for writing the diagnostics T/RH file

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_grape_loggernet/"
    path_out_1 = "/data/amof-netCDF/ncas-temperature-rh-1/ncas-temperature-rh-1_cao_"
    path_out_2 = "/data/amof-netCDF/ncas-pressure-1/ncas-pressure-1_cao_"
    path_out_3 = "/data/amof-netCDF/ncas-anemometer-2/ncas-anemometer-2_cao_"
    path_out_4 = "/data/amof-netCDF/ncas-rain-gauge-1/ncas-rain-gauge-1_cao_"
    path_out_5 = "/data/amof-netCDF/ncas-rain-gauge-2/ncas-rain-gauge-2_cao_"
    path_out_6 = "/data/amof-netCDF/ncas-rain-gauge-3/ncas-rain-gauge-3_cao_"
    path_out_7 = "/data/amof-netCDF/ncas-rain-gauge-5/ncas-rain-gauge-5_cao_"
    path_out_8 = "/data/amof-netCDF/ncas-rain-gauge-9/ncas-rain-gauge-9_cao_"
    path_out_9 = "/data/amof-netCDF/diagnostics/ncas-nocorr-temperature-rh-1/ncas-nocorr-temperature-rh-1_cao_"
    daily_path_out = "/data/range/daily_met/cr1000x_rxcabin_1/"
    daily_file = 'CR1000XSeries_Chilbolton_Rxcabinmet1_'
    daily_file   = os.path.join(daily_path_out,daily_file+datestring+'.dat')

    template_file_path_1 = "/home/jla/python/metsensors/ncas-temperature-rh-1_metadata-template.yaml"	#YAML template
    template_file_path_2 = "/home/jla/python/metsensors/ncas-pressure-1_metadata-template.yaml"	#YAML template
    template_file_path_3 = "/home/jla/python/metsensors/ncas-anemometer-2_metadata-template.yaml"	#YAML template
    template_file_path_4 = "/home/jla/python/metsensors/ncas-rain-gauge-1_metadata-template.yaml"	#YAML template rg001dc_ch
    template_file_path_5 = "/home/jla/python/metsensors/ncas-rain-gauge-2_metadata-template.yaml"	#YAML template rg006dc_ch
    template_file_path_6 = "/home/jla/python/metsensors/ncas-rain-gauge-3_metadata-template.yaml"	#YAML template rg008dc_ch
    template_file_path_7 = "/home/jla/python/metsensors/ncas-rain-gauge-5_metadata-template.yaml"	#YAML template rg004tb_ch
    template_file_path_8 = "/home/jla/python/metsensors/ncas-rain-gauge-9_metadata-template.yaml"	#YAML template rg009dc_ch

    tfp = [template_file_path_1, template_file_path_2, template_file_path_3, template_file_path_4, template_file_path_5, template_file_path_6, template_file_path_7, template_file_path_8, template_file_path_1 ]
    path_out = [path_out_1, path_out_2, path_out_3, path_out_4, path_out_5, path_out_6, path_out_7, path_out_8, path_out_9 ]
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-temperature-rh-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_2 = "/data/amof-netCDF/graphs/ncas-pressure-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_3 = "/data/amof-netCDF/graphs/ncas-anemometer-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_4 = "/data/amof-netCDF/graphs/ncas-rain-gauge-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_5 = "/data/amof-netCDF/graphs/ncas-rain-gauge-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_6 = "/data/amof-netCDF/graphs/ncas-rain-gauge-3" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_7 = "/data/amof-netCDF/graphs/ncas-rain-gauge-5" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_8 = "/data/amof-netCDF/graphs/ncas-rain-gauge-9" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1, graph_path_out_2, graph_path_out_3, graph_path_out_4, graph_path_out_5, graph_path_out_6, graph_path_out_7, graph_path_out_8]

    #data_product = "_surface-met"

    infiles = "CR1000XSeries_Chilbolton_Rxcabinmet1.dat"
    #infiles = "CR1000XSeries_Chilbolton_Rxcabinmet1_20210118.dat"
    path_in_file = path_in + infiles
    #infiles = "CR1000XSeries_Chilbolton2_Rxcabinmet_array.dat"

    # ---------------------------------------------------
    # Process data files
    # Read the cumulative file from the datalogger first
    # ---------------------------------------------------
    print('Starting to process files')

    #Check if a daily file already exists. If it does, we want to avoid overwriting it by opening it for writing
    if os.path.isfile(daily_file):	#A daily file exists
        daily_file_exist = 1


    file_vals = read_cr1000x_chilbolton1(path_in_file, daily_file, daily_file_read, daily_file_exist, nday, num_today,nvar)
    n = file_vals[0]
    timesecs = file_vals[1]
    vals = file_vals[2]

    print('no. values  = ', n)

    # ---------------------------------------------------
    # If you don't find any data from that day in the cumulative file, 
    # look for a daily file instead.
    # This case is most likely if you are processing older data which have
    # been removed from the cumulative file to reduce its size.
    # In this case, look for a previously written daily file
    # and don't try to write a daily file again. 
    # ---------------------------------------------------

    if n == 0: #No data found so far
        if os.path.isfile(daily_file):	#A daily file exists
            print('Reading from a daily file')
            daily_file_read = 1
            file_vals = read_cr1000x_chilbolton1(daily_file, daily_file, daily_file_read, daily_file_exist, nday, num_today, nvar)
            n = file_vals[0]
            timesecs = file_vals[1]
            vals = file_vals[2]

            print('no. values from daily file = ', n)




    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)


        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])
        print(vals[0:9,2])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]


        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        #qc_flag = qc_flag + input_missing	#Add value of 1 from any values that were set as missing when reading data

        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        #valid_min_max = corr_details[1]
        print('qc_flag shape = ',qc_flag.shape)
        print('valid_min_max = ',valid_min_max)
        print('qc_flag shape = ',qc_flag.shape)

        raw_vals = vals	#Keep a copy of uncorrected values
        #Set to missing value if correction file showed BADDATA
        vals = np.where(qc_flag != 2, vals, missing_value)


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])
        print(first_last_datetime)

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]


        #substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "T_min": valid_min_max[0,0], "T_max": valid_min_max[1,0], "RH_min": valid_min_max[0,1], "RH_max": valid_min_max[1,1], "P_min": valid_min_max[0,2], "P_max": valid_min_max[1,2], "WS_min": valid_min_max[0,3], "WS_max": valid_min_max[1,3], "WD_min": valid_min_max[0,4], "WD_max": valid_min_max[1,4], "min_thick": np.amin(vals[:,5]), "max_thick": np.amax(vals[:,5]), "uname": user_id, "codename": __file__, "machinename": computer_id}

        for nn in range(nn_max):	#Each data file

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "T_min": valid_min_max[0,0], "T_max": valid_min_max[1,0], "RH_min": valid_min_max[0,1], "RH_max": valid_min_max[1,1], "P_min": valid_min_max[0,2], "P_max": valid_min_max[1,2], "WS_min": valid_min_max[0,3], "WS_max": valid_min_max[1,3], "WD_min": valid_min_max[0,4], "WD_max": valid_min_max[1,4], "min_thick": valid_min_max[0,(nn+2)], "max_thick": valid_min_max[1,(nn+2)], "uname": user_id, "codename": __file__, "machinename": computer_id}
            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            if nn >= 3 and nn <= 7:
                data_product = "_precipitation"
            else: 
                data_product = "_surface-met"

            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]

            if nn == 0:
                chids = ["oatnew_ch", "rhnew_ch"]
                data_object["variables"]["air_temperature"]["values"][:] = vals[:,0]
                data_object["variables"]["relative_humidity"]["values"][:] = vals[:,1]
                data_object["variables"]["qc_flag_air_temperature"]["values"][:] = qc_flag[:,0]
                data_object["variables"]["qc_flag_relative_humidity"]["values"][:] = qc_flag[:,1]
                variables_to_plot = ["air_temperature", "relative_humidity"]
            if nn == 1:
                chids = ["QFE_ch"]
                data_object["variables"]["air_pressure"]["values"][:] = vals[:,2] 
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,2]
                variables_to_plot = ["air_pressure"]
            if nn == 2:
                chids = ["ws_ch", "wd_ch"]
                data_object["variables"]["wind_speed"]["values"][:] = vals[:,3] 
                data_object["variables"]["wind_from_direction"]["values"][:] = vals[:,4] 
                data_object["variables"]["qc_flag_wind_speed"]["values"][:] = qc_flag[:,3]
                data_object["variables"]["qc_flag_wind_from_direction"]["values"][:] = qc_flag[:,4]
                variables_to_plot = ["wind_speed", "wind_from_direction"]
            if nn == 3:
                chids = ["rg001dc_ch"]
                data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,5] 
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,5]
                variables_to_plot = ["thickness_of_rainfall_amount"]
            if nn == 4:
                chids = ["rg006dc_ch"]
                data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,6] 
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,6]
                variables_to_plot = ["thickness_of_rainfall_amount"]
            if nn == 5:
                chids = ["rg008dc_ch"]
                data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,7] 
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,7]
                variables_to_plot = ["thickness_of_rainfall_amount"]
            if nn == 6:
                chids = ["rg004tb_ch"]
                data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,8] 
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,8]
                variables_to_plot = ["thickness_of_rainfall_amount"]
            if nn == 7:
                chids = ["rg009dc_ch"]
                data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,9] 
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,9]
                variables_to_plot = ["thickness_of_rainfall_amount"]
            if nn == 8:
                chids = ["oatnew_ch", "rhnew_ch"]
                data_object["variables"]["air_temperature"]["values"][:] = raw_vals[:,0]
                data_object["variables"]["relative_humidity"]["values"][:] = raw_vals[:,1]
                data_object["variables"]["qc_flag_air_temperature"]["values"][:] = qc_flag[:,0]
                data_object["variables"]["qc_flag_relative_humidity"]["values"][:] = qc_flag[:,1]
            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print(exit_code)

            if nn < 8:	#Don't plot T/RH file with no corrections to data - it's just for quality control
                generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)


# ------------------------------------------------------------------------
# Read data from Chilbolton1 Campbell datalogger
# ------------------------------------------------------------------------
def read_cr1000x_chilbolton1(in_file, write_file, daily_file_read, daily_file_exist, nday, num_today, nvar):

    # -------------------------
    # Define various parameters
    # -------------------------
    wscal = 0.01
    wdcal = 0.045
    wdoff = 180 
    tcal = 0.02
    toff = -40.0
    rhcal = 0.02
    stdrgcal = 0.0033	#Nominal cal. for standard rate raingauge
    lowrgcal = 0.00189	#Nominal cal. for low rate raingauge
    missing_value = -1.0E+20
    n = 0	#Times
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nvar))
    write_daily_ok = 0	#Variable which shows whether we write a daily .dat file, 1 if we do write one

    f = open(in_file, 'r')
    if daily_file_read == 0 and daily_file_exist == 0 and (num_today - nday) >= 1:
        write_daily_ok = 1
        fd = open(write_file, 'w')


    for z in range(4):      #4 header lines, may eventually want some of them
        line = f.readline()
        if daily_file_read == 0 and daily_file_exist == 0:	#Keep the header lines to write out later unless you're reading a daily file
            if z == 0:
                header_list = [line]
            else:
                header_list.append(line)

    # ---------------------------
    # Repeated data for each time
    # ---------------------------

    while True:

        line = f.readline()
        if not line: break
        #startnum = int(date2num(datetime.datetime(int(line[1:5]),int(line[6:8]),int(line[9:11]),0,0,0)))
        startnum = date2num(datetime.datetime(int(line[1:5]),int(line[6:8]),int(line[9:11]),int(line[12:14]),int(line[15:17]),int(line[18:20])))
        if startnum > nday and startnum <= (nday + 1):
            #if n == 0 and daily_file_read == 0 and daily_file_exist == 0:	#Write header to file if this is first line for today's date
            if n == 0 and write_daily_ok == 1:	#Write header to file if this is first line for today's date
                for z in range(4):
                    fd.write(header_list[z])	    

            timesecs[n] = 3600.0*float(line[12:14]) + 60.0*float(line[15:17]) + float(line[18:20])
            if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                timesecs[n] = timesecs[n] + 86400.0
 
            fdata = line.split(',')
            #If a data value is missing, the string has 5 characters, "NAN". Not all values will create such a long string, so just look for the first N
            if fdata[6][1:2] != "N":	#Was 4
                vals[n,2] = float(fdata[6])	#Pressure (Pa)	#was n, not m
            if fdata[7][1:2] != "N":	#Was 5
                vals[n,3] = wscal * float(fdata[7])	#Wind speed (mV)
            if fdata[8][1:2] != "N":	#Was 6
                if np.absolute(float(fdata[8])) <= 5000.0:	#If it's a number, then check its magnitude is <=5V, to be valid
                    vals[n,4] = wdcal * float(fdata[8]) + wdoff	#Wind direction (mV)
                    if vals[n,4] >= 360.0:
                        vals[n,4] = vals[n,4] - 360.0
                    if vals[n,4] < 0.0:
                        vals[n,4] = vals[n,4] + 360.0
                #else:	#Data out of range. These 4 lines taken out to test assigning input_missing values from missing_value in vals
                    #input_missing[n,4] = 1
            #else:
                #input_missing[n,4] = 1
            #if timesecs[n] < purge_limit[0] or timesecs[n] > purge_limit[1]:	#Only use this if want to leave these values as missing
            if fdata[4][1:2] != "N":
                vals[n,0] = tcal*float(fdata[4]) + toff + 273.15 #Temperature
            if fdata[5][1:2] != "N":
                vals[n,1] = rhcal * float(fdata[5]) #RH
            if fdata[10][1:2] != "N":
                vals[n,5] = stdrgcal * float(fdata[10]) #Drop counts
            if len(fdata) > 11: #Prior to partway through 20200729 these 3 columns weren't present
                if fdata[11][1:2] != "N":
                    vals[n,6] = stdrgcal * float(fdata[11]) #Drop counts
                if fdata[12][1:2] != "N":
                    vals[n,7] = lowrgcal * float(fdata[12]) #Drop counts
                if startnum > 737824.46:	#03/02/21 11:05, when channel was added for rg009dc_ch
                    if fdata[14][1:2] != "N":
                        vals[n,8] = float(fdata[14]) #Tipping bucket
                else:
                    if fdata[13][1:2] != "N":
                        vals[n,8] = float(fdata[13]) #Tipping bucket
                if startnum >= 737867.0:
                    if fdata[13][1:2] != "N":
                        vals[n,9] = stdrgcal * float(fdata[13]) #rg009dc_ch
            #Write line to daily file
            if write_daily_ok == 1:
                fd.write(line)
            n += 1

    f.close()
    if write_daily_ok == 1:
        fd.close()	#Writes daily files
    
    print('no. values  = ', n)

    return n, timesecs, vals


# ------------------------------------------------------------------------
# Chilbolton Campbell datalogger netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_sparsholt_raingauge(nday,num_today):


    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['rg001dc_sp', 'rg002tb_sp']
    nvar = np.size(chids)
    print('nvar = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    #stdrgcal = 0.0033	#Nominal cal. for standard rate raingauge
    #tbrgcal = 0.2	#Nominal cal. for tipping bucket raingauge 
    daily_file_read = 0	#Flag showing whether we are reading data from an ongoing logger file (= 0) or an already-written daily file (= 1)
    daily_file_exist = 0	#Flag showing whether a daily file already exists(= 1)
 
    missing_value = -1.0E+20

    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    #timesecs    = np.zeros((10000))
    #vals        = missing_value*np.ones((10000,nvar))
    input_missing =  np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    path_in = "/data/range/mirror_grape_loggernet/"
    path_out_1 = "/data/amof-netCDF/ncas-rain-gauge-7/ncas-rain-gauge-7_cao_"
    path_out_2 = "/data/amof-netCDF/ncas-rain-gauge-6/ncas-rain-gauge-6_cao_"

    template_file_path_1 = "/home/jla/python/metsensors/ncas-rain-gauge-7_metadata-template.yaml"	#YAML template
    template_file_path_2 = "/home/jla/python/metsensors/ncas-rain-gauge-6_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1, template_file_path_2]
    path_out = [path_out_1, path_out_2]
    daily_path_out = "/data/range/daily_met/cr1000x_sparsholt_1/"
    daily_file = 'CR1000XSeries_Sparsholt1_'
    daily_file   = os.path.join(daily_path_out,daily_file+datestring+'.dat')
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-rain-gauge-7" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_2 = "/data/amof-netCDF/graphs/ncas-rain-gauge-6" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1, graph_path_out_2]
    data_product = "_precipitation"

    infiles = "CR1000XSeries_Sparsholt_dc_tb.dat"
    #infiles = "CR1000XSeries_Sparsholt_dc_tb_20200623.dat"
    path_in_file = path_in + infiles

    # ---------------------------------------------------
    # Process data files
    # Read the cumulative file from the datalogger first
    # ---------------------------------------------------
    print('Starting to process files')

    #Check if a daily file already exists. If it does, we want to avoid overwriting it by opening it for writing
    if os.path.isfile(daily_file):	#A daily file exists
        daily_file_exist = 1


    file_vals = read_cr1000x_sparsholt(path_in_file, daily_file, daily_file_read, daily_file_exist, nday, num_today, nvar)
    n = file_vals[0]
    timesecs = file_vals[1]
    vals = file_vals[2]
    
    print('no. values  = ', n)

    # ---------------------------------------------------
    # If you don't find any data from that day in the cumulative file, 
    # look for a daily file instead.
    # This case is most likely if you are processing older data which have
    # been removed from the cumulative file to reduce its size.
    # In this case, look for a previously written daily file
    # and don't try to write a daily file again. 
    # ---------------------------------------------------

    if n == 0: #No data found so far
        if os.path.isfile(daily_file):	#A daily file exists
            print('Reading from a daily file')
            daily_file_read = 1
            file_vals = read_cr1000x_sparsholt(daily_file, daily_file, daily_file_read, daily_file_exist, nday, num_today, nvar)
            n = file_vals[0]
            timesecs = file_vals[1]
            vals = file_vals[2]

            print('no. values from daily file = ', n)

    print('no. values  = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 3
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        #valid_min_max = corr_details[1]
        print('qc_flag shape = ',qc_flag.shape)
        print('valid_min_max = ',valid_min_max)

        print(qc_flag[0:49,0])

        vals = np.where(qc_flag != 2, vals, missing_value)


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        for nn in range(nvar):	#Each data file

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "min_thick": valid_min_max[0,nn], "max_thick": valid_min_max[1,nn], "uname": user_id, "codename": __file__, "machinename": computer_id}


            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]

            data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,nn] 
            data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,nn]
            variables_to_plot = ["thickness_of_rainfall_amount"]
            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print(exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)


# ------------------------------------------------------------------------
# Read data from Sparsholt Campbell datalogger
# ------------------------------------------------------------------------
def read_cr1000x_sparsholt(in_file, write_file, daily_file_read, daily_file_exist, nday, num_today, nvar):

    # -------------------------
    # Define various parameters
    # -------------------------
    stdrgcal = 0.0033	#Nominal cal. for standard rate raingauge
    tbrgcal = 0.2	#Nominal cal. for tipping bucket raingauge 
    missing_value = -1.0E+20
    n = 0	#Times
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nvar))
    write_daily_ok = 0	#Variable which shows whether we write a daily .dat file, 1 if we do write one

    f = open(in_file, 'r')
    if daily_file_read == 0 and daily_file_exist == 0 and (num_today - nday) >= 1:
        write_daily_ok = 1
        fd = open(write_file, 'w')


    for z in range(4):      #4 header lines, may eventually want some of them
        line = f.readline()
        if daily_file_read == 0 and daily_file_exist == 0:	#Keep the header lines to write out later unless you're reading a daily file
            if z == 0:
                header_list = [line]
            else:
                header_list.append(line)

    # ---------------------------
    # Repeated data for each time
    # ---------------------------

    while True:

        line = f.readline()
        if not line: break
        #startnum = int(date2num(datetime.datetime(int(line[1:5]),int(line[6:8]),int(line[9:11]),0,0,0)))
        startnum = date2num(datetime.datetime(int(line[1:5]),int(line[6:8]),int(line[9:11]),int(line[12:14]),int(line[15:17]),int(line[18:20])))
        #if startnum == nday:
        if startnum > nday and startnum <= (nday + 1):
            #if n == 0 and daily_file_read == 0 and daily_file_exist == 0:	#Write header to file if this is first line for today's date
            if n == 0 and write_daily_ok == 1:	#Write header to file if this is first line for today's date
                for z in range(4):
                    fd.write(header_list[z])	    

            timesecs[n] = 3600.0*float(line[12:14]) + 60.0*float(line[15:17]) + float(line[18:20])
            if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                timesecs[n] = timesecs[n] + 86400.0
            fdata = line.split(',')
            #If a data value is missing, the string has 5 characters, "NAN". Not all values will create such a long string, so just look for the first N
            if fdata[5][1:2] != "N":
                vals[n,0] = stdrgcal * float(fdata[5])	#rg001dc_sp
            if fdata[6][1:2] != "N":
                vals[n,1] = float(fdata[6])	#rg002tb_sp	#Already calibrated in logger file
            #Write line to daily file
            if write_daily_ok == 1:
                fd.write(line)
            n += 1

    f.close()
    if write_daily_ok == 1:
        fd.close()	#Writes daily files

    print('no. values  = ', n)

    return n, timesecs, vals



# ------------------------------------------------------------------------
# Raingauge from chan* file netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_rain_f5(nday):

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['rg001dc_ch', 'rg006dc_ch', 'rg008dc_ch', 'rg004tb_ch']	#rain-gauge-1, rain-gauge-2, rain-gauge-3, rain-gauge-5
    nvar = np.size(chids)
    print('nvar = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    stdrgcal = 0.0033	#Nominal cal. for standard rate raingauge
    lowrgcal = 0.00189	#Nominal cal. for low rate raingauge 
    tbrgcal = 0.2	#Nominal cal. for tipping bucket raingauge 
    missing_value = -1.0E+20

    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nvar))
    input_missing    = np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))


    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    #path_in = "/home/jla/loggernet/"
    path_in = "/data/range/mirror_marvin_home2/ranged/analog/"
    path_out_1 = "/data/amof-netCDF/ncas-rain-gauge-1/ncas-rain-gauge-1_cao_"
    path_out_2 = "/data/amof-netCDF/ncas-rain-gauge-2/ncas-rain-gauge-2_cao_"
    path_out_3 = "/data/amof-netCDF/ncas-rain-gauge-3/ncas-rain-gauge-3_cao_"
    path_out_4 = "/data/amof-netCDF/ncas-rain-gauge-5/ncas-rain-gauge-5_cao_"

    template_file_path_1 = "/home/jla/python/metsensors/ncas-rain-gauge-1_ulink_metadata-template.yaml"	#YAML template
    template_file_path_2 = "/home/jla/python/metsensors/ncas-rain-gauge-2_ulink_metadata-template.yaml"	#YAML template
    template_file_path_3 = "/home/jla/python/metsensors/ncas-rain-gauge-3_ulink_metadata-template.yaml"	#YAML template
    template_file_path_4 = "/home/jla/python/metsensors/ncas-rain-gauge-5_ulink_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1, template_file_path_2, template_file_path_3, template_file_path_4]
    path_out = [path_out_1, path_out_2, path_out_3, path_out_4]
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-rain-gauge-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_2 = "/data/amof-netCDF/graphs/ncas-rain-gauge-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_3 = "/data/amof-netCDF/graphs/ncas-rain-gauge-3" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_4 = "/data/amof-netCDF/graphs/ncas-rain-gauge-5" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1, graph_path_out_2, graph_path_out_3, graph_path_out_4]
    data_product = "_precipitation"

    # -------------------------------------------------------------------------
    # Loop to process raw data from previous, current and next day's data files
    # Data for the current day could be found in each of these format5 files.
    # contains data entirely from the previous day.
    # Hence, need to open files from all 3 days.
    # -------------------------------------------------------------------------

    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to process files')
    n = 0       #Times

    for day_incr in range(3):

        nday_file = nday + day_incr - 1
        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        print(nday_file, datevals)

        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "chan" + datestring_now[2:] + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print(infiles)
            nfiles = len(infiles)
            if nfiles == 0:
                print('No files on this day')

        else:
            nfiles = 0
            print('No directory found')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            for z in range(17):      #4 header lines, may eventually want some of them
                line = f.readline()

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:

                line = f.readline()
                if not line: break
                #May need to change to handle change of year. Or maybe easier to just lose 1 point once a year?! Stems from not having year in the f5 file
                startnum = date2num(datetime.datetime(int(year_now),int(line[0:2]),int(line[3:5]),int(line[6:8]),int(line[9:11]),int(line[12:14])))
                if startnum > nday and startnum <= (nday + 1):

                #startnum = int(date2num(datetime.datetime(int(datestring[0:4]),int(line[0:2]),int(line[3:5]),0,0,0)))
                #if startnum == nday:

                    timesecs[n] = 3600.0*float(line[6:8]) + 60.0*float(line[9:11]) + float(line[12:14])
                    if n == 0:
                        print('First point = ', nday, timesecs[0], startnum, datestring_now)
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0

                    fdata = line.split()	#Need to not specify ' ' or it treats sucessive delimiters as individuals
                    vals[n,1] = stdrgcal * int(fdata[1])	#rg006dc_ch
                    vals[n,2] = lowrgcal * int(fdata[2])	#rg008dc_ch
                    vals[n,3] = tbrgcal * int(fdata[4])	#rg004tb_ch
                    vals[n,0] = stdrgcal * int(fdata[15])	#rg001dc_ch
                    n += 1

            f.close()

        print('no. values  = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])
        print(vals[0:9,2])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs,vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 3
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        valid_min_max = corr_details[1]
        print('qc_flag shape = ',qc_flag.shape)

        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)	#Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 3, vals, missing_value)


        print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values 


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        #substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],sc_arr[0]), "end_time": datetime.datetime(year,month,day,hr_arr[n-1],mi_arr[n-1],sc_arr[n-1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "T_min": np.amin(vals[:,0]), "T_max": np.amax(vals[:,0]), "RH_min": np.amin(vals[:,1]), "RH_max": np.amax(vals[:,1]), "P_min": np.amin(vals[:,2]), "P_max": np.amax(vals[:,2]), "WS_min": np.amin(vals[:,3]), "WS_max": np.amax(vals[:,3]), "WD_min": np.amin(vals[:,4]), "WD_max": np.amax(vals[:,4]), "min_thick": np.amin(vals[:,5]), "max_thick": np.amax(vals[:,5]), "uname": user_id, "codename": __file__, "machinename": computer_id}


        for nn in range(nvar):	#Each data file

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "min_thick": valid_min_max[0,nn], "max_thick": valid_min_max[1,nn], "uname": user_id, "codename": __file__, "machinename": computer_id}

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]


            #chids = ["rg001dc_ch"]
            data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = vals[:,nn] 
            data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,nn]
            variables_to_plot = ["thickness_of_rainfall_amount"]

            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print(exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)


# ------------------------------------------------------------------------
# RD-80 impact disdrometer chds* or spds* file netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_disdro_f5(nday,sensor):


    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals)
    if sensor == 6:
        chids = ['disdrom_ch']  #Chilbolton disdrometer
        tempstr = 'chds' + datestring[2:] + '.00'	#Root of filename to search for in directory
    else:
        chids = ['disdrom_sp']	#Sparsholt disdrometer 
        tempstr = 'spds' + datestring[2:] + '.00'
    nvar = np.size(chids)
    print('nvar = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    missing_value = -1.0E+20
    nbins = 127
    disdro_A=5.0265E-03 #m2
    #Minimum diameters for each bin
    rd80_thresholds=[0.313, 0.318, 0.324, 0.331,
                    0.337, 0.343, 0.35, 0.357, 0.364, 0.371,
                    0.379, 0.387, 0.396, 0.405, 0.414, 0.423,
                    0.433, 0.443, 0.453, 0.464, 0.474, 0.484,
                    0.495, 0.505, 0.515, 0.526, 0.537, 0.547,
                    0.558, 0.57, 0.583, 0.596, 0.611, 0.628,
                    0.644, 0.662, 0.679, 0.696, 0.715, 0.735,
                    0.754, 0.771, 0.787, 0.806, 0.827, 0.845,
                    0.862, 0.879, 0.895, 0.912, 0.928, 0.944,
                    0.96, 0.978, 0.999, 1.024, 1.051, 1.08,
                    1.11, 1.14, 1.171, 1.202, 1.232, 1.262,
                    1.289, 1.318, 1.346, 1.374, 1.402, 1.429,
                    1.456, 1.483, 1.509, 1.533, 1.558, 1.582,
                    1.606, 1.631, 1.657, 1.683, 1.715, 1.748,
                    1.793, 1.841, 1.897, 1.955, 2.013, 2.077,
                    2.139, 2.2, 2.262, 2.321, 2.381, 2.441,
                    2.499, 2.558, 2.616, 2.672, 2.727, 2.781,
                    2.836, 2.893, 2.949, 3.011, 3.08, 3.155,
                    3.23, 3.306, 3.385, 3.466, 3.545, 3.625,
                    3.704, 3.784, 3.864, 3.945, 4.028, 4.127,
                    4.231, 4.34, 4.456, 4.573, 4.686, 4.801,
                    4.915, 5.03, 5.145]



    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nbins))	#Could/should expand to 10000,nbins,nvar
    input_missing    = np.zeros((10000,nvar))
    #valid_min_max = np.zeros((2))
    mean_diam    = np.zeros((127))
    mean_vol     = np.zeros((127))
    ds_bin_width = np.zeros((127))



    # ---------------------------------------------------
    # Calculate mean drop diameter, volume and bin widths 
    # ---------------------------------------------------
    for n in range(126):
        mean_diam[n]=np.mean([rd80_thresholds[n],rd80_thresholds[n+1]])
        mean_vol[n]=0.5236*np.mean([rd80_thresholds[n]**3,rd80_thresholds[n+1]**3])
        ds_bin_width[n]=rd80_thresholds[n+1]-rd80_thresholds[n]
    mean_diam[126] = 5.203      #arbitrary value based on spacing of previous bins
    mean_vol[126] = 74.0
    ds_bin_width[126] = ds_bin_width[125]       #arbitrary value based on spacing of previous bins

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    #path_in = "/home/jla/loggernet/"
    data_version = "_v1.0"
    if sensor == 6:
        path_in = "/data/range/mirror_marvin_home2/ranged/distrom/"
        if nday < 737916:	#20210506
        #if nday < 738000:	#fake date in future
            path_out_1 = "/data/amof-netCDF/ncas-disdrometer-1/ncas-disdrometer-1_cao_"
            template_file_path_1 = "/home/jla/python/metsensors/ncas-disdrometer-1_metadata-template.yaml"	#YAML template
            graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-disdrometer-1" + data_version + "/" + str(year) + '/' + strmon + '/'
        else:
            path_out_1 = "/data/amof-netCDF/ncas-disdrometer-2/ncas-disdrometer-2_cao_"
            template_file_path_1 = "/home/jla/python/metsensors/ncas-disdrometer-2-chilbolton_metadata-template.yaml"	#YAML template
            graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-disdrometer-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    else:	#Should only be 6 or 7 if you get to this function
        path_in = "/data/sparsholt/mirror_rhubarb_home2/ranged/distrom/"
        path_out_1 = "/data/amof-netCDF/ncas-disdrometer-2-sparsholt/ncas-disdrometer-2_cao-sparsholt_"
        graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-disdrometer-2-sparsholt" + data_version + "/" + str(year) + '/' + strmon + '/'
        template_file_path_1 = "/home/jla/python/metsensors/ncas-disdrometer-2_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1]
    path_out = [path_out_1]
    graph_path_out = [graph_path_out_1]
    data_product = "_precipitation"

    if os.access(path_in, os.F_OK):     #A data directory exists for the day
        files = os.listdir(path_in)
        expr = re.compile(tempstr+'[0-9]')
        infiles = [elem for elem in files if expr.search(elem)]
        infiles.sort()
        print(infiles)
        nfiles = len(infiles)
        if nfiles == 0:
            print('No files on this day')

    else:
        nfiles = 0
        print('No directory found')


    #nfiles = 1
    #infiles = path_in + "chan" + datestring[2:] + ".000"	#Might need to add .001 etc. files

    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to process files')
    n = 0	#Times

    for nf in range(nfiles):    #Indent from here

        source_file_path = path_in + infiles[nf]
        disdro_handler = module_distrometer_format5.Handler(source_file_path)

        number_of_records = disdro_handler.return_number_of_records()
        extracted_spectral_data = disdro_handler.return_spectral_data(0)       #0 is the index. An index is needed
        #For some reason, when I read values into vals, the value of extracted_spectral_data is set to zero in the print statement below.
        extracted_raw_data = disdro_handler.return_raw_data(0)
        datetimes = disdro_handler.return_datetimes_for_records()
        numtimes = np.rint((date2num(datetimes) - nday) * 86400.0)	#Time in seconds since midnight 
        timesecs[n : (n + number_of_records)] = numtimes 

        for nn in range(number_of_records):
            vals[n,:] = disdro_handler.return_spectral_data(nn)
            n = n + 1	#Counter for total number of points across all files for the days



        print('Number of records = ', number_of_records)
        print('Spectral data = ', extracted_spectral_data)
        print('Raw data = ', extracted_raw_data)
        print('Datetimes[0] = ', datetimes[0])
        print('Numeric time[0] = ', date2num(datetimes[0]))
        print('Numeric times = ', date2num(datetimes))
        print(timesecs)
        print(vals[0,:])
        print(vals[1,:])

        print('n = ', n)
        print('timesecs[n-1] = ',timesecs[n-1])

        if timesecs[n-1] >= 86400.0: 
            n = n - 1


    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])
        print(vals[0:9,2])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        ch_count_vol_product = mean_vol*vals
        print('ch_count_vol_product min, max = ', np.amin(ch_count_vol_product), np.amax(ch_count_vol_product))
        ch_accum_array = (np.sum(ch_count_vol_product,axis=1))/(disdro_A*1.0e+06)
        print( ch_accum_array)
        print( 'ch_accum_array min, max = ', np.amin(ch_accum_array), np.amax(ch_accum_array))
        print( ch_accum_array.shape)

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs,vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 3
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            #if len(vals_ok) > 0:
                #valid_min_max[0,k] = np.amin(vals_ok)
                #valid_min_max[1,k] = np.amax(vals_ok)

        valid_min_max = corr_details[1]
        print('valid min, max = ',valid_min_max)
        print('valid min, max shape = ',valid_min_max.shape)
        valid_min_counts = int(valid_min_max[0])
        valid_max_counts = int(valid_min_max[1])
        print('qc_flag shape = ',qc_flag.shape)
        qc_flag_thickness = qc_flag[:,0]	#1d qc_flag for thickness_of_rainfall_amount
        qc_flag = np.broadcast_to(qc_flag,(n,nbins))
        print('qc_flag shape = ',qc_flag.shape)

        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)	#Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 3, vals, missing_value)


        print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values 


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        #lengths_of_dimensions = time_details[8]
        #substitutions = time_details[9]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        #substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],sc_arr[0]), "end_time": datetime.datetime(year,month,day,hr_arr[n-1],mi_arr[n-1],sc_arr[n-1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "T_min": np.amin(vals[:,0]), "T_max": np.amax(vals[:,0]), "RH_min": np.amin(vals[:,1]), "RH_max": np.amax(vals[:,1]), "P_min": np.amin(vals[:,2]), "P_max": np.amax(vals[:,2]), "WS_min": np.amin(vals[:,3]), "WS_max": np.amax(vals[:,3]), "WD_min": np.amin(vals[:,4]), "WD_max": np.amax(vals[:,4]), "min_thick": np.amin(vals[:,5]), "max_thick": np.amax(vals[:,5]), "uname": user_id, "codename": __file__, "machinename": computer_id}


        for nn in range(nvar):	#Each data file

            substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(year,month,day,hr_arr[0],mi_arr[0],sc_arr[0]), "end_time": datetime.datetime(year,month,day,hr_arr[n-1],mi_arr[n-1],sc_arr[n-1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "min_diam": np.amin(rd80_thresholds), "max_diam": np.amax(rd80_thresholds), "min_counts": valid_min_counts, "max_counts": valid_max_counts, "min_thick": np.amin(ch_accum_array), "max_thick": np.amax(ch_accum_array), "uname": user_id, "codename": __file__, "machinename": computer_id}

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n, "diameter": nbins}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]
            data_object["variables"]["diameter"]["values"][:] = rd80_thresholds[:]

            data_object["variables"]["number_of_hydrometeors_per_size_channel"]["values"][:] = vals[:,:] 
            data_object["variables"]["thickness_of_rainfall_amount"]["values"][:] = ch_accum_array[:] 
            data_object["variables"]["qc_flag_hydrometeors"]["values"][:] = qc_flag[:,:]
            data_object["variables"]["qc_flag_thickness_of_rainfall_amount"]["values"][:] = qc_flag_thickness[:]
            variables_to_plot = ["number_of_hydrometeors_per_size_channel", "thickness_of_rainfall_amount"]
            #variables_to_plot = ["thickness_of_rainfall_amount"]

            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print(exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)


# ------------------------------------------------------------------------
# Broadband radiometer from chpy* file netcdf generation function
# ------------------------------------------------------------------------
def generate_netcdf_bbrad_f5(nday):

    datevals=generate_netcdf_common(nday)
    year=datevals[0]
    month=datevals[1]
    strmon=str(month)
    if month <= 9:
        strmon = '0'+strmon
    day=datevals[2]
    datestring=datevals[3]
    #print "datestring = ",datestring
    print(datevals) 
    chids = ['pyrCM21_ch','pyr_CMP21_ch', 'pyrCG4_ch', 'pyrCG4_ch', 'pyrCP1_ch', 'pyrCP1_T_ch']	#total pyranometer, total pyrgeo, pyrgeo T (duplicate same chids as no separate corr file), diffuse pyran, direct pyrhel, direct pyrhel T
    nvar = np.size(chids)
    print('nvar = ',nvar)

    # -------------------------
    # Define various parameters
    # -------------------------
    A = -4.56582625e-7
    B =  8.97922514e-5
    C = -6.95640241e-3
    D =  2.74163515e-1
    E = -6.23224724
    F = 66.1824897

    A1=0.0010295
    B1=2.391e-4
    G1=1.568e-7
    Tk = 273.15
    SBconst = 5.67e-8
    scal = load_bbrad_calibrations(nday)    
    #scal = np.zeros((nvar-1))
    #scal[0] = 10.89e-6	#CM21
    #scal[1] = 8.7e-6	#CMP21
    #scal[2] = 13.96e-6	#CG4
    #scal[3] = 7.89e-6	#CHP1 
    missing_value = -1.0E+20

    # ----------------------
    # Initialize data arrays
    # This is general so could go in a function, but as it's short, leave it here for now
    # ----------------------
    timesecs    = np.zeros((10000))
    vals        = missing_value*np.ones((10000,nvar))
    input_missing    = np.zeros((10000,nvar))
    valid_min_max = np.zeros((2,nvar))

    # -----------
    # Setup paths
    # -----------
    print('Setting up paths')
    #path_in = "/home/jla/loggernet/"
    path_in = "/data/range/mirror_moe_home2/range/broadband_radiometers/pyr/"
    path_out_1 = "/data/amof-netCDF/ncas-radiometer-1/ncas-radiometer-1_cao_"
    path_out_2 = "/data/amof-netCDF/ncas-radiometer-2/ncas-radiometer-2_cao_"
    path_out_3 = "/data/amof-netCDF/ncas-radiometer-3/ncas-radiometer-3_cao_"
    path_out_4 = "/data/amof-netCDF/ncas-radiometer-4/ncas-radiometer-4_cao_"

    template_file_path_1 = "/home/jla/python/metsensors/ncas-radiometer-1_metadata-template.yaml"	#YAML template
    template_file_path_2 = "/home/jla/python/metsensors/ncas-radiometer-2_metadata-template.yaml"	#YAML template
    template_file_path_3 = "/home/jla/python/metsensors/ncas-radiometer-3_metadata-template.yaml"	#YAML template
    template_file_path_4 = "/home/jla/python/metsensors/ncas-radiometer-4_metadata-template.yaml"	#YAML template

    tfp = [template_file_path_1, template_file_path_2, template_file_path_3, template_file_path_4]
    path_out = [path_out_1, path_out_2, path_out_3, path_out_4]
    data_version = "_v1.0"
    graph_path_out_1 = "/data/amof-netCDF/graphs/ncas-radiometer-1" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_2 = "/data/amof-netCDF/graphs/ncas-radiometer-2" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_3 = "/data/amof-netCDF/graphs/ncas-radiometer-3" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out_4 = "/data/amof-netCDF/graphs/ncas-radiometer-4" + data_version + "/" + str(year) + '/' + strmon + '/'
    graph_path_out = [graph_path_out_1, graph_path_out_2, graph_path_out_3, graph_path_out_4]
    data_product = "_radiation"

    # -------------------------------------------------------------------------
    # Loop to process raw data from previous, current and next day's data files
    # Data for the current day could be found in each of these format5 files.
    # contains data entirely from the previous day.
    # Hence, need to open files from all 3 days.
    # -------------------------------------------------------------------------

    # ----------------------
    # Process data files
    # ----------------------
    print('Starting to process files')
    n = 0       #Times

    for day_incr in range(3):
        n_wrong = -1

        nday_file = nday + day_incr - 1

        # -----------------------------
        # Date-time numbers and strings
        # -----------------------------

        datevals=generate_netcdf_common(nday_file)
        year_now=datevals[0]
        month_now=datevals[1]
        day_now=datevals[2]
        datestring_now=datevals[3]
        #print "datestring_now = ",datestring_now
        print(nday_file, datevals)


        if os.access(path_in, os.F_OK):     #A data directory exists for the day
            files = os.listdir(path_in)
            tempstr = "chpy" + datestring_now[2:] + '.00'
            expr = re.compile(tempstr+'[0-9]')
            infiles = [elem for elem in files if expr.search(elem)]
            infiles.sort()
            print(infiles)
            nfiles = len(infiles)
            print("nfiles = ",nfiles)
            if nfiles == 0:
                print('No files on this day')
            #Check for any files whose name start with . and remove them from list. A rare problem but it happened on 03/03/21
            if nfiles > 0:
                for m in range(nfiles):
                    if infiles[m][0] == ".":
                        n_wrong = m
                        print("File ", n_wrong, " found with wrong format on ", datestring_now)
                        nfiles = nfiles - 1
                        print("Number of files reduced to ", nfiles)
                if n_wrong >= 0:
                    del infiles[n_wrong]
                    print("Modified infiles = ",infiles)


        else:
            nfiles = 0
            print('No directory found')

        for nf in range(nfiles):    #Indent from here
            f = open(path_in+infiles[nf], 'r')

            for z in range(8):      #8 header lines, may eventually want some of them
                line = f.readline()

            # ---------------------------
            # Repeated data for each time
            # ---------------------------

            while True:
 
                line = f.readline()
                if not line: break
                #startnum = int(date2num(datetime.datetime(int(datestring[0:4]),int(line[0:2]),int(line[3:5]),0,0,0)))
                #print 'startnum, nday = ',startnum, nday
                #if startnum == nday:
                #May need to change to handle change of year. Or maybe easier to just lose 1 point once a year?!
                startnum = date2num(datetime.datetime(int(year_now),int(line[0:2]),int(line[3:5]),int(line[6:8]),int(line[9:11]),int(line[12:14])))
                if startnum > nday and startnum <= (nday + 1):

                    timesecs[n] = 3600.0*float(line[6:8]) + 60.0*float(line[9:11]) + float(line[12:14])
                    if n == 0:
                        print('First point = ', nday, timesecs[0], startnum, datestring_now)
                    if startnum == (nday + 1):  #If we're reading the first timestamp from the next day, add 86400 secs to time or it will be zero
                        timesecs[n] = timesecs[n] + 86400.0

                    fdata = line.split()	#Need to not specify ' ' or it treats sucessive delimiters as individuals
                    vals[n,0] = float(fdata[1])/scal[0]	#CM21
                    vals[n,1] = float(fdata[4])/scal[1]	#CMP21
                    vals[n,2] = float(fdata[2])/scal[2]	#CG4
                    vals[n,4] = float(fdata[5])/scal[3]	#CHP1
                    #Temperatures and corrections
                    R_CG4 = float(fdata[3])/1000.	#kohms. 10^5 ohms is sensible highest
                    if R_CG4 <= 100.0:
                        T_CG4 = (A*R_CG4**5 + B*R_CG4**4 + C*R_CG4**3 + D*R_CG4**2 + E*R_CG4 + F) + Tk
                        corr_CG4 = SBconst * T_CG4 ** 4               
                        vals[n,2] = vals[n,2] + corr_CG4
                        vals[n,3] = T_CG4
                    else:
                        vals[n,2] = missing_value	#This it should be missing as can't give it a value? Set a value for QC flag too?
                        vals[n,3] = missing_value
                    R_CHP1 = float(fdata[6])
                    if R_CHP1 <= 100000.0:
                        lnR_CHP1 = np.log(R_CHP1);
                        vals[n,5] = 1.0/(A1 + (B1 * lnR_CHP1) + (G1 * lnR_CHP1 ** 3)) 
                    else:
                        vals[n,5] = missing_value

                    n += 1

            f.close()
    
        print('no. values  = ', n)

    if n > 0:
        # -------------------------------------
        # Trim arrays to remove unused elements
        # -------------------------------------
        vals = vals[0:n,:]
        print('vals.shape = ',vals.shape)
        timesecs   = timesecs[0:n]
        input_missing = input_missing[0:n,:]

        #Assign values of input_missing depending on location of missing_value in vals array
        input_missing = np.where(vals == missing_value,1,0)

        print(timesecs[0:9])
        print(vals[0:9,0])
        print(vals[0:9,1])
        print(vals[0:9,2])
        print(vals[0:9,3])
        print(vals[0:9,4])
        print(vals[0:9,5])

        sampling_interval = str(timesecs[1]-timesecs[0])+' seconds'
        t_interval = timesecs[1]-timesecs[0]

        #--------------------------------------------------------
        #Read corrections file
        #Returns qc_flag compatible with amof netCDF
        #--------------------------------------------------------

        #qc_flag = load_netcdf_corrections(nday, chids, n, t_interval, vals)
        corr_details = load_netcdf_corrections(nday, chids, n, t_interval, timesecs, vals)
        qc_flag = corr_details[0]
        for k in range(nvar):
            sub_vals = vals[:,k]	#Data for each instrument in file
            for m in range(n):
                if input_missing[m,k] == 1:
                    qc_flag[m,k] = 2
            vals_ok = sub_vals[np.where(qc_flag[:,k] == 1)]
            if len(vals_ok) > 0:
                valid_min_max[0,k] = np.amin(vals_ok)
                valid_min_max[1,k] = np.amax(vals_ok)

        #valid_min_max = corr_details[1]
        print('qc_flag shape = ',qc_flag.shape)

        #vals_qc = np.where(qc_flag != 2, vals, 0)
        #vals_qc = np.where(qc_flag != 3, vals_qc, missing_value)	#Will over-ride HOLDCAL=0 set in previous line
        vals = np.where(qc_flag != 2, vals, missing_value)


        print(qc_flag[0:49,0])

        #Need to decide how we treat missing values. If missing from all datasets do we not report those times?
        #Affects how we quality flag those values 


    #--------------------------------------------
    # Sorting out day/time information 
    #--------------------------------------------

        time_details = generate_netcdf_datetimeinfo(year,month,day,n,timesecs)

        #return yr_arr,mn_arr,day_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs
        yr_arr = time_details[0]
        mn_arr = time_details[1]
        dy_arr = time_details[2]
        dyyr_arr = time_details[3]
        hr_arr = time_details[4]
        mi_arr = time_details[5]
        sc_arr = time_details[6]
        epoch_timesecs = time_details[7]
        first_last_datetime = time_details[8]
        print(yr_arr[0:9])
        print(mn_arr[0:9])
        print(dy_arr[0:9])
        print(epoch_timesecs[0:9])
        print(hr_arr[0:9])
        print(mi_arr[0:9])
        print(sc_arr[0:9])
        print(dyyr_arr[0:9])

        uname              = os.uname()
        nodename           = uname[1]
        os_name            = uname[0]
        os_release         = uname[2]
        computer_id        = nodename + ' under ' + os_name + ' ' + os_release
        user_id            = pwd.getpwuid(os.getuid())[0]

        substitutions = {"last_revised_date": datetime.datetime.now(), "start_time": datetime.datetime(first_last_datetime[0,0],first_last_datetime[1,0],first_last_datetime[2,0],first_last_datetime[3,0],first_last_datetime[4,0],first_last_datetime[5,0]), "end_time": datetime.datetime(first_last_datetime[0,1],first_last_datetime[1,1],first_last_datetime[2,1],first_last_datetime[3,1],first_last_datetime[4,1],first_last_datetime[5,1]), "conversion_time": datetime.datetime.utcnow(), "observation_date": datetime.date(year,month,day), "start_unix_time": np.amin(epoch_timesecs), "end_unix_time": np.amax(epoch_timesecs), "start_dyyr": np.amin(dyyr_arr),"end_dyyr": np.amax(dyyr_arr), "start_yr": np.amin(yr_arr), "end_yr": np.amax(yr_arr), "start_mon": np.amin(mn_arr), "end_mon": np.amax(mn_arr), "start_day": np.amin(dy_arr), "end_day": np.amax(dy_arr), "start_hour": np.amin(hr_arr), "end_hour": np.amax(hr_arr), "start_min": np.amin(mi_arr), "end_min": np.amax(mi_arr), "start_sec": np.amin(sc_arr), "end_sec":np.amax(sc_arr), "CM21_min": valid_min_max[0,0], "CM21_max": valid_min_max[1,0], "CMP21_min": valid_min_max[0,1], "CMP21_max": valid_min_max[1,1], "CG4_min": valid_min_max[0,2], "CG4_max": valid_min_max[1,2], "CG4_T_min": valid_min_max[0,3], "CG4_T_max": valid_min_max[1,3], "CHP1_min": valid_min_max[0,4], "CHP1_max": valid_min_max[1,4], "CHP1_T_min": valid_min_max[0,5], "CHP1_T_max": valid_min_max[1,5], "uname": user_id, "codename": __file__, "machinename": computer_id}


        for nn in range(nvar-2):	#Each data file

            #----------------------------
            # Setting up YAML data object
            #----------------------------
            #Need export PYTHONPATH=/home/jla/python/global_modules
            handler = module_data_object_python3.Handler()
            data_object_id = handler.load_template(tfp[nn])
            print(data_object_id)
            #handler.show_requirements_for_template(data_object_id)


            # -------------------------------
            # Open new netCDF file for output
            # -------------------------------
            #out_file   = os.path.join(path_out[nn],datestring+data_version+'.nc')
            out_file   = path_out[nn]+datestring+data_product+data_version+'.nc'
            print('Opening new NetCDF file ' + out_file)

            lengths_of_dimensions = {"time": n}

            data_object = handler.return_data_object_from_template(data_object_id,lengths_of_dimensions,substitutions)
            data_object["variables"]["time"]["values"][:] = epoch_timesecs[:]
            data_object["variables"]["day_of_year"]["values"][:] = dyyr_arr[:]	#Need to include fraction of day
            data_object["variables"]["year"]["values"][:] = yr_arr[:]
            data_object["variables"]["month"]["values"][:] = mn_arr[:]
            data_object["variables"]["day"]["values"][:] = dy_arr[:]
            data_object["variables"]["hour"]["values"][:] = hr_arr[:]
            data_object["variables"]["minute"]["values"][:] = mi_arr[:]
            data_object["variables"]["second"]["values"][:] = sc_arr[:]


            if nn == 0:
                chids = ["pyrCM21_ch"]
                #Call corrections, but may want to re-order so only add that code once
                data_object["variables"]["downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,0]
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,0]
                variables_to_plot = ["downwelling_shortwave_flux_in_air"]
            if nn == 1:
                chids = ["pyr_CMP21_ch"]
                data_object["variables"]["diffuse_downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,1]	#Should diffuse be added to name?
                data_object["variables"]["qc_flag"]["values"][:] = qc_flag[:,1]
                variables_to_plot = ["diffuse_downwelling_shortwave_flux_in_air"]
            if nn == 2:
                chids = ["pyrCG4_ch, pyrCG4_ch"]
                data_object["variables"]["downwelling_longwave_flux_in_air"]["values"][:] = vals[:,2] 
                data_object["variables"]["body_temperature"]["values"][:] = vals[:,3] 
                data_object["variables"]["qc_flag_downwelling_longwave_flux_in_air"]["values"][:] = qc_flag[:,2]
                data_object["variables"]["qc_flag_body_temperature"]["values"][:] = qc_flag[:,3]
                variables_to_plot = ["downwelling_longwave_flux_in_air", "body_temperature"]
            if nn == 3:
                chids = ["pyrCP1_ch, pyrCP1_T_ch"]
                data_object["variables"]["direct_downwelling_shortwave_flux_in_air"]["values"][:] = vals[:,4] 
                data_object["variables"]["body_temperature"]["values"][:] = vals[:,5] 
                data_object["variables"]["qc_flag_direct_downwelling_shortwave_flux_in_air"]["values"][:] = qc_flag[:,4]
                data_object["variables"]["qc_flag_body_temperature"]["values"][:] = qc_flag[:,5]
                variables_to_plot = ["direct_downwelling_shortwave_flux_in_air", "body_temperature"]
            exit_code = handler.write_data_object_to_netcdf_file(out_file, data_object)
            print(exit_code)

            generate_netcdf_graphs(out_file, graph_path_out[nn], variables_to_plot, datestring)



# ------------------------------------------------------------------------
# Define the general netcdf generation function, common to all sensors
# ------------------------------------------------------------------------
#def generate_netcdf_common(nday,year,month,day,datestring):
def generate_netcdf_common(nday):

    # -----------------------------
    # Date-time numbers and strings
    # -----------------------------
    ndat             = num2date(nday)
    year             = ndat.year
    month            = ndat.month
    day              = ndat.day
    dat              = datetime.datetime(year,month,day,0,0,0)
    nowdat           = datetime.datetime.now()
    nowstring        = nowdat.strftime("%Y-%m-%d %H:%M:%S")
    datestring       = dat.strftime("%Y%m%d")
    start_of_day_str = dat.strftime("%Y-%m-%d %H:%M")

    print(year,month,day,datestring)
    return year,month,day,datestring,start_of_day_str,nowstring



# ------------------------------------------------------------------------
# Define the general netcdf generation function, common to all sensors
# ------------------------------------------------------------------------
def generate_netcdf_inputfiles(path_in,datestring):


    if os.access(path_in, os.F_OK):	#A data directory exists for the day
        files = os.listdir(path_in)
        tempstr = "pldc_"+datestring+'.00'
        expr = re.compile(tempstr+'[0-9]')
        infiles = [elem for elem in files if expr.search(elem)]
        infiles.sort()
        print(infiles)
        nfiles = len(infiles)
        if nfiles == 0:
            print('No files on this day')

    else:
        nfiles = 0
        print('No directory found')

    return nfiles,infiles



#--------------------------------------------
# Sorting out day/time information 
#--------------------------------------------
def generate_netcdf_datetimeinfo(year,month,day,n,timesecs):

    first_last_datetime = np.zeros((6,2),dtype = np.int32)      #Values to use in datetime.datetime, derived from time.gmtime
    tt = (year,month,day,0,0,0,0,0,0)
    tt_upper = (year,month,(day+1),0,0,0,0,0,0)
    day_start_epoch = time.mktime(tt)
    yr_arr=np.zeros(n)+year
    mn_arr=np.zeros(n)+month
    dy_arr=np.zeros(n)+day
    dyyr_arr=np.zeros(n)
    epoch_timesecs=np.zeros(n)
    hr_arr=np.zeros(n,dtype=np.int32)
    mi_arr=np.zeros(n,dtype=np.int32)
    sc_arr=np.zeros(n,dtype=np.float32)
    epoch_timesecs = day_start_epoch + timesecs
    #In cases where parts of previous or next day are found in file, it would be best to derive the above from the data timestamp
    for mm in range(n):
        tt = time.gmtime(epoch_timesecs[mm])
        hr_arr[mm] = np.int32(tt[3])
        mi_arr[mm] = np.int32(tt[4])
        sc_arr[mm] = timesecs[mm]-3600*hr_arr[mm]-60*mi_arr[mm]	#Keeps digits after decimal place, unlike tt[5] 
        dyyr_arr[mm] = np.float32(tt[7]+timesecs[mm]/86400.0)
        #If the last time is exactly midnight, we need different in
        if timesecs[mm] == 86400.0:     #Point from midnight of the following day
            sc_arr[mm] = sc_arr[mm] - 86400.0
            hr_arr[mm] = 24
            dyyr_arr[mm] = np.float32(tt[7])    #Will be number for next day
        #These values are used for first, last times in the global attributes
        #It's safer to derive them from the tt array, as using the arrays can cause problems with the datetime.datetime function
        #if values aren't compatible, e.g. last hr_arr value being 24
        if mm == 0:     #First point
            for nn in range(6):
                first_last_datetime[nn,0] = np.int32(tt[nn])
        if mm == n-1:   #Last point
            for nn in range(6):
                first_last_datetime[nn,1] = np.int32(tt[nn])

    print(np.amin(hr_arr),np.amax(hr_arr),np.amin(mi_arr),np.amax(mi_arr),np.amin(sc_arr),np.amax(sc_arr))
    #print(sc_arr)
    print(first_last_datetime)

    return yr_arr,mn_arr,dy_arr,dyyr_arr,hr_arr,mi_arr,sc_arr,epoch_timesecs,first_last_datetime#,lengths_of_dimensions,substitutions



#--------------------------------------------
# Reading corrections file 
#--------------------------------------------
def load_netcdf_corrections(nday, chids, n_values, t_interval, timesecs, vals):

    n_chids = len(chids)
    qualflag = np.ones((n_values,n_chids), dtype = int)
    valid_min_max = np.zeros((2,n_chids))

    for n_inst in range(n_chids):
        corr_file = '/data/netCDF/corrections/' + chids[n_inst] + '.corr'
        print('n_inst, corr_file = ', n_inst, corr_file)

        sub_vals = vals[:, n_inst]	#Data for each instrument in file
        f = open(corr_file, 'r')

        #Check for HMP155 purgetimes first, as want any points flagged as purge time to be superseded by baddata from a correction file if it exists
        #Get hmp155 purge times for oatnew_ch and rhnew_ch
        if np.char.find(chids[n_inst], "oatnew") >= 0 or np.char.find(chids[n_inst], "rhnew") >= 0:
            print('HMP155 purge times')
            purge_limit = load_hmp155_purgetime(nday)
            start_purge = int(round(purge_limit[0]/t_interval))
            end_purge = int(round(purge_limit[1]/t_interval))
            print('start_purge, end_purge = ', start_purge, end_purge)
            qualflag[start_purge:(end_purge+1),n_inst] = 3



        while True:	#Reading correction file

            line = f.readline().strip()
            if not line: break
            if len(line) == 30:	#Standard length for correction file
            #if len(line) == 31:	#Standard length for correction file
                startnum = int(date2num(datetime.datetime(int(line[0:4]),int(line[4:6]),int(line[6:8]),0,0,0)))
                #print 'nday, startnum = ',nday, startnum
                if startnum == nday:
                    print(len(line), line)
                    print(len(line.strip()), line.strip())
                    t_start = 3600.0*float(line[9:11]) + 60.0*float(line[11:13]) + float(line[13:15])
                    t_end   = 3600.0*float(line[16:18]) + 60.0*float(line[18:20]) + float(line[20:22])
               	    print(t_start, t_end, t_interval)
               	    #start_index = int(round(t_start/t_interval))	#Time range in file is 0 to 86390 
                    #end_index = int(round(t_end/t_interval))	#Don't think we need +1 here to round up.
               	    start_index = np.argmin(np.absolute(timesecs-t_start))	#Time range in file is 0 to 86390 
                    end_index = np.argmin(np.absolute(timesecs-t_end))	#Time range in file is 0 to 86390 
                    #end_index = np.where(t_end == timesecs)		#Don't think we need +1 here to round up.
                    print(start_index, end_index, n_values)
                    if end_index >= (n_values-2):       #May be rounding errors with max. index so reduce if necessary. Also avoid last point in day being missed.
                        end_index = n_values - 1
                    if start_index == 1:        #Avoid having 1 single point remaining at start of day
                        start_index = 0
                    print(start_index, end_index, n_values)
                    instring = line[23:]
                    if np.char.equal(instring,'HOLDCAL'):
                        if np.char.find(chids[n_inst],"rg") >= 0:
                            print("qualflag for a raingauge")	#For a raingauge, only flag values as HOLDCAL if they're > 0
                            for flag_ind in range((end_index-start_index+1)):
                                #Set qualflag to 3 if the data value is greater than 0 and the flag isn't already set to 2 (BADDATA) from a previous line
                                if sub_vals[flag_ind+start_index] > 0 and qualflag[(flag_ind+start_index),n_inst] <= 1:
                                    qualflag[(flag_ind+start_index),n_inst] = 3
                                    print(flag_ind+start_index, sub_vals[flag_ind+start_index])

                    if np.char.equal(instring,'BADDATA'):
                        qualflag[start_index:(end_index+1),n_inst] = 2	#Need end_index+1 to make change all points from start_index to (and including) end_index
                    print('t_start, t_end, start_index, end_index, instring = ', t_start, t_end, start_index, end_index, instring)

        f.close()

        vals_ok = sub_vals[np.where(qualflag[:,n_inst] == 1)]
        if len(vals_ok) > 0:
            valid_min_max[0,n_inst] = np.amin(vals_ok)
            valid_min_max[1,n_inst] = np.amax(vals_ok)

    return qualflag, valid_min_max


#--------------------------------------------
# Reading hmp155 purgetimes file 
#--------------------------------------------
def load_hmp155_purgetime(nday):

    purge_limit = np.zeros((2))

    hmp_file = '/data/netCDF/corrections/hmp155_purgetime.txt'

    f = open(hmp_file, 'r')

    while True:

        line = f.readline()
        if not line: break
        fdata = line.split(' ')
        date_start = int(date2num(datetime.datetime(int(line[0:4]),int(line[4:6]),int(line[6:8]),0,0,0)))
        date_end = int(date2num(datetime.datetime(int(line[16:20]),int(line[20:22]),int(line[22:24]),0,0,0)))
        if nday > date_start and nday <= date_end:
            print(line)
            purge_limit[0] = 3600.0*float(fdata[4])
            purge_limit[1] = 3600.0*float(fdata[5])
            print('purge_limit = ', purge_limit)

    f.close()

    return purge_limit 


#-----------------------------------------------
# Reading broadband radiometer calibration files 
#-----------------------------------------------
def load_bbrad_calibrations(nday):

    scal = np.zeros((4))

    cm21_file = '/data/netCDF/corrections/CM21_calibration.txt'
    cmp21_file = '/data/netCDF/corrections/CMP21_calibration.txt'
    cg4_file = '/data/netCDF/corrections/CG4_calibration.txt'
    chp1_file = '/data/netCDF/corrections/CHP1_calibration.txt'
    file_in = [cm21_file, cmp21_file, cg4_file, chp1_file]

    for n in range(4):

        f = open(file_in[n], 'r')

        while True:

            line = f.readline()
            if not line: break
            fdata = line.split(' ')
            date_start = int(date2num(datetime.datetime(int(line[0:4]),int(line[4:6]),int(line[6:8]),0,0,0)))
            date_end = int(date2num(datetime.datetime(int(line[9:13]),int(line[13:15]),int(line[15:17]),0,0,0)))
            print(date_start, date_end)
            if nday >= date_start and nday < date_end:
                print(line)
                scal[n] = float(fdata[2])

        f.close()
    print('scal = ',scal)
    return scal 





#--------------------------------------------
# Generate a plot of the data 
#--------------------------------------------
def generate_netcdf_graphs(out_file, graph_path_out, variables_to_plot, datestring):

    n_plots = len(variables_to_plot)
    
    tt = (int(datestring[0:4]),int(datestring[4:6]),int(datestring[6:8]),0,0,0,0,0,0)
    day_start = time.mktime(tt)

    ncfile = nc4.Dataset(out_file, 'r',format='NETCDF4_CLASSIC')
    splitstr = out_file.split('/')	#Split data filename at /, want last string of this array (filename)
    out_name = splitstr[len(splitstr)-1]
    point_pos = out_name.rfind('.')
    uscore_pos = out_name.find('_')
    print(out_name, uscore_pos)
    #plot_file_root = out_name[0:point_pos]
    plot_file_root = out_name[0:uscore_pos+1]
    temp_time = ncfile.variables['time']
    var_time = temp_time[:]
    #var_time = (var_time - var_time[0])/3600.0
    var_time = (var_time - day_start)/3600.0
    n_points = var_time.size
    print(n_points,var_time[0:9])
    x_axis_title = 'Time (hours since midnight)' 

    for n_meas in range(n_plots):

        plot_2d = 0	#Default is a 1d plot, but e.g. disdrometer counts data require a 2d plot
        #Derive output plot filename and create directory for it if it doesn't exist
        #out_plot_file = '/home/jla/pythonplot' + str(n_meas) + '.png'
        out_plot_path = graph_path_out + variables_to_plot[n_meas] + '/'
        #If this directory doesn't exist, create it
        if not os.path.isdir(out_plot_path):
            print("Creating plot directory ", out_plot_path)
            os.makedirs(out_plot_path)
            oscommand = "chgrp -R netcdf " + out_plot_path
            os.system(oscommand) 
            oscommand = "chmod -R g+w " + out_plot_path
            os.system(oscommand) 
        #out_plot_file = out_plot_path + plot_file_root + '_' + variables_to_plot[n_meas] + '.png'
        out_plot_file = out_plot_path + datestring + '_' + plot_file_root + variables_to_plot[n_meas] + '.png'  #Shorter name, starting with date, better for quicklooks and thumbnails
        print('out_plot_file = ',out_plot_file)

        #Mark any flagged points with a separate marker. Needs qc_flag name and values.
        qc_flag_read = 'qc_flag'	#Most are called qc_flag, but not all
        if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' and out_file.find('rain-gauge-4') >= 0:
            qc_flag_read = 'qc_flag_thickness_of_rainfall_amount'
        if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' and out_file.find('disdrometer') >= 0:
            qc_flag_read = 'qc_flag_thickness_of_rainfall_amount'
        if variables_to_plot[n_meas] == 'number_of_hydrometeors_per_size_channel' and out_file.find('disdrometer') >= 0:
            plot_2d = 1	#Requires a 2d plot 
            qc_flag_read = 'qc_flag_hydrometeors'
        if variables_to_plot[n_meas] == 'rainfall_rate':
            qc_flag_read = 'qc_flag_rainfall_rate'
        if variables_to_plot[n_meas] == 'air_temperature':
            qc_flag_read = 'qc_flag_air_temperature'
        if variables_to_plot[n_meas] == 'relative_humidity':
            qc_flag_read = 'qc_flag_relative_humidity'
        if variables_to_plot[n_meas] == 'wind_speed':
            qc_flag_read = 'qc_flag_wind_speed'
        if variables_to_plot[n_meas] == 'wind_from_direction':
            qc_flag_read = 'qc_flag_wind_from_direction'
        if variables_to_plot[n_meas] == 'downwelling_longwave_flux_in_air':
            qc_flag_read = 'qc_flag_downwelling_longwave_flux_in_air'
        if variables_to_plot[n_meas] == 'direct_downwelling_shortwave_flux_in_air':
            qc_flag_read = 'qc_flag_direct_downwelling_shortwave_flux_in_air'
        if variables_to_plot[n_meas] == 'body_temperature':
            qc_flag_read = 'qc_flag_body_temperature'



        #Read variable data from file
        temp_var_y = ncfile.variables[variables_to_plot[n_meas]]
        var_y = temp_var_y[:]
        y_axis_title = ncfile.variables[variables_to_plot[n_meas]].long_name + ' (' + ncfile.variables[variables_to_plot[n_meas]].units + ')'
        if plot_2d == 1:
            temp_sizes = ncfile.variables["diameter"]
            y_ds_scale = temp_sizes[:]
            n_sizes = y_ds_scale.size
            y_axis_title = ncfile.variables["diameter"].long_name + ' (' + ncfile.variables["diameter"].units + ')'
            cbar_label = ncfile.variables[variables_to_plot[n_meas]].long_name + ' (' + ncfile.variables[variables_to_plot[n_meas]].units + ')'
        temp_qc_flag = ncfile.variables[qc_flag_read]
        qc_flag_in = temp_qc_flag[:]
        print(x_axis_title, y_axis_title)
        fig_title = ncfile.source + ' ' + datestring
        
        if plot_2d == 0:
            flagged_time = var_time[np.where(qc_flag_in > 1)]
            flagged_var_y = var_y[np.where(qc_flag_in > 1)]
            missing_value = ncfile.variables[variables_to_plot[n_meas]]._FillValue
            #y_axis_title = ncfile.variables[variables_to_plot[n_meas]].long_name + ' (' + ncfile.variables[variables_to_plot[n_meas]].units + ')'
            #fig_title = ncfile.source + ' ' + datestring
            #print(var_y[0:9])
            #print(x_axis_title, y_axis_title)
            #Calculate % missing values
            suspect_vals = 100.0 * flagged_var_y.size/n_points
            suspect_text = "Suspect: %.3f %%" % (suspect_vals) 
            print('suspect_text = ', suspect_text)

        else:	#2d plot
            flagged_var_y = var_y[np.where(qc_flag_in > 1)]
            missing_value = ncfile.variables[variables_to_plot[n_meas]]._FillValue
            suspect_vals = 100.0 * flagged_var_y.size/(n_points * n_sizes)
            suspect_text = "Suspect: %.3f %%" % (suspect_vals) 
            print('suspect_text = ', suspect_text)
            

        plt.ion()
        plt.ioff()
        fig1 = plt.subplots(nrows=1)
        if plot_2d == 0:
        #plt.figure()
            plt.plot(var_time, var_y, linewidth=1, label = suspect_text)	#see https://matplotlib.org/2.0.0/api/pyplot_api.html#matplotlib.pyplot.legend and example code for more than 1 plot on figure
            plt.plot(flagged_time, flagged_var_y, 'ro', fillstyle = 'none', markersize = 4)	#ro- joins points together, not wanted
            plt.xlabel(x_axis_title)
            plt.ylabel(y_axis_title)
            plt.xlim(0, 24)
            if variables_to_plot[n_meas] == 'wind_from_direction':
                plt.ylim(0, 360)
            if variables_to_plot[n_meas] == 'thickness_of_rainfall_amount' or variables_to_plot[n_meas] == 'rainfall_rate':
                print('Rainfall plot')
                max_y = np.amax(var_y)
                print('max_y = ', max_y)
                if max_y <= 0.001:
                    print('Using minimum y-axis = 0')
                    min_y = -0.0001	#So that can see line at zero
                    max_y = 0.001    #Arbitrary low value
                    plt.ylim(min_y, max_y)
 
            plt.xticks( 4.0 * arange(7) )
            plt.suptitle(fig_title)
            plt.minorticks_on()
            plt.grid()
            plt.legend(loc=9)

        else:	#plot_2d == 1
            #levels = MaxNLocator(nbins=15).tick_values(0.0, np.ceil(var_y.max()))    #Otherwise min is -Inf when plotting logs
            print('np.amin, np.amax counts = ', np.amin(var_y), np.amax(var_y))
            print('var_y.shape = ',var_y.shape)
            #Set points where qc_flag != 1 to missing value
            #Can't use var_y = np.where(qc_flag_in == 1 , var_y, missing_value_array) as the output is an array, not values
            var_y[np.where(qc_flag_in > 1)] = missing_value	
            #for m in range(n_points):	#Using long-winded, clunky way instead!
                #for k in range(n_sizes):
                   #if qc_flag_in[m, k] != 1:
                       #var_y[m, k] = missing_value
            print('After qc_flag applied, np.amin, np.amax counts = ', np.amin(var_y), np.amax(var_y))

            #print('Ch disdro levels = ',levels)
            cmap = plt.get_cmap('jet')
            cmap.set_under("white")
            #Set up range for colorbar
            #For disdrometer, I'm using steps of 10 counts for vmax, but might want to make this more general for other insts.
            min_cbar = 0	#Can set this to e.g. 0.01, then zero values are white, but then they can't be distinguished from missing values
            max_cbar = 6
            #print('max(var_y), max_cbar = ',np.amax(var_y), max_cbar)
            #if max_cbar <= 1.0:
                #max_cbar = 5.0
            #print('max(var_y), max_cbar = ',np.amax(var_y), max_cbar)
            #A custom colour map option is available. See /home/jla/matlab/hoganjet.dat for details of ours. Has anyone coded it for Python

            fig, ax = plt.subplots(figsize=[12,7])

            #plt.pcolormesh(var_time, y_ds_scale, np.transpose(var_y), cmap=cmap, vmin = min_cbar, vmax = max_cbar)       #Remove norm, so that vmin,vmax are used
            cax = ax.pcolormesh(var_time, y_ds_scale, np.transpose(var_y), cmap=cmap, vmin = min_cbar, vmax = max_cbar)       #Remove norm, so that vmin,vmax are used
            cbar = fig.colorbar(cax, ticks = [0, 1, 2, 3, 4, 5, 6])
            #cbar = colorbar(ticks = [0,1,2,3,4,5])
            cbar.ax.set_yticklabels(['$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$', '$10^{4}$', '$10^{5}$', '$10^{6}$'])
            cbar.set_label(cbar_label)
            ax.set_xlabel(x_axis_title)
            ax.set_ylabel(y_axis_title)
            ax.set_xlim(0, 24)
            ax.set_xticks( 4.0 * arange(7) )
            ax.set_title(fig_title)
            ax.minorticks_on()
            #ax.xaxis.set_minor_locator(MultipleLocator(4))
            

        plt.savefig(out_plot_file,format='png')


    ncfile.close()

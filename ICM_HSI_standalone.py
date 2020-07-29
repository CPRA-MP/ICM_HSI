def daily2ave(all_sd,ave_sd,ave_ed,input_file,input_nrows=-9999): 
    # this function reads a portion of the ICM-Hydro daily timeseries file into a numpy array and then computes the average for the time slice read in
    # this function returns a dictionary 'comp_ave' that has ICM-Hydro compartment ID as the key and a temporal average for each compartment as the value
    # the key is of type integer and the values are of type float
    
    # if looping through the whole file and batch generating averages for a bunch of timeslices it will be faster to read in the whole file to a numpy array and iterating over the whole array rather than iteratively calling this function
    
    # 'all_sd' is the start date of all data included in the daily timeseries file - all_sd is a datetime.date object   
    # 'ave_sd' is the start date of averaging window, inclusive - ave_sd is a datetime.date object                       
    # 'ave_ed' is the end date of averaging window, inclusive - ave_ed is a datetime.date object                           

    # check if number of rows in input file was passed into function
    if input_nrows == -9999:                # if not passed in, calculate number of rows by calling file_len function
        all_rows = file_len(input_file)
    except:                                 # if passed in, do not call file_len function
        all_rows = input_nrows

    ave_n = (ave_ed - ave_sd).days + 1  # number of days to be used for averaging
    skip_head = (ave_sd - all_sd).days  # number of rows at top of daily timeseries to skip until start date for averaging window is met
    skip_foot = all_rows - skip_head - ave_n
    data = np.genfromtxt(input_file,dtype='str',delimiter=',',skip_header=skip_head,skip_footer=skip_foot)
    comp_ave = {}
    nrow = 0
    for row in data:
        if nrow == 0:
            for comp in range(1,len(row)+1):
                comp_ave[comp]=0.0
        for col in range(0,len(row)):
            comp = col + 1
            val = float(row[col])
            comp_ave[comp] += val/ave_n
        nrow += 1
    return comp_ave


def file_len(fname):
    # this function counts the number of lines in a text file
    # this function returns an integer value that is the number of lines in the file 'fname'
    
    # 'fname' is a string variable that contains the full path to a text file with an unknown number of lines
    
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def comp2grid(comp_data_dict,grid_comp_dict):
    # this function maps ICM-Hydro compartment level data to the 500-m grid
    # this function returns a dictionary 'grid_data' that has grid ID as the key and the respective compartment-level data as the value
    
    # 'comp_data_dict' is a dictionary with ICM-Hydro compartment as the key and some value to be mapped to the grid as the value
    # 'grid_comp_dict' is a dictionary with grid ID as the key and the corresponding ICM-Hydro compartment number as the value
    
        
    grid_data = {}
    for gid in grid_comp_dict.keys():
        cid = grid_comp_dict[gid]
        grid_data[gid] = comp_data_dict[cid]
    
    return grid_data


def compout2dict(input_file,import_column): 
    # this function reads the compartment-based summary output file 'input_file' into a dictionary
    # the first row in input_file must contain header text - this will be skipped on import
    # the first column in 'input_file' must contain ICM-Hydro compartment ID numbers (this will be used as keys in the dict)
    # the key is of type integer and the values are of type float
    
    # 'import_column'  is the zero-indexed column number that contains the data to be imported and mapped
    # import_column = 7 will import the the maximum 2-week mean salinity in compartment_out_YYYY.csv
    
    # this function returns a dictionary 'comp_ave' that has ICM-Hydro compartment ID as the key and a single average for each compartment as the value

    data = np.genfromtxt(input_file,dtype='str',delimiter=',',skip_header=1)
    comp_ave = {}
    nrow = 0
    for row in data:
        comp = int(float(row[0]))
        val = float(row[import_column])
        comp_ave[comp] = val
        nrow += 1

##################################
####
####
####
####    Run ICM standalone for HSI
####
####
####
##################################

import os
import sys
import exceptions
import shutil
import math
import time
import errno
import numpy as np
import random

inputs = np.genfromtxt('ICM_control.csv',dtype=str,comments='#',delimiter=',')

# Parent directory locations for various ICM components
# These directories must exist in the model folder
# Other directories are created throughout the ICM run - but they are all based on these parent directories
ecohydro_dir = os.path.normpath(inputs[1,1].lstrip().rstrip())
wetland_morph_dir = os.path.normpath(inputs[2,1].lstrip().rstrip())
vegetation_dir = os.path.normpath(inputs[3,1].lstrip().rstrip())
bimode_dir = os.path.normpath(inputs[4,1].lstrip().rstrip())
HSI_dir = os.path.normpath(inputs[5,1].lstrip().rstrip())
ewe_dir = os.path.normpath(inputs[6,1].lstrip().rstrip())
# Configuration files used by various ICM components
VegConfigFile = inputs[7,1].lstrip().rstrip()
WMConfigFile = inputs[8,1].lstrip().rstrip()
EHConfigFile = inputs[9,1].lstrip().rstrip()
EHCellsFile = inputs[10,1].lstrip().rstrip()
EHLinksFile = inputs[11,1].lstrip().rstrip()
BIPrismFile = inputs[12,1].lstrip().rstrip()
BIMHWFile = inputs[13,1].lstrip().rstrip()
EHInterfaceFile = inputs[14,1].lstrip().rstrip()
BMInterfaceFile = inputs[15,1].lstrip().rstrip()
compartment_output_file = 'compartment_out.csv'
grid_output_file = 'grid_500m_out.csv'

# Filenames for Veg model input
WaveAmplitudeFile = inputs[16,1].lstrip().rstrip()
MeanSalinityFile = inputs[17,1].lstrip().rstrip()
SummerMeanWaterDepthFile = inputs[18,1].lstrip().rstrip()
SummerMeanSalinityFile = inputs[19,1].lstrip().rstrip()
SummerMeanTempFile = inputs[20,1].lstrip().rstrip()
TreeEstCondFile = inputs[21,1].lstrip().rstrip()
HtAbvWaterFile = inputs[22,1].lstrip().rstrip()
PerLandFile = inputs[23,1].lstrip().rstrip()

## Simulation Settings
startyear = int(inputs[24,1].lstrip().rstrip())
endyear = int(inputs[25,1].lstrip().rstrip())
nvegtype = int(inputs[26,1].lstrip().rstrip())
inputStartYear = int(inputs[27,1].lstrip().rstrip())
hotstart_year = int(inputs[28,1].lstrip().rstrip())
elapsed_hotstart = hotstart_year - startyear
update_hydro_attr = int(inputs[29,1].lstrip().rstrip())

## grid information for Veg ASCII grid files
n500grid= int(inputs[30,1].lstrip().rstrip())
# n500gridveg = int(inputs[25,1].lstrip().rstrip()) #total number of grid cells in Veg model - including NoData cells
n500rows = int(inputs[31,1].lstrip().rstrip())
n500cols = int(inputs[32,1].lstrip().rstrip())
xll500 = int(inputs[33,1].lstrip().rstrip())
yll500 = int(inputs[34,1].lstrip().rstrip())

## grid information for EwE ASCII grid files
n1000grid = int(inputs[35,1].lstrip().rstrip())
n1000rows = int(inputs[36,1].lstrip().rstrip())
n1000cols = int(inputs[37,1].lstrip().rstrip())
xll1000 = inputs[38,1].lstrip().rstrip()
yll1000 = inputs[39,1].lstrip().rstrip()

# file naming settings
mpterm = inputs[40,1].lstrip().rstrip()
sterm = inputs[41,1].lstrip().rstrip()
gterm = inputs[42,1].lstrip().rstrip()
cterm = inputs[43,1].lstrip().rstrip()
uterm = inputs[44,1].lstrip().rstrip()
vterm = inputs[45,1].lstrip().rstrip()
rterm = inputs[46,1].lstrip().rstrip()
runprefix = '%s_%s_%s_%s_%s_%s_%s' % (mpterm,sterm,gterm,cterm,uterm,vterm,rterm)

EHtemp_path = os.path.normpath(r'%s\\TempFiles' % ecohydro_dir)

## read Wetland Morph parameters csv file into array (first column is descriptor, second column is variable)                            
WM_params = np.genfromtxt(WMConfigFile,dtype=str,delimiter=',',usecols=1)   

# read in grid-to-compartment lookup table into a dictionary
# key is grid ID and value is compartment
grid_lookup_file = r'%s\hydro\grid_lookup_500m.csv' % input_directory
grid_lookup = np.genfromtxt(grid_lookup_file,skip_header=1,delimiter=',',dtype='int',usecols=[0,1])
grid_comp = {row[0]:row[1] for row in grid_lookup}
grid500_res = 500.0

# Save list of GridIDs into an array for use in some loops later # ultimately can replace with grid_comp.keys() in loops
gridIDs=grid_comp.keys():

for year in range(startyear,endyear+1):
    
    elapsedyear = year - startyear + 1
    
    veg_output_file = '%s_O_%02d_%02d_V_vegty.asc+' % (runprefix,elapsedyear,elapsedyear)
    veg_output_filepath = os.path.normpath(vegetation_dir + '\\' + veg_output_file)
        
    dom = {}
    dom[1]=31
    if year in range(2000,4000,4):
        dom[2]=29
    else:
        dom[2]=28
    dom[3]=31
    dom[4]=30
    dom[5]=31
    dom[6]=30
    dom[7]=31
    dom[8]=31
    dom[9]=30
    dom[10]=31
    dom[11]=30
    dom[12]=31
    
    # update name of grid data file generated by Morph to include current year in name
    new_grid_file = 'grid_data_500m_end%s.csv' % (year)  # this must match name set in "WM.CalculateEcohydroAttributes" with the exception of (year) here instead of CurrentYear
    new_grid_filepath = os.path.normpath('%s\\%s' % (EHtemp_path,new_grid_file)) # location of grid file after it is generated in "WM.CalculateEcohydroAttributes"

    ##############################################
    ##    HABITAT SUITABILITY INDICES ~ HSIs    ##
    ##############################################
    print('\n--------------------------------------------------')
    print( '  RUNNING HABITAT SUITABILITY INDICES - Year %s' % year)
    print('--------------------------------------------------\n')
    os.chdir(ecohydro_dir)

    # read in Morph output file
    print(' Reading in Morphology output files to be used for HSIs.')
    
    # import grid summary file (percent land, elevations) generated by Morphology
    griddata = np.genfromtxt(new_grid_filepath,delimiter=',',skip_header=1)
    
    # landdict is a dictionary of percent land (0-100) in each 500-m grid cell, key is gridID
    landdict = dict((griddata[n][0],griddata[n][3]) for n in range(0,n500grid))
    waterdict = dict((griddata[n][0],griddata[n][5]) for n in range(0,n500grid))
    melevdict = dict((griddata[n][0],griddata[n][2]) for n in range(0,n500grid))
    wetlanddict = dict((griddata[n][0],griddata[n][4]) for n in range(0,n500grid))
    
    # Post-process Ecohydro output for HSI calculations
    print(' Reading in Ecohydro output files to be used for HSIs.')

    # import annual open water sediment deposition (mass/area) data by ICM-Hydro compartment (OW_sed_dep is 11th column in compartment_out_YYYY.csv)
    comp_summary_file = r'%s\hydro\TempFiles\compartment_out_%4d.csv' % (input_directory,year)
    comp_summary_dict = compout2dict(comp_summary_file,10)
    OWseddep_mass_dict = comp2grid(comp_summary_dict,grid_comp)
    OWseddep_depth_mm_dict = {}
    BDWaterVal = float(WM_params[31].lstrip().rstrip())
    # convert sediment deposition loading (kg/m2) to depth (mm) using bulk density of open water area (kg/m3)
    # if deposition is negative, that indicates erosion
    for n in grid_comp.keys():
        OWseddep_depth_mm_dict[n] = max(0.0,grid500_res^2*OWseddep_mass_dict[n]*1000./BDWaterVal)
    
    del(OWseddep_mass_dict,comp_summary_dict)
    
    # import annual Ecohydro output that is summarized by grid ID (Column 0 corresponds to 500m ID#, Column 7 is percent sand, and  Column 17 is average depth)    
    EH_grid_out = np.genfromtxt(EH_grid_results_filepath,delimiter=',',skip_header=1)
    depthdict = dict((EH_grid_out[n][0],EH_grid_out[n][15:17]) for n in range(0,n500grid)) # column 15 is mean summer depth, column 16 is mean annual depth
    stagedict = dict((EH_grid_out[n][0],EH_grid_out[n][12]) for n in range(0,n500grid))
    pctsanddict = dict((EH_grid_out[n][0],EH_grid_out[n][7]) for n in range(0,n500grid))

    del(EH_grid_out)

    # Calculate monthly averages from Hydro output daily timeseries 
    # read in daily timeseries and calculate averages for each ICM-Hydro compartment for a variety of variables
    
    # check length of timeseries.out files once for year and save total length of file
    daily_timeseries_file = os.path.normpath(r'%s\\STG.out' % ecohydro_dir)
    ndays_run = file_len(daily_timeseries_file)
      
    # build empty dictionaries that will be filled with monthly average values
    # key will be grid ID
    # value will be an array of 12 monthly values
    saldict = {}
    tmpdict = {}
    tssdict = {}
    ChlAdict = {}
    
    for n in grid_comp.keys():
        saldict[n] = []
        tmpdict[n] = []
        tssdict[n] = []
        ChlAdict[n]  =[]
        
    # calculate monthly averages for compartment
    print(' calculating monthly averages of ICM-Hydro output for import into HSI.py')
    for mon in range(1,13):
        data_start = dt.date(startyear,1,1)          # start date of all data included in the daily timeseries file (YYYY,M,D)
        ave_start = dt.date(year,mon,dom)          # start date of averaging window, inclusive (YYYY,M,D)
        ave_end = dt.date(year,mon,dom[mon])          # end date of averaging window, inclusive (YYYY,M,D)
        
        ##############
        # Salinity   
        ##############
        # read in daily salinity and calculate monthly mean for compartment then map to grid using daily2ave and comp2grid functions
        daily_timeseries_file = os.path.normpath(r'%s\\SAL.out' % ecohydro_dir)
        comp_month_ave_dict = daily2ave(data_start,ave_start,ave_end,daily_timeseries_file,ndays_run)
        mon_ave = comp2grid(comp_month_ave_dict,grid_comp)
        # loop through monthly average and append to array of monthly averages in dictionary to be passed into HSI.py
        for n in grid_comp.keys():
            saldict[n].append(mon_ave[n])
    
        ##############
        # Temperature 
        ##############
        # read in daily temperature and calculate monthly mean for compartment then map to grid using daily2ave and comp2grid functions
        daily_timeseries_file = os.path.normpath(r'%s\\TMP.out' % ecohydro_dir)
        comp_month_ave_dict = daily2ave(data_start,ave_start,ave_end,daily_timeseries_file,ndays_run)
        mon_ave = comp2grid(comp_month_ave_dict,grid_comp)
        # loop through monthly average and append to array of monthly averages in dictionary to be passed into HSI.py
        for n in grid_comp.keys():
            tmpdict[n].append(mon_ave[n])
    
        ##############
        # TSS 
        ##############
        # read in daily TSSand calculate monthly mean for compartment then map to grid using daily2ave and comp2grid functions
        daily_timeseries_file = os.path.normpath(r'%s\\TSS.out' % ecohydro_dir)
        comp_month_ave_dict = daily2ave(data_start,ave_start,ave_end,daily_timeseries_file,ndays_run)
        mon_ave = comp2grid(comp_month_ave_dict,grid_comp)
        # loop through monthly average and append to array of monthly averages in dictionary to be passed into HSI.py
        for n in grid_comp.keys():
            tssdict[n].append(mon_ave[n])
    
        ##########################
        # TKN as a proxy for ChlA
        ##########################
        # read in daily TKN and calculate monthly mean for compartment then map to grid using daily2ave and comp2grid functions
        daily_timeseries_file = os.path.normpath(r'%s\\TKN.out' % ecohydro_dir)
        comp_month_ave_dict = daily2ave(data_start,ave_start,ave_end,daily_timeseries_file,ndays_run)
        mon_ave = comp2grid(comp_month_ave_dict,grid_comp)
        # loop through monthly average and append to array of monthly averages in dictionary to be passed into HSI.py
        for n in grid_comp.keys():
            ChlAdict[n].append(mon_ave[n])        





    # run HSI function (run in HSI directory so output files are saved there)
    os.chdir(HSI_dir)

    # import percent edge output from geomorph routine that is summarized by grid ID
    pctedge_file = os.path.normpath('%s\\%s_N_%02d_%02d_W_pedge.csv'% (HSI_dir,runprefix,elapsedyear,elapsedyear)) # this must match name set in "WM.CalculateEcohydroAttributes" with the exception of (year) here instead of CurrentYear
    pedge = np.genfromtxt(pctedge_file,delimiter = ',',skiprows = 1)
    pctedgedict = dict((pedge[n][0],pedge[n][1]) for n in range(0,n500grid))
    del(pedge)
    
    
    try:
        HSI.HSI(gridIDs,stagedict,depthdict,melevdict,saldict,tmpdict,TSSdict,ChlAdict,veg_output_filepath,nvegtype,landdict,waterdict,pctsanddict,OWseddep_depth_mm_dict,pctedgedict,n500grid,n500rows,n500cols,yll500,xll500,year,elapsedyear,HSI_dir,WM_params,vegetation_dir,wetland_morph_dir,runprefix)
    except:
        print '******ERROR******'
        print '\n HSI model run failed - Year %s.' % year
        
import os
import ICM_HelperFunctions as hf



def HSI_AdultMenhaden(HSI_dir,csv_outprefix,asc_outprefix,gridIDs,land_mult,fresh_for_mult,bare_mult,saldict,tmpdict,wetlndict,map2grid,ascii_grid_lookup,ascii_header,ascii_header_nrows,n500rows,n500cols):
    #########################################
    ###     Adult Gulf Menhaden HSI        ##
    #########################################

    print( ' Calculating Adult Gulf Menhaden HSI')

    #months for averaging salinity and temperature data
    months2use = {3:'mar', 4:'apr', 5:'may', 6:'jun', 7:'jul', 8:'aug', 9:'sep', 10:'octb', 11:'nov'}

    gmena_gamm_gilln = {}
    gmena_gilln_file = os.path.normpath('%s/gillnet_gulfmenhaden_gamm_table_1dec.txt' % HSI_dir)
    gamm_table_delimiter = '\t' #','
    with open(gmena_gilln_file) as tf:
        nline = 0
        for line in tf: 
            if nline > 0:
                linesplit = line.split(gamm_table_delimiter)
                s = float(linesplit[0])
                t = float(linesplit[1])
                cpue_sc = float(linesplit[6])
                try:
                    gmena_gamm_gilln[s][t] = cpue_sc    # if sal is already a key in the gamm dictionary, add the temp as another key and save cpue as the value
                except:
                    gmena_gamm_gilln[s] = {}            # if sal is not already a key in gamm dictionary, add sal as key and the value will be an empty dictionary
                    gmena_gamm_gilln[s][t] = cpue_sc    # populate dictionary with temp as key and cpue as value
            nline +=1      

    HSIcsv = r'%sGMENA.csv' % csv_outprefix
    HSIasc = r'%sGMENA.asc' % asc_outprefix

    sal2use = hf.monthly_avg_dict(saldict,months2use)
    tmp2use = hf.monthly_avg_dict(tmpdict,months2use)

    with open(HSIcsv,'w') as fGMA:
        headerstring = 'GridID,HSI,sa,sa_1,ta,ta_1,v2a\n'
        fGMA.write(headerstring)

        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            sa = sal2use[gridID]
            ta = tmp2use[gridID]
            v2a = max(0.0,min(wetlndict[gridID],100.0))
            dayva = 180.1

        # truncate salinity and temperature to max values in GAMM lookup tables - temp also is truncated at a minimum value
            sa_1 = min(sa,36.8)
            ta_1 = max(3.4,min(ta,35.9))
            
        # use sal & temp to lookup scaled CPUE value from GAMM lookup table (imported above)
        # if sal/temp combination does not exist in lookup table, set term to error flag which is used later to skip HSI calculation
            
            try:
                S1a = gmena_gamm_gilln[round(sa_1,1)][round(ta_1,1)]
                if S1a < 0.0:
                    S1a = 0.0
            except:
                S1a = -9999  
                
            if v2a <= 30.:
                S2a = 1.
            else:
                S2a = 1.43-0.0143*v2a
            
        # check for error in imported GAMM lookup table values - if sal/temp combination was not in table (or there is a precision mismatch), do not calculate HSI and report out error flag of -9999
            if S1a == -9999:    
                HSI_adltMenh = -9999
            else:
                HSI_adltMenh = zero_mult*(S1a*S2a)**(1./2.)
            
            writestring = '%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_adltMenh,sa,sa_1,ta,ta_1,v2a)
            fGMA.write(writestring)

    # NOTE this was commented out; map2grid flag added (2029MP_Task06)
    hf.HSIascii_grid(HSIcsv,HSIasc,map2grid,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    # delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(sal_MarNov_ave,tmp_MarNov_ave)

    # delete temporary variables so they do not accidentally get used in other HSIs
    # del(sj,sj_1,tj,tj_1,dayvj,v2j,gmenj_gamm_seine,S1j,S2j,sa,sa_1,ta,ta_1,dayva,v2a,gmena_gamm_gilln,S1a,S2a)
    del(sa,sa_1,ta,ta_1,dayva,v2a,gmena_gamm_gilln,S1a,S2a)
    
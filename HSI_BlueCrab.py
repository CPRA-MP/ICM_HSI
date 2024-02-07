import os
import ICM_HelperFunctions as hf



def HSI_BlueCrab(HSI_dir,csv_outprefix,asc_outprefix,gridIDs,land_mult,fresh_for_mult,bare_mult,saldict,tmpdict,wetlndict,watsavdict,cultchdict,map2grid,ascii_grid_lookup,ascii_header,ascii_header_nrows,n500rows,n500cols):
    ########################################
    ##       Juvenile Blue Crab HSI       ##
    ########################################
    print( ' Calculating Blue Crab HSI')
    
    months2use = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'octb', 'nov', 'dec']

    # read in GAMM lookup table for sal/temp combinations
    blucj_gamm_seine = {}
    blucj_seine_file = os.path.normpath('%s/seine_bluecrab_gamm_table_1dec.txt' % HSI_dir)
    gamm_table_delimiter = '\t' #','
    with open(blucj_seine_file) as tf:
        nline = 0
        for line in tf: 
            if nline > 0:
                linesplit = line.split(gamm_table_delimiter)
                s = float(linesplit[0])
                t = float(linesplit[1])
                cpue_sc = float(linesplit[6])
                try:
                    blucj_gamm_seine[s][t] = cpue_sc    # if sal is already a key in the gamm dictionary, add the temp as another key and save cpue as the value
                except:
                    blucj_gamm_seine[s] = {}            # if sal is not already a key in gamm dictionary, add sal as key and the value will be an empty dictionary
                    blucj_gamm_seine[s][t] = cpue_sc    # populate dictionary with temp as key and cpue as value
            nline +=1        

    HSIcsv = r'%sBLUCJ.csv' % csv_outprefix
    HSIasc = r'%sBLUCJ.asc' % asc_outprefix

    sal2use = hf.monthly_avg_dict(saldict,months2use)
    tmp2use = hf.monthly_avg_dict(tmpdict,months2use)

    with open(HSIcsv,'w') as fBC:
        
        headerstring = 'gridID,HSI,s,s_1,t,t_1,v2,oysc,savc,S1\n'
        fBC.write(headerstring)

        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            
            s = sal2use[gridID]
            t = tmp2use[gridID]
            v2 =  max(0.0,min(wetlndict[gridID],100.0))
            savc = max(0.0,min(watsavdict[gridID],100.0))  
            oysc = max(0.0,min(cultchdict[gridID],1.0))
            dayv = 6.56

            # truncate salinity and temperature to max values in GAMM lookup tables - temp also is truncated at a minimum value
            s_1 = min(s,36.8)
            t_1 = max(3.2,min(t,35.2))
                
            # use sal & temp to lookup scaled CPUE value from GAMM lookup table (imported above)
            # if sal/temp combination does not exist in lookup table, set term to error flag which is used later to skip HSI calculation
            try:
                S1 = blucj_gamm_seine[round(s_1,1)][round(t_1,1)]        # this will lookup the scaled cpue value from the imported blcrab seine GAMM lookup table that has precision to the tenths place #.#
                if S1 < 0.0:
                    S1 = 0.0
            except:
                S1 = -9999 
                
            if v2 < 25.0:
                S2 = 0.03*v2+0.25         # note three different functions for when v2s is less than 25.
                if oysc >= 0.5:           # if oyster HSI greater than 0.5 or sav cover greater than 20. then use different S2s function
                    S2 = 0.02*v2+0.5     
                if savc >= 20.:             
                    S2 = 0.008*v2+0.8
            elif v2 >= 25.0 and v2 <= 80.:
                S2 = 1.0
            else:
                S2 = max(0.0,(5.-0.05*v2))

            # check for error in imported GAMM lookup table values - if sal/temp combination was not in table (or there is a precision mismatch), do not calculate HSI and report out error flag of -9999
            if S1 == -9999:
                HSI_juvBlueCrab = -9999
            else:
                HSI_juvBlueCrab = zero_mult*max(0.0,(S1*S2))**(1./2.)

            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_juvBlueCrab,s,s_1,t,t_1,v2,oysc,savc,S1)
            fBC.write(writestring)

    # NOTE this was commented out; map2grid flag added (2029MP_Task06)
    hf.HSIascii_grid(HSIcsv,HSIasc,map2grid,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    # delete temporary variables so they do not accidentally get used in other HSIs
    del(s,s_1,t,t_1,v2,oysc,savc,dayv,blucj_gamm_seine,S1,S2)
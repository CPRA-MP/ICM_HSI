def HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header):
    import numpy as np
    import arcpy
    print ' - mapping HSI to ASCII grid'    
# read HSI csv file into numpy array - usecol = column of csv file that has HSI value
    newHSI = np.genfromtxt(HSIcsv,delimiter=',',usecols=[0,1],skiprows=1)
    newHSIdict = dict((newHSI[n][0],newHSI[n][1])for n in range(0,len(newHSI)))
    # prepare zero array in same shape of original Veg output ASCII grid
    newHSIgrid=np.zeros([n500rows,n500cols])
                    
    for m in range(0,n500rows):
        for n in range(0,n500cols):
            cellID = ascii_grid_lookup[m][n]
            if cellID == -9999:
                newHSIgrid[m][n] = -9999
            else:
                try:
                    newHSIval = newHSIdict[cellID] 
                    if np.isnan(newHSIval):
                        newHSIgrid[m][n] = -1.0
                    elif np.isinf(newHSIval):
                        newHSIgrid[m][n] = -1.0
                    else:
                        newHSIgrid[m][n] = newHSIval
                except:   # if cellID is not a key in the newLULCdictionay - assign cell to NoData
                    newHSIgrid[m][n] = -9999
    print " - saving new HSI ASCII raster file"
    # save formatted LULC grid to ascii file with appropriate ASCII raster header
    np.savetxt(HSIasc,newHSIgrid,fmt='%.2f',delimiter=' ',header=ascii_header,comments='')
    
    newHSI = 0
    newHSIdict = {}
    
    
def HSI(gridIDs,stagedict,depthdict,melevdict,saldict,tmpdict,TSSdict,ChlAdict,veg_output_filepath,nvegtype,landdict,waterdict,pctsanddict,pctedgedict,n500grid,n500rows,n500cols,yll500,xll500,year,elapsedyear,HSI_dir,WM_params,vegetation_dir,wetland_morph_dir,runprefix):

    import numpy as np
    from dbfpy import dbf
    import os
    import csv
    import code

# set some general variables
    print ' Setting up HSI runs.'

    asc_outprefix = r'%s\\output_%02d\\Deliverables\\%s_O_%02d_%02d_X_' % (wetland_morph_dir,elapsedyear,runprefix,elapsedyear,elapsedyear)
    csv_outprefix = r'%s\\%s_O_%02d_%02d_X_' % (HSI_dir,runprefix,elapsedyear,elapsedyear)
    

    e = 2.718281828
    jan,feb,mar,apr,may,jun,jul,aug,sep,octb,nov,dec = 0,1,2,3,4,5,6,7,8,9,10,11
    ChlA = ChlAdict.keys()
    sal = saldict.keys()
    tmp = tmpdict.keys()
    TSS = TSSdict.keys()
    dept = depthdict.keys()
    
#    veg_ascii_grid = WM_params[61].lstrip().rstrip()
#    veg_grid_ascii_file = os.path.normpath(vegetation_dir + '\\' + veg_ascii_grid)    
    grid_ascii_file = os.path.normpath(HSI_dir + '\\hsi_grid_ecoregion.asc')    
    print ' Reading in ASCII grid template.'

    ascii_grid_lookup = np.genfromtxt(grid_ascii_file,delimiter=' ',skiprows=6)
    ascii_header='nrows %s \nncols %s \nyllcorner %s \nxllcorner %s \ncellsize 500.0 \nnodata_value -9999.00' % (n500rows,n500cols,yll500,xll500)
    
    
    
# generate some dictionaries from the Hydro output files - these combine the monthly output from hydro into mean values for various time frames (e.g. annual, April-July, etc)
# other input values that are used by specific HSIs will be written and deleted within the respective HSIs to minimize memory requirements
# GridIDs from Ecohydro output file - this could be in a different order than n500grid, so use it as the dictionary key
    ChlA_JanDec_ave = dict( (ChlA[n-1],1000.0*np.mean([ChlAdict[n][jan],ChlAdict[n][feb],ChlAdict[n][mar],ChlAdict[n][apr],ChlAdict[n][may],ChlAdict[n][jun],ChlAdict[n][jul],ChlAdict[n][aug],ChlAdict[n][sep],ChlAdict[n][octb],ChlAdict[n][nov],ChlAdict[n][dec]])) for n in range(1,n500grid+1))
    sal_JanDec_ave = dict((sal[n-1],np.mean([saldict[n][jan],saldict[n][feb],saldict[n][mar],saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb],saldict[n][nov],saldict[n][dec]]))for n in range(1,n500grid+1))
    tmp_JanDec_ave = dict((tmp[n-1],np.mean([tmpdict[n][jan],tmpdict[n][feb],tmpdict[n][mar],tmpdict[n][apr],tmpdict[n][may],tmpdict[n][jun],tmpdict[n][jul],tmpdict[n][aug],tmpdict[n][sep],tmpdict[n][octb],tmpdict[n][nov],tmpdict[n][dec]]))for n in range(1,n500grid+1))
    sal_AprJul_ave = dict((sal[n-1],np.mean([saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul]]))for n in range(1,n500grid+1))
    dep_JanDec_ave = dict((dept[n-1],depthdict[n][1])for n in range(1,n500grid+1))
    dep_JulSep_ave = dict((dept[n-1],depthdict[n][0]) for n in range(1,n500grid+1)) #this is summertime depth value from hydro output
    dep_OctJun_ave = dict((dept[n-1],(0.2*depthdict[n][0]+0.8*depthdict[n][1])) for n in range(1,n500grid+1))

    
    #dep_JanDec_ave = dict((dept[n-1],np.mean([depthdict[n][jan],depthdict[n][feb],depthdict[n][mar],depthdict[n][apr],depthdict[n][may],depthdict[n][jun],depthdict[n][jul],depthdict[n][aug],depthdict[n][sep],depthdict[n][octb],depthdict[n][nov],depthdict[n][dec]]))for n in range(1,n500grid+1))
    #dep_OctJun_ave = dict((dept[n-1],np.mean([depthdict[n][jan],depthdict[n][feb],depthdict[n][mar],depthdict[n][apr],depthdict[n][may],depthdict[n][jun],depthdict[n][octb],depthdict[n][nov],depthdict[n][dec]]))for n in range(1,n500grid+1))
    #dep_JulSep_ave = dict((dept[n-1],np.mean([depthdict[n][jul],depthdict[n][aug],depthdict[n][sep]]))for n in range(1,n500grid+1))
    
    
# read in Veg output file - this is the same code that is used in WM.ImportVegResults()
    print ' Reading in LAVegMod output files to be used for HSIs.'
        
    # skipvalue is the number of rows contained in the header and the grid array located at the start of the Veg output file
    skipvalue = n500rows + 7
    
    # generate zeros array that will be filled with Veg results
    vegcolumns = nvegtype + 1   #veg columns is the number of vegetation types modeled plus one for the grid ID column
    new_veg = np.zeros((n500grid,vegcolumns))
    veg_missing = 0
    # open Vegetation output file
    with open(veg_output_filepath,'r') as vegfile:
    # skip ASCII header rows and ASCII grid at start of output file
        for n in range(0,skipvalue-1):
            dump=vegfile.readline()
    # read in header of Vegetation output at bottom of ASCII grid    
        vegtypenames = vegfile.readline().split(',')
    # remove any leading or trailing spaces in veg types
        for n in range(0,len(vegtypenames)):
            vegtypenames[n] = vegtypenames[n].lstrip().rstrip()
    # loop through rest of Vegetation file         
        for nn in range(0,n500grid):
    # split each line of file based on comma delimiter (any spaces will be removed with rstrip,lstrip)
            vline = vegfile.readline().split(",")
    # if all columns have data in output file, assign veg output to veg_ratios array
            if (len(vline) == vegcolumns):
                for nnn in range(0,len(vline)):
                    new_veg[nn,nnn]=float(vline[nnn].lstrip().rstrip())
    # if there are missing columns in line, set first column equal to grid cell, and set all other columns equal to 0.
            else:
                for nnn in range(1,vegcolumns):
                    new_veg[nn,0]=nn+1
                    new_veg[nn,nnn] = 0.0
                veg_missing += 1
    if (veg_missing > 0):
        print ( ' Some Vegetation output was not written correctly to Veg output file.')
        print ('  - %s 500m grid cells did not have complete results in Veg Output file.' % veg_missing)
    
    # read lookup table into dictionary that converts output VegType to new LULC value
    LULC_Lookup = r'%s\\%s' % (wetland_morph_dir,WM_params[32].lstrip().rstrip())           #this is the same code as in WM.main()
    LULC_OldValField = WM_params[33].lstrip().rstrip()      #this is the same code as in WM.main()
    LULC_NewValField = WM_params[34].lstrip().rstrip()      #this is the same code as in WM.main()

    VegLULC_lookup={}
    with open(LULC_Lookup, mode='r') as infile:
        for i,row in enumerate(infile):
            if i == 0:
                hds = row.split(',')
                for n in range(0,len(hds)):
                    hds[n] = hds[n].lstrip().rstrip()   #remove any leading or trailing spaces
                old = hds.index(LULC_OldValField)
                new = hds.index(LULC_NewValField)
            else:
                VegLULC_lookup[row.split(',')[old].rstrip().lstrip()]=row.split(',')[new].rstrip().lstrip()    

    print ' Reclassifying Veg species output into general LULC types used by HSI equations.'
# generate some blank dictionaries that will be filled with Veg output
#    waterdict = {}
    wetlndict = {}
    frattdict = {}
    frfltdict = {}
    interdict = {}
    brackdict = {}
    salmardict = {}
    swfordict = {}
    btfordict = {}
    watsavdict = {}
    baldcypdict = {}
    blackmangrovedict = {}
    marshelderdict = {}
    baredict = {}
    bare_mult = {}
    fresh_for_mult = {}
    land_mult = {}

    # determine portion of cell that is covered by water, land, and different wetland types
    for n in range(0,len(new_veg)):
        gridID = new_veg[n][0]
        
        # use landdict to assign portion of cell that is water - this value is the updated land/water ratio AFTER the Morph run, 'WATER' value from Veg output is not needed here
        # Floating marsh Veg output was used to update land/water in WM.ImportVegResults()
        # Dead floating marsh output is classified as water in landdict
        # Live floating marsh output is included in the total percent land in landdict
       
#        try:
            # check that percent land is Data (-9999 if NoData), if NoData, set water area to zero
#            if landdict[gridID] >= 0:
#                waterdict[gridID] = 100 - landdict[gridID]
#            else:
#                waterdict[gridID] = 0
#        except:
#            waterdict[gridID] = 0
        
        # set initial dictionary values to zero for current gridID
        wetlndict[gridID] = 0.0     # percent of cell that has wetland (all non-floating types)
        frattdict[gridID] = 0.0     # percent of cell that has fresh wetland (non-floating types)
        frfltdict[gridID] = 0.0     # percent of cell that has fresh floating marsh
        interdict[gridID] = 0.0     # percent of cell that has intermediate wetland
        brackdict[gridID] = 0.0     # percent of cell that has brackish wetland
        salmardict[gridID] = 0.0    # percent of cell that has saline wetland
        swfordict[gridID] = 0.0     # percent of cell that has swamp forest wetland (same as swamp forest - LULC reclass doesn't differentiate between the two)
        btfordict[gridID] = 0.0     # percent of cell that has bottomland  forest (same as swamp forest - LULC reclass doesn't differentiate between the two)
        watsavdict[gridID] = 0.0    # percent of cell that has subaquatic vegetation
        baldcypdict[gridID] = 0.0 
        marshelderdict[gridID] = 0.0
        blackmangrovedict[gridID] = 0.0
        baredict[gridID] = 0.0
        bare_mult[gridID] = 1.0
        fresh_for_mult[gridID] = 1.0
        
        # set land multiplier to zero for grid cells that are 100% land
        if waterdict[gridID] == 0.0:
            land_mult[gridID] = 0.0
        else:
            land_mult[gridID] = 1.0
        
       
        # loop through Veg output and add percent cover to the respective wetland type dictionaries
        # veg output is portion of cell covered by species (0-1) HSI equations use percentages - therefore multiply Veg data by 100
        for k in VegLULC_lookup.keys():         # k is veg type as dictionary key in VegLULC_lookup
            vc = vegtypenames.index(k)          # vc is column in new_veg array that matches veg type, k
            nlulc = int(VegLULC_lookup[k])      # nlulc is new LULC number to convert veg type,k, t
            
            # Fresh forested LULC
            if nlulc == 1:
                #wetlndict[gridID] += new_veg[n][vc] #do not consider fresh forested as a 'wetland' type
                btfordict[gridID] += new_veg[n][vc]
                swfordict[gridID] += new_veg[n][vc]
                
            # Fresh herbaceous wetland LULC
            if nlulc == 2:
                wetlndict[gridID] += new_veg[n][vc]
                frattdict[gridID] += new_veg[n][vc]
            
            # Intermediate herbaceous wetland LULC    
            if nlulc == 3:
                wetlndict[gridID] += new_veg[n][vc]
                interdict[gridID] += new_veg[n][vc]

            # Brackish herbaceous wetland LULC                
            if nlulc == 4:
                wetlndict[gridID] += new_veg[n][vc]
                brackdict[gridID] += new_veg[n][vc]
            # Saline herbaceous wetland LULC                
            if nlulc == 5:
                wetlndict[gridID] += new_veg[n][vc]
                salmardict[gridID] += new_veg[n][vc]
            
            # Subaqauatic vegetation LULC
            if vc == vegtypenames.index('SAV'):
                watsavdict[gridID] += new_veg[n][vc]

            # Live floating marsh LULC                
            if nlulc == 8:              
                frfltdict[gridID] += new_veg[n][vc]
                wetlndict[gridID] += frfltdict[gridID]      # add floating marsh percent cover to the total wetland percent cover

            if nlulc == 9:
                baredict[gridID] += new_veg[n][vc]
                

            # Bald cypress
            if vc == vegtypenames.index('TADI2'): 
                baldcypdict[gridID] = new_veg[n][vc]
                
            # Black mangrove
            if vc == vegtypenames.index('AVGE'): 
                blackmangrovedict[gridID] = new_veg[n][vc]
                
            # Marsh elder
            if vc == vegtypenames.index('IVFR'):
                marshelderdict[gridID] = new_veg[n][vc]
# Check for bareground                 
# if there is no wetland or forest type, but there is bareground, set bareground multiplier to zero

        if baredict[gridID] > 0.0:
            if wetlndict[gridID] == 0.0:
                if btfordict[gridID] == 0.0:
                    if watsavdict[gridID] == 0.0:
                        bare_mult[gridID] = 0.0
# if there is wetland area, and it is greater than forested area, add bareground to wetland area
            elif wetlndict[gridID] > btfordict[gridID]:
                wetlndict[gridID] += baredict[gridID]
           # if forest is greater than wetland area, add bareground to foreseted areas (both swamp forest and bottom hardwood - since they are set equal)
            else:
                btfordict[gridID] += baredict[gridID]
                swfordict[gridID] += baredict[gridID]

# if fresh forest is present and greater than wetland area, set fresh forest multiplier to zero
        if btfordict[gridID] > 0.0:
            if btfordict[gridID] > wetlndict[gridID]:
                fresh_for_mult[gridID] = 0.0


# convert marsh/land type dictionaries from ratio to percentage    
    for gridID in gridIDs:
        wetlndict[gridID] =  max(0.0,min(100.0,100.0*wetlndict[gridID]))
        btfordict[gridID] =  max(0.0,min(100.0,100.0*btfordict[gridID]))
        swfordict[gridID] =  max(0.0,min(100.0,100.0*swfordict[gridID]))
        frattdict[gridID] =  max(0.0,min(100.0,100.0*frattdict[gridID]))
        interdict[gridID] =  max(0.0,min(100.0,100.0*interdict[gridID]))
        brackdict[gridID] =  max(0.0,min(100.0,100.0*brackdict[gridID]))
        salmardict[gridID] = max(0.0,min(100.0,100.0*salmardict[gridID]))
        watsavdict[gridID] = max(0.0,min(100.0,100.0*watsavdict[gridID]))
        frfltdict[gridID] =  max(0.0,min(100.0,100.0*frfltdict[gridID]))
        baldcypdict[gridID]= max(0.0,min(100.0,100.0*baldcypdict[gridID]))

   ########################################
    ##       Juvenile Blue Crab HSI       ##
    ########################################
    print ' Calculating Blue Crab HSI'
    
    HSIcsv = r'%sBLUCJ.csv' % csv_outprefix
    HSIasc = r'%sBLUCJ.asc' % asc_outprefix
    

    # save input values that are only used in this HSI
    sal_JanMarAugDec_ave = dict((sal[n-1],np.mean([saldict[n][jan],saldict[n][feb],saldict[n][mar],saldict[n][aug],saldict[n][sep],saldict[n][octb],saldict[n][nov],saldict[n][dec]]))for n in range(1,n500grid+1))
    tmp_JanMarAugDec_ave = dict((tmp[n-1],np.mean([tmpdict[n][jan],tmpdict[n][feb],tmpdict[n][mar],tmpdict[n][aug],tmpdict[n][sep],tmpdict[n][octb],tmpdict[n][nov],tmpdict[n][dec]]))for n in range(1,n500grid+1))

    with open(HSIcsv,'w') as fBC:
        
        headerstring = 'GridID,HSI,s,t,v2,\n'
        fBC.write(headerstring)

        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            
            s = sal_JanMarAugDec_ave[gridID]
            t = tmp_JanMarAugDec_ave[gridID]
            v2 =  max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))
            dayv = 2.088
            if s < 35.0:
                if t < 35.0:
                    lnCPUE1=0.8587-0.2451*dayv+0.07012*dayv**2.-0.03677*s+0.06561*t+0.000312*s**2.-0.00182*t**2.
                    S1 =(e**lnCPUE1 - 1.)/2.47
                else:
                    S1 = 0.0
            else:
                S1 = 0.0
                
            if v2 < 25.:                       
                S2 = 0.028*v2+0.3                     
            elif v2 <= 80.:                       
                S2 = 1.                       
            else:                       
                S2 = max(0,(5.-0.05*v2))
            
            S2 = S2
            
            HSI_juvBlueCrab = zero_mult*max(0,(S1*S2))**(1./2.)

#           writestring = '%s,%s\n' % (gridID,HSI_juvBlueCrab)
            writestring = '%s,%s,%s,%s,%s\n' %(gridID,HSI_juvBlueCrab,s,t,v2)
            fBC.write(writestring)
    
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
    
    # delete temporary variables so they do not accidentally get used in other HSIs
    del(s,t,v2,dayv,lnCPUE1,S1,S2)

    ########################################
    ##        Largemouth Bass HSI         ##
    ########################################
    print ' Calculating Largemouth Bass HSI'

    sal_MarNov_ave = dict((sal[n-1],np.mean([saldict[n][mar],saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb],saldict[n][nov]]))for n in range(1,n500grid+1))
    tmp_MarNov_ave = dict((tmp[n-1],np.mean([tmpdict[n][mar],tmpdict[n][apr],tmpdict[n][may],tmpdict[n][jun],tmpdict[n][jul],tmpdict[n][aug],tmpdict[n][sep],tmpdict[n][octb],tmpdict[n][nov]]))for n in range(1,n500grid+1))
    TSS_MarNov_ave = dict((TSS[n-1],np.mean([TSSdict[n][mar],TSSdict[n][apr],TSSdict[n][may],TSSdict[n][jun],TSSdict[n][jul],TSSdict[n][aug],TSSdict[n][sep],TSSdict[n][octb],TSSdict[n][nov]]))for n in range(1,n500grid+1))
    
    HSIcsv = r'%sLMBAS.csv' % csv_outprefix
    HSIasc = r'%sLMBAS.asc' % asc_outprefix
    
    
    with open(HSIcsv,'w') as fLMB:

        headerstring = 'GridID,HSI,s,t,tss,ntu,v2,c,S1,S2,S3,lmb_zero_mult\n'
        fLMB.write(headerstring)

        for gridID in gridIDs:
            # zero multiplier for LMB will set all HSIs for grid cells 
            lmb_zero_mult = land_mult[gridID]*bare_mult[gridID]
            s = sal_MarNov_ave[gridID]
            t = tmp_MarNov_ave[gridID]
            tss = TSS_MarNov_ave[gridID]
            v2 = max(0.0,min(wetlndict[gridID] + watsavdict[gridID] + btfordict[gridID] ,100.0))
            c = ChlA_JanDec_ave[gridID]


            dayv = 1.99
            dayv2 = 4.808
            # convert TSS to turbidity - current conversion based median conversion from Lit review (
            # tss = m*NTU + b
            #m = 0.8
            #b = -7.37
            #ntu = max((tss - b)/m,0) # removed from equation #Edw - 8/3/2015
            ntu = 0.0
            
            lnCPUE1 = 0.8752-1.7125*dayv+0.3768*dayv2-0.2759*s+0.007203*s**2.+0.3328*t-.0406*ntu-0.00764*t**2.+0.000632*ntu**2.

            if s <= 20:
                S1 = (e**lnCPUE1 - 1.)/12.214
            else:
                S1 = 0.0

            if v2 < 20.:
                S2 = 0.01
            elif v2 < 30.:
                S2 = 0.099*v2-1.97
            elif v2 < 50.:
                S2 = 1.
            elif v2 < 85.:
                S2 = -0.0283*v2+2.414
            elif v2 < 100.:
                S2 = 0.01
            else:
                S2 = 0.

            sdenom = (1+e**(-(c-35.75)/7.8414))
            if sdenom <> 0.0:
                S3 = 0.24+0.85/sdenom
            else:
                S3 = 0.0

# ChlA not currently  calibrated, hardset S3 to 1.0 and c to -9999 for write statements
            c = -9999.0
            S3 = 1.0

            HSI_LMBass = lmb_zero_mult*(S1*S2*S3)**(1./3.)

            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_LMBass,s,t,tss,ntu,v2,c,S1,S2,S3,lmb_zero_mult)

            fLMB.write(writestring)
    
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)            
            
    # delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(sal_MarNov_ave,tmp_MarNov_ave,TSS_MarNov_ave)
    
    # delete temporary variables so they do not accidentally get used in other HSIs
    del(s,t,tss,c,dayv,v2,lnCPUE1,S1,S2,S3)

########################################
##          Bay Anchovy HSIs          ##
########################################
##      Juvenile Bay Anchovy HSI      ##
########################################
    print ' Calculating Juvenile and Adult Bay Anchovy HSIs'

    HSIcsv = r'%sBAYAJ.csv' % csv_outprefix
    HSIasc = r'%sBAYAJ.asc' % asc_outprefix

    HSIcsv2 = r'%sBAYAA.csv' % csv_outprefix
    HSIasc2 = r'%sBAYAA.asc' % asc_outprefix


    with open(HSIcsv,'w') as fBAJ, open(HSIcsv2,'w') as fBAA:
        
        headerstring = 'GridID,HSI,s,t,c,v2\n'
        fBAJ.write(headerstring)
        fBAA.write(headerstring)
        
        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            s = sal_JanDec_ave[gridID]
            t = tmp_JanDec_ave[gridID]
            c = ChlA_JanDec_ave[gridID]
            v2j =  max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))
            v2a = v2j
            dayv = 2.
# Calculate variables for Juvenile Bay Anchovy
            lnCPUE1j = -2.6496+0.8946*dayv-0.1896*dayv**2.-0.00678*s+0.4324*t-0.0003*s**2.+0.000008*(s*t)**2.-0.00023*t*s**2.-0.00924*t**2.
            S1j = (e**lnCPUE1j - 1.)/30.14

            if v2j < 25:
                S2j = 0.028*v2j+0.3
            elif v2j <= 80.:
                S2j = 1.
            else:
                S2j = 5.-0.05*v2j

            S3j = 4.18*e**(-4.59*e**(-0.02*c))/3.82

# ChlA not currently  calibrated, hardset S3 to 1.0 and c to -9999 for write statements
            c = -9999.0
            S3j = 1.0

            HSI_juvAnch = zero_mult*(S1j*S2j*S3j)**(1./3.)

            writestring = '%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_juvAnch,s,t,c,v2j)
            fBAJ.write(writestring)
            
########################################
##       Adult Bay Anchovy HSI        ##
########################################
# Calculate variables for Adult Bay Anchovy
            lnCPUE1a = 4.3195-0.363*s-0.3057*t+0.0108*s**2.+0.00872*t**2.+0.0632*s*t+0.000045*(s*t)**2.-0.00162*s*t**2.-0.00183*s**2.*t

            S1a = (e**lnCPUE1a - 1.)/42.92

            if v2a <= 30.:
                S2a = 1.
            else:
                S2a = 1.43-0.0143*v2a

            S3a = S3j

            HSI_adltAnch = zero_mult*(S1a*S2a*S3a)**(1./3.)

            writestring = '%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_adltAnch,s,t,c,v2a)
            fBAA.write(writestring)

# map juvenile bay anchovy HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
# map adult bay anchovy HSI to Ascii grid
    HSIascii_grid(HSIcsv2,HSIasc2,ascii_grid_lookup,n500cols,n500rows,ascii_header)
    
# delete temporary variables so they do not accidentally get used in other HSIs
    del(s,t,c,dayv,v2j,lnCPUE1j,S1j,S2j,S3j,v2a,lnCPUE1a,S1a,S2a,S3a)


########################################
##         Gulf Menhaden HSIs         ##
########################################
##      Juvenile Gulf Menhaden HSI    ##
########################################
    print ' Calculating Juvenile and Adult Gulf Menhaden HSIs'
# save input values that are only used in this HSI
    sal_JanJul_ave = dict((sal[n-1],np.mean([saldict[n][jan],saldict[n][feb],saldict[n][mar],saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul]]))for n in range(1,n500grid+1))
    tmp_JanJul_ave = dict((tmp[n-1],np.mean([tmpdict[n][jan],tmpdict[n][feb],tmpdict[n][mar],tmpdict[n][apr],tmpdict[n][may],tmpdict[n][jun],tmpdict[n][jul]]))for n in range(1,n500grid+1))
    sal_MarOct_ave = dict((sal[n-1],np.mean([saldict[n][mar],saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb]]))for n in range(1,n500grid+1))
    tmp_MarOct_ave = dict((tmp[n-1],np.mean([tmpdict[n][mar],tmpdict[n][apr],tmpdict[n][may],tmpdict[n][jun],tmpdict[n][jul],tmpdict[n][aug],tmpdict[n][sep],tmpdict[n][octb]]))for n in range(1,n500grid+1))

    HSIcsv = r'%sGMENJ.csv' % csv_outprefix
    HSIasc = r'%sGMENJ.asc' % asc_outprefix

    HSIcsv2 = r'%sGMENA.csv' % csv_outprefix
    HSIasc2 = r'%sGMENA.asc' % asc_outprefix

    with open(HSIcsv,'w') as fGMJ, open(HSIcsv2,'w') as fGMA:

        headerstring = 'GridID,HSI,s,t,c,v2\n'
        fGMJ.write(headerstring)
        fGMA.write(headerstring)

        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            sj = sal_JanJul_ave[gridID]
            tj = tmp_JanJul_ave[gridID]
            c = ChlA_JanDec_ave[gridID]
            v2j = max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))
            dayvj = 1.04

# Calculate variables for Juvenile Gulf Menhaden
            lnCPUE1j = -0.4572+2.6189*dayvj-1.3848*dayvj**2.-0.06918*sj+0.1778*tj-0.00331*tj**2.

            S1j = (e**lnCPUE1j - 1.)/22.48

            if v2j < 25.:
                S2j = 0.028*v2j+0.3
            elif v2j <= 80.:
                S2j = 1.
            else:
                S2j = 5.-0.05*v2j

            S3j = 4.18*e**(-4.59*e**(-0.02*c))/3.82
            
# ChlA not currently  calibrated, hardset S3 to 1.0 and c to -9999 for write statements
            c = -9999.0
            S3j = 1.0

            

            HSI_juvMenh = zero_mult*(S1j*S2j*S3j)**(1./3.)
            
            writestring = '%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_juvMenh,sj,tj,c,v2j)
            fGMJ.write(writestring)
########################################
##     Adult Gulf Menhaden HSI        ##
########################################
# Calculate variables for Adult Gulf Menhaden
            sa = sal_MarOct_ave[gridID]
            ta = tmp_MarOct_ave[gridID]
            v2a = v2j
            dayva = 1.8
            lnCPUE1a = -0.9567+0.3062*dayva-0.1123*dayva**2.+0.01829*sa+0.1109*ta-0.00018*sa**2.-0.00008*ta*sa**2.-0.00263*ta**2.+0.000112*sa*ta**2.

            S1a = (e**lnCPUE1a - 1.)/5.81

            if v2a <= 30.:
                S2a = 1.
            else:
                S2a = 1.43-0.0143*v2a

            S3a = S3j

            HSI_adltMenh = zero_mult*(S1a*S2a*S3a)**(1./3.)
            writestring = '%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_adltMenh,sa,ta,c,v2a)
            fGMA.write(writestring)

# map juvenile menhaden HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
# map adult menhaden HSI to Ascii grid
    HSIascii_grid(HSIcsv2,HSIasc2,ascii_grid_lookup,n500cols,n500rows,ascii_header)

# delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(sal_JanJul_ave,tmp_JanJul_ave,sal_MarOct_ave,tmp_MarOct_ave)
    
# delete temporary variables so they do not accidentally get used in other HSIs
    del(tj,c,dayvj,v2j,lnCPUE1j,S1j,S2j,S3j,sa,ta,dayva,v2a,lnCPUE1a,S1a,S2a,S3a)


##############################################
##          Brown Shrimp HSIs               ##
##############################################

    print ' Calculating Brown Shrimp HSIs'
# save input values that are only used in this HSI
    sal_AprJun_ave = dict((sal[n-1],np.mean([saldict[n][apr],saldict[n][may],saldict[n][jun]]))for n in range(1,n500grid+1))
    tmp_AprJun_ave = dict((tmp[n-1],np.mean([tmpdict[n][apr],tmpdict[n][may],tmpdict[n][jun]]))for n in range(1,n500grid+1))
    tmp_AprJul_ave = dict((tmp[n-1],np.mean([tmpdict[n][apr],tmpdict[n][may],tmpdict[n][jun],tmpdict[n][jul]]))for n in range(1,n500grid+1))

    HSIcsv = r'%sBSHRS.csv' % csv_outprefix
    HSIasc = r'%sBSHRS.asc' % asc_outprefix

    HSIcsv2 = r'%sBSHRL.csv' % csv_outprefix
    HSIasc2 = r'%sBSHRL.asc' % asc_outprefix


    with open(HSIcsv,'w') as fBSS, open(HSIcsv2,'w') as fBST:
        
        headerstring = 'GridID,HSI,s,t,v2\n'

        fBSS.write(headerstring)
        fBST.write(headerstring)

        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            ss = sal_AprJun_ave[gridID]
            ts = tmp_AprJun_ave[gridID]
            v2s =  max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))
            dayvs = 1.35

##############################################
##  Small Juvenile Brown Shrimp HSI - Seine ##
##############################################            
            lnCPUE1s = -15.9328+24.9838*dayvs-9.0311*dayvs**2.+0.2203*ss+0.02229*ts-0.00629*ss**2.+0.000544*ts**2.-0.00007*ss*ts**2.
            S1s = (e**lnCPUE1s - 1.)/32.17

            if v2s < 25.:
                S2s = 0.028*v2s+0.3
            elif v2s <= 80.:
                S2s = 1.
            else:
                S2s = 5.-0.05*v2s

            HSI_BrShrSeine = zero_mult*max(0,(S1s*S2s))**(1./2.)

            writestring = '%s,%s,%s,%s,%s\n' %(gridID,HSI_BrShrSeine,ss,ts,v2s)
            fBSS.write(writestring)

#################################################
##     Large Juvenile Brown Shrimp HSI - Trawl ##
#################################################


            st = sal_AprJul_ave[gridID]
            tt = tmp_AprJul_ave[gridID]
            v2t = v2s
            dayvt = 1.4578
            gear =  0.0
            
            lnCPUE1t = -8.931-0.1434*st-0.1801*tt+0.003639*st**2.+0.006205*tt**2.+0.04524*st*tt+0.000034*st**2.*tt**2.-0.00126*st*tt**2.-0.00125*tt*st**2.+15.973*dayvt-5.3793*dayvt**2.+0.0676*gear
            
            S1t = (e**lnCPUE1t - 1.)/50.55

# salinity conditional filter on S1a
            smult = 1.0
#            if st < 4.:
#                smult = 0.
#            else:
#                smult = 1.

# temperature conditional filter on S1a
            if tt < 4.:
                tmult = 0.
            elif tt> 32.:
                tmult = 0.
            else:
                tmult = 1.

            S1t = S1t*smult*tmult

            if v2t <= 30.:
                S2t = 1.
            else:
                S2t = 1.43-0.0143*v2t

            HSI_BrShrTrawl = zero_mult*max(0,(S1t*S2t))**(1./2.)

            writestring = '%s,%s,%s,%s,%s\n' %(gridID,HSI_BrShrTrawl,st,tt,v2t)
            fBST.write(writestring)

# map juvenile shrimp HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
# map adult shrimp HSI to Ascii grid
    HSIascii_grid(HSIcsv2,HSIasc2,ascii_grid_lookup,n500cols,n500rows,ascii_header)

# delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(sal_AprJun_ave,tmp_AprJun_ave,tmp_AprJul_ave)

# delete temporary variables so they do not accidentally get used in other HSIs
    del(ss,ts,dayvs,v2s,lnCPUE1s,S1s,S2s,st,tt,dayvt,v2t,lnCPUE1t,S1t,S2t)


########################################
##        Eastern Oyster HSI          ##
########################################
    print ' Reading in cultch map for Oyster HSI'

 # 2023 Update - SES 6/18/20 commented out cultch inputs, setting cultch and just left S1=1.0 for the HSI equation.  No longer using cultch variable in HSI.   
 #   cultchdict = {}
 #   cnp = np.genfromtxt('OysterCultch.csv',skip_header=True,usecols=(0,5),delimiter=',')
 #   for row in cnp:
 #       gid = row[0]
 #       cultchdict[gid] = row[1]
 
 # 2023 Update - SES 6/18/20 added sedimentation dictionary (Eric will build) in same format as cultch input file
 # Note Eric will need to build us sediments.csv file from ICM outputs
    seddict = {}
    cnp = np.genfromtxt('sediments.csv',skip_header=True,usecols=(0,5),delimiter=',')
    for row in cnp:
        gid = row[0]
        seddict[gid] = row[1]

    print ' Calculating Eastern Oyster HSI'
# save input values that are only used in this HSI
#    sal_JanDec_min = dict((sal[n-1],min(saldict[n][jan],saldict[n][feb],saldict[n][mar],saldict[n][apr],saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb],saldict[n][nov],saldict[n][dec])) for n in range(1,n500grid+1))
# 2023 Update -- SES 6/18/20 minimum monthly salinity changed to minimum_AprSept and min_OctMar
    sal_AprSep_min = dict((sal[n-1],min(saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep])) for n in range(1,n500grid+1))
    sal_OctMar_min = dict((sal[n-1],min(saldict[n][jan],saldict[n][feb],saldict[n][mar],saldict[n][octb],saldict[n][nov],saldict[n][dec])) for n in range(1,n500grid+1))
#    sal_MaySep_ave = dict((sal[n-1],np.mean([saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep]]))for n in range(1,n500grid+1))
# 2023 Update -- SES 6/18/20 expanded mean spawning salinity from April-November
    sal_AprNov_ave = dict((sal[n-1],np.mean([saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb],saldict[n][nov]]))for n in range(1,n500grid+1))
   
    HSIcsv = r'%sOYSTE.csv' % csv_outprefix
    HSIasc = r'%sOYSTE.asc' % asc_outprefix


    with open(HSIcsv,'w') as fEO:
        
#        headerstring = 'GridID,HSI,s_maysep,s_min,s_mean,cultch,pct_land\n'
        headerstring = 'GridID,HSI,s_aprnov,s_minw,s_minc,s_mean,pct_land,sedim\n'
        fEO.write(headerstring)

# 2023 Update --- SES 6/18/20 sav2_w and _c added, sedim added and cultch turned off
        for gridID in gridIDs:
            sav1 = sal_AprNov_ave[gridID]
            sav2_w = sal_AprSep_min[gridID]
            sav2_c = sal_OctMar_min[gridID]
            sav3 = sal_JanDec_ave[gridID]
 #           cultch = cultchdict[gridID]  turned off for 2023 Update
            pland = landdict[gridID]/100.0
            sedim = seddict[gridID]
            
            
 #           if cultch <= 10.:
 #               S1 = 0.04*cultch
 #           elif cultch <= 30.:
 #               S1 = 0.02*cultch + 0.2
 #           elif cultch <= 50.:
 #               S1 = 0.01*cultch + 0.5
 #           else:
                S1 = 1.0

#            if cultch < 10.:
#                S1 = 0.
#            elif cultch < 20.:
#                S1 = 0.4
#            elif cultch < 30.:
#                S1 = 0.6
#            elif cultch < 40.:
#                S1 = 0.8
#            elif cultch < 50.:
#                S1 = 0.9
#            else:
#                S1 = 1.0



 #AH - 2017 i revised all of these below based on email from Mandy with attachment
            if sav1 < 5.:
                S2 = 0.
            elif sav1 < 10.:
                S2 = 0.06*sav1-0.3
            elif sav1 < 15.:
                S2 = 0.07*sav1-0.4
            elif sav1 < 18.:
                S2 = 0.1167*sav1-1.1
            elif sav1 < 22.:
                S2 = 1.0
            elif sav1 < 30.:
                S2 = -0.0875*sav1+2.925
            elif sav1 < 35.:
                S2 = -0.04*sav1+1.5
            elif sav1 < 40.:
                S2 = -0.02*sav1+0.8
            else:
                S2 = 0.0
# 2023 Update -- SES 6/18/20 revised for warm and cold in a combined S3
            if sav2_w < 2.:
                S3w = 0.0
            elif sav2_w < 8.:
                S3w = 0.1668*sav2_w-0.33
            elif sav2_w < 10.:
                S3w = 1.0
            elif sav2_w < 15.:
                S3w = -0.16*sav2_w+2.6
            elif sav2_w < 20.:
                S3w = -0.04*sav2_w+0.8
            else:
                S3w = 0.0

            if sav2_c < 1.:
                S3c = 0.0
            elif sav2_c < 8.:
                S3c = 0.1429*sav2_c-0.1429
            elif sav2_c < 10.:
                S3c = 1.0
            elif sav2_c < 15.:
                S3c = -0.16*sav2_c+2.6
            elif sav2_c < 20.:
                S3c = -0.04*sav2_c+0.8
            else:
                S3c = 0.0

            S3 = (S3w*S3c)**(1./2.)
# 2023 Update -- SES 6/18/20 mean annual salinity function revised
            if sav3 < 5.:
                S4 = 0.
            elif sav3 < 10.:
                S4 = 0.2*sav3-1.0
            elif sav3 < 15.:
                S4 = 1.
            elif sav3 < 20.:
                S4 = -0.16*sav3+3.4
            elif sav3 < 25.:
                S4 = -0.04*sav3+1.0
            else:
                S4 = 0.0

            S5 = 1.0-pland

            
            if sedim < 35.:
                S6 = 1.0
            elif sedim < 40.:
                S6 = -0.2*sedim+8.0
            else:
                S6 = 0.0

 #           HSI_EOys = (S1*S2*S3*S4*S5)**(1./5.)
            HSI_EOys = (S1*S2*S3*S4*S5*S6)**(1./6.)
            
 #           writestring = '%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_EOys,sav1,sav2,sav3,cultch,pland)
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_EOys,sav1,sav2_w,sav2_c,sav3,pland,sedim)
 
            fEO.write(writestring)

    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

# delete any dictionaries that aren't used in any other HSIs - frees up memory
#    del(sal_JanDec_min,sal_MaySep_ave,cultchdict)
    del(sal_AprSep_min,sal_OctMar_min,sal_AprNov_ave,seddict)

# delete temporary variables so they do not accidentally get used in other HSIs
# 2023 Update -- SES 6/18/20 removed cultch and added sav2_w, _c, sedim, and S6 variables
    del(sav1,sav2_w,sav2_c,sav3,sedim,S1,S2,S3,S4,S5,S6,cnp)


########################################
##       Spotted Seatrout HSIs        ##
########################################
##    Juvenile Spotted Seatrout HSI   ##
########################################
    print ' Calculating Juvenile and Adult Spotted Seatrout HSIs'
# save input values that are only used in this HSI
    sal_SepNov_ave = dict((sal[n-1],np.mean([saldict[n][sep],saldict[n][octb],saldict[n][nov]]))for n in range(1,n500grid+1))
    tmp_SepNov_ave = dict((tmp[n-1],np.mean([tmpdict[n][sep],tmpdict[n][octb],tmpdict[n][nov]]))for n in range(1,n500grid+1))

    HSIcsv = r'%sSPSTJ.csv' % csv_outprefix
    HSIasc = r'%sSPSTJ.asc' % asc_outprefix

    HSIcsv2 = r'%sSPSTA.csv' % csv_outprefix
    HSIasc2 = r'%sSPSTA.asc' % asc_outprefix

    
    with open(HSIcsv,'w') as fSSJ, open(HSIcsv2,'w') as fSSA:
        
        headerstring = 'GridID,HSI,s,t,v2\n'
        fSSJ.write(headerstring)
        fSSA.write(headerstring)

        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            sj = sal_SepNov_ave[gridID]
            tj = tmp_SepNov_ave[gridID]
            v2j =  max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))
            dayvj = 2.9006

# Calculate variables for Juvenile Spotted Seatrout
            lnCPUE1j = -8.6532+6.2748*dayvj-1.1591*dayvj**2.+0.0251*sj+0.07216*tj-0.00077*sj**2.-0.00000085*sj**2.*tj**2.-0.00168*tj**2.

            S1j = (e**lnCPUE1j - 1.)/1.029

            if v2j < 25.:
                S2j = 0.5 + 0.02*v2j
            elif v2j <= 80:
                S2j = 1.0
            else:
                S2j = 5.0 - 0.05*v2j
            
            HSI_juvSStr = zero_mult*(max(0,S1j)**(2./3.))*(max(0,S2j)**(1./3.))
            
            writestring = '%s,%s,%s,%s,%s\n' %(gridID,HSI_juvSStr,sj,tj,v2j)
            fSSJ.write(writestring)
            

########################################
##      Adult Spotted Seatrout HSI    ##
########################################
            sa = sal_JanDec_ave[gridID]
            ta = tmp_JanDec_ave[gridID]
            v2a = v2j
            dayva = 1.805

# Calculate variables for Adult Spotted Seatrout
            lnCPUE1a = -0.2433-0.00983*dayva-0.0109*dayva**2.-0.02731*sa+0.0904*ta+0.00357*sa*ta+0.00144*sa**2.+0.000007*sa**2.*ta**2.-0.00027*ta*sa**2.-0.00233*ta**2.

            S1a = (e**lnCPUE1a - 1.)/2.869

            if v2a < 25.:
                S2a = 0.7 +0.012*v2a
            elif v2a <= 70.:
                S2a = 1.0
            elif v2a < 100.:
                S2a = 3.33 - 0.033*v2a
            else:
                S2a = 0.0
                
            HSI_adltSStr = zero_mult*max(0,(S1a*S2a))**(1./2.)
            
            writestring = '%s,%s,%s,%s,%s\n' %(gridID,HSI_adltSStr,sa,ta,v2a)
            fSSA.write(writestring)

# map juvenile seatrout HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
# map adult seatrout HSI to Ascii grid
    HSIascii_grid(HSIcsv2,HSIasc2,ascii_grid_lookup,n500cols,n500rows,ascii_header)

# delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(sal_SepNov_ave,tmp_SepNov_ave)

# delete temporary variables so they do not accidentally get used in other HSIs
    del(sj,tj,v2j,dayvj,lnCPUE1j,S1j,S2j,sa,ta,v2a,dayva,lnCPUE1a,S1a,S2a)


########################################
##          White Shrimp HSIs         ##
########################################

    print ' Calculating White Shrimp HSIs'
# save input values that are only used in this HSI
    sal_JunNov_ave = dict((sal[n-1],np.mean([saldict[n][jun],saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb],saldict[n][nov]]))for n in range(1,n500grid+1))
    tmp_JunNov_ave = dict((tmp[n-1],np.mean([tmpdict[n][jun],tmpdict[n][jul],tmpdict[n][aug],tmpdict[n][sep],tmpdict[n][octb],tmpdict[n][nov]]))for n in range(1,n500grid+1))
    sal_JulOct_ave = dict((sal[n-1],np.mean([saldict[n][jul],saldict[n][aug],saldict[n][sep],saldict[n][octb]]))for n in range(1,n500grid+1))
    tmp_JulOct_ave = dict((tmp[n-1],np.mean([tmpdict[n][jul],tmpdict[n][aug],tmpdict[n][sep],tmpdict[n][octb]]))for n in range(1,n500grid+1))
    
    HSIcsv = r'%sWSHRS.csv' % csv_outprefix
    HSIasc = r'%sWSHRS.asc' % asc_outprefix

    HSIcsv2 = r'%sWSHRL.csv' % csv_outprefix
    HSIasc2 = r'%sWSHRL.asc' % asc_outprefix


    ##do not calculate Large Juvenile White Shrimp Trawl HSI            
    #with open(HSIcsv,'w') as fWSS, open(HSIcsv2,'w') as fWST:
    ##do not calculate Large Juvenile White Shrimp Trawl HSI            
    
    with open(HSIcsv,'w') as fWSS:
        headerstring = 'GridID,HSI,s,t,v2\n'
        fWSS.write(headerstring)
 
 ##do not calculate Large Juvenile White Shrimp Trawl HSI            
 #       fWST.write(headerstring)
 ##do not calculate Large Juvenile White Shrimp Trawl HSI                   
 
        for gridID in gridIDs:
            zero_mult = land_mult[gridID]*fresh_for_mult[gridID]*bare_mult[gridID]
            ss = sal_JunNov_ave[gridID]
            ts = tmp_JunNov_ave[gridID]
            v2s =  max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))
            dayvs = 2.543

################################################
##   Small Juvenile White Shrimp HSI - Seine  ##
################################################

            lnCPUE1s = -7.9150+8.1556*dayvs-1.6101*dayvs**2.-0.03471*ss+0.009624*ts+0.01385*ss*ts-0.00264*ss**2.-0.00115*ts**2.-0.00034*ss*ts**2.

            S1s = (e**lnCPUE1s - 1.)/24.68

            if v2s < 25.:
                S2s = 0.028*v2s+0.3
            elif v2s <= 80.:
                S2s = 1.
            else:
                S2s = 5.-0.05*v2s

            HSI_juvWhShrS = zero_mult*max(0,(S1s*S2s))**(1./2.)
            
            writestring = '%s,%s,%s,%s,%s\n'  %(gridID,HSI_juvWhShrS,ss,ts,v2s)
            fWSS.write(writestring)

##do not calculate Large Juvenile White Shrimp Trawl HSI            
##
##################################################
####   Large Juvenile White Shrimp HSI - Trawl  ##
##################################################
##
##
##            st = sal_JulOct_ave[gridID]
##            tt = tmp_JulOct_ave[gridID]
##            v2t = v2s
##            dayvt = 2.394
##
##            lnCPUE1t = -08.2435+1.248*st+0.3732*tt-0.03476*st**2.-0.00751*tt**2.-0.08411*st*tt-0.00004*st**2.*tt**2.+0.00141*st*tt**2.+0.002221*tt*st**2.+4.6448*dayvt-0.8693*dayvt**2.
##
##            S1t = (e**lnCPUE1t - 1.)/53.93
##
##            if v2t > 30.:
##                S2t = -0.0143*v2t+1.429
##            else:
##                S2t = 1.
##
##            HSI_juvWhShrT = zero_mult*max(0,(S1t*S2t))**(1./2.)
##            
##            writestring = '%s,%s,%s,%s,%s\n'  %(gridID,HSI_juvWhShrT,st,tt,v2t)
##            fWST.write(writestring)
##do not calculate Large Juvenile White Shrimp Trawl HSI            

# map white shrimp seine HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

##do not calculate Large Juvenile White Shrimp Trawl HSI            
### map white shrimp trawl HSI to Ascii grid
##    HSIascii_grid(HSIcsv2,HSIasc2,ascii_grid_lookup,n500cols,n500rows,ascii_header)
##do not calculate Large Juvenile White Shrimp Trawl HSI            

# delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(sal_JunNov_ave,tmp_JunNov_ave,sal_JulOct_ave,tmp_JulOct_ave)

# delete temporary variables so they do not accidentally get used in other HSIs
    del(ss,ts,v2s,dayvs,lnCPUE1s,S1s,S2s)
    
    ##do not calculate Large Juvenile White Shrimp Trawl HSI            
    #del(st,tt,v2t,dayvt,lnCPUE1t,S1t,S2t)
    ##do not calculate Large Juvenile White Shrimp Trawl HSI            


######################################
##         MOTTLED DUCK HSI         ##
######################################
# This HSI assumes is coded based on the depth tables generated by WM.HSIreclass provide depth in centimeters
    print ' Converting Mottled Duck Depth DBF to CSV file.'
    MotDuckDepDBF = os.path.normpath("%s\\MotDuckDepths_cm.dbf" % HSI_dir) # must match name of file set in WM.HSIreclass
    MotDuckCSV = os.path.normpath("%s\\MotDuckDepths_cm_%s.csv" % (HSI_dir,year))
    
    with open(MotDuckCSV,'w') as csvfile:
        idb = dbf.Dbf(MotDuckDepDBF)
        ocsv = csv.writer(csvfile)
        names =[]
        for field in idb.header.fields:
            names.append(field.name)
        ocsv.writerow(names)
        for rec in idb:
            ocsv.writerow(rec.fieldData)
        idb.close()
    
    print ' - removing DBF file.'
    os.remove(MotDuckDepDBF)
    
    print ' Calculating Mottled Duck HSI'
    MotDuck = np.genfromtxt(MotDuckCSV,delimiter = ',',skiprows = 1)
    MotDuckDepdict = {}
    MDDdict ={}
    MDDmissing = [0,0,0,0,0,0,0,0,0]        # values to assign to grid cells with no GridID key in depth dictionary
    
    # convert MotDuck depths array into dictionary, GridID is key (only grid cells overlaying geospatial extent in WM.CalculateEcohydroAttributes() will have a key and values
    for n in range(0,len(MotDuck)):
        gridIDinD = MotDuck[n,0]
        MDDdict[gridIDinD]=MotDuck[n,1:10]
    
    # generate dictionary of various depths for all gridID values        
    for gridID in gridIDs:
        try:
            MotDuckDepdict[gridID] = MDDdict[gridID]
        except:
            MotDuckDepdict[gridID] = MDDmissing

    HSIcsv = r'%sMOTDK.csv' % csv_outprefix
    HSIasc = r'%sMOTDK.asc' % asc_outprefix

    with open(HSIcsv,'w') as fMD:
        
        headerstring = 'GridID,HSI,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v4\n'
        fMD.write(headerstring)

        for gridID in gridIDs:
            s = sal_JanDec_ave[gridID]
            w = waterdict[gridID]
            warea = max(wetlndict[gridID],0.001)       #total wetland area of 500 m cell
            
            v1a = frattdict[gridID]
            v1b = frfltdict[gridID]
            v1c = interdict[gridID]
            v1d = brackdict[gridID]
            v1e = salmardict[gridID]
            v1f = swfordict[gridID]
            v1g = btfordict[gridID]
            
            
            vegland = v1a + v1b + v1c + v1d + v1e + v1f  # percent land as summarized by veg output, exclude v1g b/c it is the same as v1f
            
            # Reclassify water area as marsh type based on salinity value
            # initialize additional marsh areas to zero
            w_fresh = 0.0
            w_inter = 0.0
            w_brack = 0.0
            w_saline = 0.0
            
            # classify water areas based on salinity to add to wetland areas in S1 equation
            # vegetation land may not exactly match percent water due to differences in morph and veg output - therefore check that percent water + percent land from veg output is not greater than 100.0
            if s < 1.5:
                w_fresh = max(0,min(w,100.0-vegland))
            elif s < 4.5:
                w_inter = max(0,min(w,100.0-vegland))
            elif s < 9.5:
                w_brack = max(0,min(w,100.0-vegland))
            else:
                w_saline = max(0,min(w,100.0-vegland))
            
            S1 = (v1a + w_fresh)/100.0 + v1b/100.0 + 0.67*(v1c + w_inter)/100.0 + 0.55*(v1d +  w_brack)/100.0 + 0.23*(v1e + w_saline)/100.0
 
            v2 = min(1.0,max(0.0,wetlndict[gridID]/100.0))      #wetlndict is in percentages, therefore divide by 100 to get V2 in correct unit
            if v2 < 0.32:
                S2 = 0.1 + 2.81*v2
            elif v2 <= 0.70:
                S2 = 1.0
            else:
                S2 = 3.1 - 3.0*v2

            area = 0.0
            for x in range(0,9):    # x is number of columns in depth dictionary (as summarized in WM.HSIreclass)
                area = area + MotDuckDepdict[gridID][x] #determine area of cell analyzed when developing depth values in morph(not exactly equal to 500x500 since the 30x30 m grid doesn't fit in the 500x500
            if area < 250000:
                area = 500*500
            less0 = MotDuckDepdict[gridID][0]/area   # portion of cell less than 0-cm deep    
            v3a = MotDuckDepdict[gridID][1]/area    # portion of cell 0-8 cm deep
            v3b = MotDuckDepdict[gridID][2]/area    # portion of cell 8-30 cm deep 
            v3c = MotDuckDepdict[gridID][3]/area    # portion of cell 30-36 cm deep
            v3d = MotDuckDepdict[gridID][4]/area    # portion of cell 36-42 cm deep
            v3e = MotDuckDepdict[gridID][5]/area    # portion of cell 42-46 cm deep
            v3f = MotDuckDepdict[gridID][6]/area    # portion of cell 46-50 cm deep
            v3g = MotDuckDepdict[gridID][7]/area    # portion of cell 50-56 cm deep
            v3h = MotDuckDepdict[gridID][8]/area    # portion of cell more than 56 cm deep

            S3 = 0.6*v3a+1.0*v3b+0.83*v3c+0.57*v3d+0.35*v3e+0.22*v3f+0.009*v3g

            v4 = sal_AprJul_ave[gridID]

            if v4 <= 9.0:
                S4 = 1.0
            elif v4 < 18:
                S4 = 2.0 - 0.11*v4
            else:
                S4 = 0.0

            HSI_MotDuck = (S1*S2*S3*S4)**(1./4.)
            
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (gridID,HSI_MotDuck,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v4)
            fMD.write(writestring)

# map mottled duck HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    # delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(MotDuckDepdict,MotDuck,MDDdict)
    
    # delete temporary variables so they do not accidentally get used in other HSIs
    del(area,warea,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v4,S1,S2,S3,S4)

######################################
##           GADWALL HSI            ##
######################################
# This HSI assumes is coded based on the depth tables generated by WM.HSIreclass provide depth in centimeters
    print ' Converting Gadwall Depth DBF to CSV file.'
    GadwallDepDBF = os.path.normpath("%s\\GadwallDepths_cm.dbf" % HSI_dir) # must match name of file set in WM.HSIreclass
    GadwallCSV = os.path.normpath("%s\\GadwallDepths_cm_%s.csv" % (HSI_dir,year))
    
    with open(GadwallCSV,'w') as csvfile:
        idb = dbf.Dbf(GadwallDepDBF)
        ocsv = csv.writer(csvfile)
        names =[]
        for field in idb.header.fields:
            names.append(field.name)
        ocsv.writerow(names)
        for rec in idb:
            ocsv.writerow(rec.fieldData)
        idb.close()
    
    print ' - removing DBF file.'
    os.remove(GadwallDepDBF)
    
    print ' Calculating Gadwall HSI'
    Gadwall = np.genfromtxt(GadwallCSV,delimiter = ',',skiprows = 1)
    GadwallDepdict = {}
    GDdict ={}
    GDmissing = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]        # values to assign to grid cells with no GridID key in depth dictionary
    
    # convert Gadwall depths array into dictionary, GridID is key (only grid cells overlaying geospatial extent in WM.CalculateEcohydroAttributes() will have a key and values
    for n in range(0,len(Gadwall)):
        gridIDinD = Gadwall[n,0]
        GDdict[gridIDinD]=Gadwall[n,1:15]
    
    # generate dictionary of various depths for all gridID values        
    for gridID in gridIDs:
        try:
            GadwallDepdict[gridID] = GDdict[gridID]
        except:
            GadwallDepdict[gridID] = GDmissing
    # initialize deepwat dictionary to be used to save deepwater from Gadwall Depths for use in alligator HSI
    deepwat={}

    HSIcsv = r'%sGADWA.csv' % csv_outprefix
    HSIasc = r'%sGADWA.asc' % asc_outprefix

    with open(HSIcsv,'w') as fG:
        
        headerstring = 'GridID,HSI,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v3i,v3j,v3k,v3l\n'
        fG.write(headerstring)
      
        for gridID in gridIDs:
            s = sal_JanDec_ave[gridID]
            w = waterdict[gridID]
            v1a = frattdict[gridID]
            v1b = frfltdict[gridID]
            v1c = interdict[gridID]
            v1d = brackdict[gridID]
            v1e = salmardict[gridID]
            v1f = swfordict[gridID]
            v1g = btfordict[gridID]

            vegland = v1a + v1b + v1c + v1d + v1e + v1f # percent land as summarized by veg output, exclude v1g b/c it is the same as v1f
            
            # Reclassify water area as marsh type based on salinity value
            # initialize additional marsh areas to zero
            w_fresh = 0.0
            w_inter = 0.0
            w_brack = 0.0
            w_saline = 0.0
            
            # classify water areas based on salinity to add to wetland areas in S1 equation
            # vegetation land may not exactly match percent water due to differences in morph and veg output - therefore check that percent water + percent land from veg output is not greater than 100.0
            if s < 1.5:
                w_fresh = max(0,min(w,100.0-vegland))
            elif s < 4.5:
                w_inter = max(0,min(w,100.0-vegland))
            elif s < 9.5:
                w_brack = max(0,min(w,100.0-vegland))
            else:
                w_saline = max(0,min(w,100.0-vegland))
            
            #S1 = 0.68*(v1a + w_fresh)/100.0 + 0.68*v1b/100.0 + (v1c + w_inter)/100.0 + 0.50*(v1d +  w_brack)/100.0 + 0.09*(v1e + w_saline)/100.0 + 0.25*v1f/100.0 + 0.25*v1g/100.0
            # 2023 Update -- SES 6/18/20
            S1 = 0.68*(v1a + w_fresh)/100.0 + 0.68*v1b/100.0 + (v1c + w_inter)/100.0 + 0.50*(v1d +  w_brack)/100.0 + 0.09*(v1e + w_saline)/100.0 + 0.05*v1f/100.0 + 0.05*v1g/100.0
            
            v2 = watsavdict[gridID]/100.0     #watsavdict is in percentages, therefore divide by 100 to get V2 in correct unit
            if v2 < 0.30:
                S2 = 0.08
            elif v2 < 0.70:
                S2 = 2.3*v2 - 0.61
            else:  
                S2 = 1.0
                   
            area = 0.0
            for x in range(1,14):
                area = area + GadwallDepdict[gridID][x]  #determine area of cell (not exactly equal to 500x500 since the 30x30 m grid doesn't fit in the 500x500
            if area < 250000:
                area = 500*500
            less0 = GadwallDepdict[gridID][0]/area       # portion of cell less than 0 cm deep
            v3a = GadwallDepdict[gridID][1]/area         # portion of cell 0-4 cm deep
            v3b = GadwallDepdict[gridID][2]/area         # portion of cell 4-8 cm deep
            v3c = GadwallDepdict[gridID][3]/area         # portion of cell 8-12 cm deep
            v3d = GadwallDepdict[gridID][4]/area         # portion of cell 12-18 cm deep
            v3e = GadwallDepdict[gridID][5]/area         # portion of cell 18-22 cm deep
            v3f = GadwallDepdict[gridID][6]/area         # portion of cell 22-28 cm deep
            v3g = GadwallDepdict[gridID][7]/area         # portion of cell 28-32 cm deep
            v3h = GadwallDepdict[gridID][8]/area         # portion of cell 32-36 cm deep
            v3i = GadwallDepdict[gridID][9]/area         # portion of cell 36-40 cm deep
            v3j = GadwallDepdict[gridID][10]/area        # portion of cell 40-44 cm deep
            v3k = GadwallDepdict[gridID][11]/area        # portion of cell 44-78 cm deep
            v3l = GadwallDepdict[gridID][12]/area        # portion of cell 78-150 cm deep
                   
       # save deep water from GadwallDepdict for use in Alligator HSI
            deepwat[gridID] = GadwallDepdict[gridID][13]/area        # portion of cell greater than 150 cm deep
                   
                   
            S3 = 0.05*v3a + 0.15*v3b + 0.35*v3c + 0.6*v3d + 0.83*v3e + 1.0*v3f + 0.86*v3g + 0.61*v3h + 0.37*v3i + 0.2*v3j + 0.1*v3k + 0.05*v3l
                   
            HSI_Gadwall = (S1*S2*S3)**(1./3.)
            
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (gridID,HSI_Gadwall,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v3i,v3j,v3k,v3l)
            fG.write(writestring)

# map gadwall HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

# delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(GadwallDepdict,Gadwall,GDdict)
    
# delete temporary variables so they do not accidentally get used in other HSIs
    del(area,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v3i,v3j,v3k,v3l,S1,S2,S3)

######################################
##       GREEN-WINGED TEAL HSI      ##
######################################
# This HSI assumes is coded based on the depth tables generated by WM.HSIreclass provide depth in centimeters
    print ' Converting Green-winged Teal depth DBF to CSV file.'
    GWTealDepDBF = os.path.normpath("%s\\GWTealDepths_cm.dbf" % HSI_dir) # must match name of file set in WM.HSIreclass
    GWTealCSV = os.path.normpath("%s\\GWTealDepths_cm_%s.csv" % (HSI_dir,year)) 
    
    with open(GWTealCSV,'w') as csvfile:
        idb = dbf.Dbf(GWTealDepDBF)
        ocsv = csv.writer(csvfile)
        names =[]
        for field in idb.header.fields:
            names.append(field.name)
        ocsv.writerow(names)
        for rec in idb:
            ocsv.writerow(rec.fieldData)
        idb.close()
    
    print ' - removing DBF file.'
    os.remove(GWTealDepDBF)
    
    print ' Calculating Green-winged Teal HSI'
    GWTeal = np.genfromtxt(GWTealCSV,delimiter = ',',skiprows = 1)
    GWTealDepdict = {}
    GTDdict ={}
    GTDmissing = [0,0,0,0,0,0,0,0,0]        # values to assign to grid cells with no GridID key in depth dictionary
    
    # convert Gadwall depths array into dictionary, GridID is key (only grid cells overlaying geospatial extent in WM.CalculateEcohydroAttributes() will have a key and values
    for n in range(0,len(GWTeal)):
        gridIDinD = GWTeal[n,0]
        GTDdict[gridIDinD]=GWTeal[n,1:10]
    
    # generate dictionary of various depths for all gridID values        
    for gridID in gridIDs:
        try:
            GWTealDepdict[gridID] = GTDdict[gridID]
        except:
            GWTealDepdict[gridID] = GTDmissing

    HSIcsv = r'%sGTEAL.csv' % csv_outprefix
    HSIasc = r'%sGTEAL.asc' % asc_outprefix

    with open(HSIcsv,'w') as fGT:
        
        headerstring = 'GridID,HSI,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g\n'
        fGT.write(headerstring)

        for gridID in gridIDs:
            warea = max(wetlndict[gridID],0.001)       #total wetland area of 500 m cell
            s = sal_JanDec_ave[gridID]
            w = waterdict[gridID]
            v1a = frattdict[gridID]
            v1b = frfltdict[gridID]
            v1c = interdict[gridID]
            v1d = brackdict[gridID]
            v1e = swfordict[gridID]
            v1f = btfordict[gridID]
            v1g = salmardict[gridID]

            vegland = v1a + v1b + v1c + v1d + v1f + v1g # percent land as summarized by veg output, exclude v1e b/c it is the same as v1f
            
            # Reclassify water area as marsh type based on salinity value
            # initialize additional marsh areas to zero
            w_fresh = 0.0
            w_inter = 0.0
            w_brack = 0.0
            w_saline = 0.0
            
            # classify water areas based on salinity to add to wetland areas in S1 equation
            # vegetation land may not exactly match percent water due to differences in morph and veg output - therefore check that percent water + percent land from veg output is not greater than 100.0
            if s < 1.5:
                w_fresh = max(0,min(w,100.0-vegland))
            elif s < 4.5:
                w_inter = max(0,min(w,100.0-vegland))
            elif s < 9.5:
                w_brack = max(0,min(w,100.0-vegland))
            else:
                w_saline = max(0,min(w,100.0-vegland))
            
            S1 = (v1a + w_fresh)/100.0 + v1b/100.0 + 0.60*(v1c + w_inter)/100.0 + 0.93*(v1d +  w_brack)/100.0 + 0.46*(v1e + w_saline)/100.0 + 0.25*v1f/100.0 + 0.25*v1g/100.0

            v2 = max(0.0,min(watsavdict[gridID]+waterdict[gridID],100.0))/100.0     #watsavdict is in percentages,  therefore divide by 100 to get V2 in correct unit
            
            if v2 < 0.35:
                S2 = 0.1 + 2.5*v2
            elif v2 <= 0.75:
                S2 = 1.0
            else:
                S2 = 3.7 - 3.6*v2

            area = 0.0
            for x in range(0,8):
                area = area + GWTealDepdict[gridID][x]   #determine area of cell (not exactly equal to 500x500 since the 30x30 m grid doesn't fit in the 500x500
            if area < 250000:
                area = 500*500
            less0 = GWTealDepdict[gridID][0]/area         # portion of cell less than 0 cm deep
            v3a = GWTealDepdict[gridID][1]/area           # portion of cell 0-6 cm deep
            v3b = GWTealDepdict[gridID][2]/area           # portion of cell 6-18 cm deep
            v3c = GWTealDepdict[gridID][3]/area           # portion of cell 18-22 cm deep
            v3d = GWTealDepdict[gridID][4]/area           # portion of cell 18-26 cm deep
            v3e = GWTealDepdict[gridID][5]/area           # portion of cell 26-30 cm deep
            v3f = GWTealDepdict[gridID][6]/area           # portion of cell 30-34 cm deep
            v3g = GWTealDepdict[gridID][7]/area           # portion of cell 34-100 cm deep
            v3h = GWTealDepdict[gridID][8]/area           # portion of cell more than 100 cm deep

            S3 = 0.8*v3a + 1.0*v3b + 0.87*v3c + 0.68*v3d + 0.43*v3e + 0.17*v3f + 0.07*v3g

            HSI_GWTeal = (S1*S2*S3)**(1./3.)

            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (gridID,HSI_GWTeal,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g)
            fGT.write(writestring)

# map green winged teal HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

# delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(GWTealDepdict,GTDdict,GWTeal)

# delete temporary variables so they do not accidentally get used in other HSIs
    del(area,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,S1,S2,S3)

######################################
##         BROWN PELICAN HSI        ##
######################################
    print ' Calculating Brown Pelican HSI'


    HSIcsv = r'%sBRWNP_noGMENA.csv' % csv_outprefix
    HSIasc = r'%sBRWNP_noGMENA.asc' % asc_outprefix
    
    HSIcsv2 = r'%sBRWNP.csv' % csv_outprefix
    HSIasc2 = r'%sBRWNP.asc' % asc_outprefix
    
    bpel_input_file = 'BrownPelican_HSI_inputs_%02d.csv' % elapsedyear 
    BP_inputs = np.genfromtxt(bpel_input_file,skiprows=1,delimiter=',')
    
    distancemultiplier = {}
    saltmarsh_islandarea_m2 = {}
    
    for n in range(0,len(BP_inputs)):
        gridIDinBP = BP_inputs[n][0]
        saltmarsh_islandarea_m2[gridIDinBP] = BP_inputs[n][1]
        distancemultiplier[gridIDinBP] = BP_inputs[n][2]
    
    del BP_inputs
           
    with open(HSIcsv,'w') as fBP:
        
        headerstring = 'GridID,HSI,area_ha,distance_multiplier,v3,s4,s5=1,s6\n'
        fBP.write(headerstring)

        for gridID in gridIDs:
            try:
                area_ha = saltmarsh_islandarea_m2[gridID]/10000.0
            except:
                area_ha = 0.0
            
            # if there is no area value, there is no small island with salt marsh present in grid cell - set HSI to zero
            #if area_ha <= 0.0: # update this lower value to 25 ha  - if island is smaller than one 500mx500m grid, HSI is zero # fix for alternative runs
            if area_ha <= 25.0:
                HSI_BP = 0.0
                S2 = -9999
                v3 = -9999
                S4 = -9999
                S5 = -9999
                S6 = -9999
            # if the island area is larger than 200 ha, island is too large - set HSI to zero
            elif area_ha > 200:
                HSI_BP = 0.0
                S2 = -9999
                v3 = -9999
                S4 = -9999
                S5 = -9999
                S6 = -9999
            # otherwise first S term is function of island size
            else:
                if area_ha <= 180:
                    S1 = 1.0
                elif area_ha <= 200:
                    S1 = 10 - 0.05*area_ha
                    
                # if there is a distance multiplier value for the grid cell, use, otherwise set to zero so no HSI is calculated
                # if small island is within 1.0 km of land, multiplier = 0
                # if small island is within 1.5 km of land, multiplier = 0.2
                # if small island is within 2.0 km of land, multiplier = 0.4
                # if small island is within 2.5 km of land, multiplier = 0.6
                # if small island is within 3.0 km of land, multiplier = 0.8
                # if small island is more than 3.0 km from land, multiplier = 1.0
                try:
                    S2 = distancemultiplier[gridID]
                except:
                    S2 = 0.0
                
                try: 
                    v3bm = blackmangrovedict[gridID]
                except:
                    v3bm = 0.0
                try:
                    v3me = marshelderdict[gridID]
                except:
                    v3me = 0.0
                    
                v3 = min(max(0,v3bm + v3me),1.0)
                
                if v3 >= 0.5:
                    S3 = 1.0
                else:
                    S3 = 1.6*v3 + 0.2
                    
                # S4 is distance to human activity - this data is not in the model, but the 500-m grid structure excludes developed areas, therefore it is not included in this HSI code 
                S4 = 1.0
                
                # S5 is menhaden HSI value withing 20km radius  - this is added after mapped to ASCII grid
                S5 = 1.0                
                # S6 is the dominant emergent vegetation type - this term is zero for all types EXCEPT salt marsh
                # when island sizes are calculated in WM.HSIpelican() function, the output areas are only calculated for islands that have some salt marsh
                # therefore this term is set to 1 in this equation since the area term will be zero for non-salt marsh areas
                S6 = 1.0
                
                HSI_BP = (S1*S2*S3*S4*S5*S6)**(1./6.)
                
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_BP,area_ha,S2,v3,S4,S5,S6) 
            fBP.write(writestring)

# map pelican HSI to Ascii grid without Menhaden habitat values
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
    
    print ' Re-calculating Brown Pelican HSI with Menhaden'
    
    with open(HSIcsv2,'w') as fBP2:
        
        headerstring = 'GridID,HSI,area_ha,distance_multiplier,v3,s4,s5,s6\n'
        fBP2.write(headerstring)

        menhaden_asci = r'%sGMENA.asc' % asc_outprefix
        
        men = np.genfromtxt(menhaden_asci,delimiter=" ",skiprows=6)
        pelcsv = np.genfromtxt(HSIcsv,delimiter=",",skiprows=1)
        peldict = dict((pelcsv[n][0],pelcsv[n][1:7])for n in range(0,len(pelcsv)))

        newHSIgrid = np.zeros([n500rows,n500cols])
        for m in range(0,n500rows):
            for n in range(0,n500cols):
                cellID = ascii_grid_lookup[m][n]
                if cellID == -9999:
                    newHSIgrid[m][n] = -9999
                else:
                    pelHSInoMen = peldict[cellID][0] 
                    area_ha = peldict[cellID][1]
                    S2 = peldict[cellID][2]
                    v3 = peldict[cellID][3]
                    S4 = peldict[cellID][4]
                    S6 = peldict[cellID][5]
                    
                    if pelHSInoMen == 0.0:
                        newHSIgrid[m][n] = 0.0
                    # if grid cell has data and a non-zero pelican HSI, loop over surrounding cells and determine average adult menhaden HSI
                    else:
                        menave = 0
                        avecells_n = 0
                        S5 = 0
                    # look at all cells 35 rows and 35 columns away from current grid cell
                    # this 71 grid cell wide surrounding area is equal to 1260.5 sq km, which is equivalent to the 20-km search radius that brown pelicans have, which is 1256.6 sq km
                        for mm in range(-35,36):
                            newrow = m + mm
                            for nn in range(-35,36):    
                                newcol = n + nn
                                surroundcell = ascii_grid_lookup[newrow][newcol]
                                if surroundcell <> -9999:
                                    menave += men[newrow][newcol]
                                    avecells_n += 1
                            if avecells_n <> 0:
                                V5 = (menave/avecells_n)**(1./6.)
                            else:
                                V5 = 0.0
                            if V5 >= 0.6:
                                S5 = 1.0
                            else:
                                S5 = V5*5./3.
                        newHSIgrid[m][n] = pelHSInoMen*(menave/avecells_n)**(1./6.)
                    
                    writestring = '%s,%s,%s,%s,%s,%s,%s,%s\n' %(cellID,newHSIgrid[m][n],area_ha,S2,v3,S4,S5,S6) 
                    fBP2.write(writestring)
    HSIascii_grid(HSIcsv2,HSIasc2,ascii_grid_lookup,n500cols,n500rows,ascii_header)       
                            
                  

############################
##      CRAWFISH HSI      ##
############################
#   sal_OctJun_ave = dict((sal[n-1],np.mean([saldict[n][jan],saldict[n][feb],saldict[n][mar],saldict[n][mar],saldict[n][apr],saldict[n][may],saldict[n][jun],saldict[n][octb],saldict[n][nov],saldict[n][dec]))for n in range(1,n500grid+1))
    print ('Calculating Crawfish HSI.')
    # 2023 Update - SES 6/18/20 - saving new inputs for just this HSI - this changed from 2017 so old inputs are still saved above (dep_OctJun_ave[gridID] and dep_JulSep_ave[gridID])
    dep_DecJul_ave = dict((dept[n-1],np.mean([depthdict[n][jan],depthdict[n][feb],depthdict[n][mar],depthdict[n][apr],depthdict[n][may],depthdict[n][jun],depthdict[n][jul],depthdict[n][dec]]))for n in range(1,n500grid+1))
    dep_AugNov_ave = dict((dept[n-1],np.mean([depthdict[n][aug],depthdict[n][sep],depthdict[n][octb],depthdic[n][nov]]))for n in range(1,n500grid+1))
    
    HSIcsv = r'%sCRAYF.csv' % csv_outprefix
    HSIasc = r'%sCRAYF.asc' % asc_outprefix

    with open(HSIcsv,'w') as fCF:
        
        #headerstring = 'GridID,HSI,s,dep_OctJun,dep_JulSep,swamp_for,fresh,water,inter,brack,sand\n'
        # 2023 Update - SES 6/18/20 - created new header b/c removed sand as variable and changed dep_Months
        headerstring = 'GridID,HSI,s,dep_DecJul,dep_AugNov,swamp_for,fresh,water,inter,brack\n'
        fCF.write(headerstring)
    
        for gridID in gridIDs:
            s = sal_JanDec_ave[gridID]
            #depoj = 100*dep_OctJun_ave[gridID]
            # 2023 Update -- SES 6/18/20 
            depdj = 100*dep_DecJul_ave[gridID]
            #depjs = 100*dep_JulSep_ave[gridID]
            # 2023 Update -- SES 6/18/20
            depan = 100*dep_AugNov_ave[gridID]
            swampfor = swfordict[gridID]
            fresh = frattdict[gridID]+ frfltdict[gridID]
            wat = waterdict[gridID]
            inter = interdict[gridID]
            brack = brackdict[gridID]
            #sand = pctsanddict[gridID] - 2023 Update -- SES 6/18/20 - set sand and S4 to dummy values not used and to delete at end
            sand = 999
            
            if s <= 1.5:
                S1 = 1
            elif s <= 3.0:
                S1 = 1.5 - s/3.0
            elif s <= 6.0:
                S1 = 1 - s/6.0
            else:
                S1 = 0.0
            
            if depdj <= 0:
                S2 = 0.0
            elif depdj <= 46:
                S2 = depdj/46.0
            elif depdj <= 91:
                S2 = 1.0
            elif depdj <= 274:
                S2 = 1.5 - 1.5*depdj/274.0
            else:
                S2 = 0.0    
            
            S3= swampfor/100.0 + 0.85*fresh/100.0 + 0.75*wat/100. + 0.6*inter/100.0 + 0.2*brack/100.0

            #if sand >= 90:
            #    S4 = 0.0
            #elif sand > 50:
            #    S4 = 2.25-2.25*sand/50.0
            #else:
            #    S4 = 1.0
            S4=0.0

            if depan <= 0:
                S5 = 1.0
            elif depan <= 15:
                S5 = 1.0 - depan/15.0
            else:
                S5 = 0.0

            #HSI_CF = ((S1*S2)**(1./6.))*(S3**(1./3.))*((S4*S5)**(1./6.))
            # 2023 Update -- SES 6/18/20 - took out sand and then removed variable from output file below 
            HSI_CF = ((S1*S2)**(1./6.))*(S3**(1./3.))*(S5**(1./3.))
            
            #writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_CF,s,depoj,depjs,swampfor,fresh,wat,inter,brack,sand)
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_CF,s,depdj,depan,swampfor,fresh,wat,inter,brack)
            fCF.write(writestring)

# map crawfish HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    del(s,depdj,depan,swampfor,fresh,wat,inter,brack,sand,S1,S2,S3,S4,S5,HSI_CF)
    


 

############################
##     ALLIGATOR HSI      ##
############################
    print ('Calculating Alligator HSI.')
    
    
    HSIcsv = r'%sALLIG.csv' % csv_outprefix
    HSIasc = r'%sALLIG.asc' % asc_outprefix

    with open(HSIcsv,'w') as fAl:
        
        headerstring = 'GridID,HSI,pct_wat,depth,bald,fresh,int,brack,deep,s,pct_edge\n'
        fAl.write(headerstring)
        
        for gridID in gridIDs:
            pwat = waterdict[gridID]
            
            # if there is no elevation for the marsh (e.g. no marsh in grid), hard code marsh-relative depth to 9.999 - this will default S2 to smallest value of 0.1
            if melevdict > -9990.0:
                deprelmar = stagedict[gridID] - melevdict[gridID]
            else:
                deprelmar = 9.999
            
            baldcyp = baldcypdict[gridID]
            fresh = max(0.0,min(frfltdict[gridID] + frattdict[gridID],1.0))
            inter = interdict[gridID]
            brack = brackdict[gridID]
            deep = 100.0*deepwat[gridID]  #deepwat is the portion of the cell that is deeper than 1.5 meters (calculated from Gadwall depths)
            s = sal_JanDec_ave[gridID]
            edge = pctedgedict[gridID]
            
            if pwat < 20:
                S1 = (4.5*pwat/100.)+0.1
            elif pwat < 40:
                S1 = 1
            else:
                S1 = (-1.6667*pwat/100.)+1.6667
            
            if deprelmar < -0.55:
                S2 = 0.1
            elif deprelmar < -0.15:
                S2 = 2.25*deprelmar + 1.3375
            elif deprelmar == -0.15:
                S2 = 1.0
            elif deprelmar < 0.25:
                S2 = -2.25*deprelmar + 0.6625
            else:
                S2 = 0.1
            
            S3 = 0.551*baldcyp/100.0 + 0.713*fresh/100.0 + inter/100.0 + 0.408*brack/100.0
            
            if edge < 22.0:
                S4 = 0.05 + 0.95*(edge/22.0)
            else:
                S4 = 1.0 
            
            if s < 10.0:
                S5 = 1.0 - 0.1*s
            else:
                S5  = 0.0

# no deepwater                
#            if deep < 10.0:
#                S6 = deep/10.0
#            elif deep < 20.0:
#                S6 = 1.0
#            elif deep < 100.0:
#                S6 = 1.25 - 0.0125*deep
#            else:
#                S6 = 0.0
#              
#            HSI_Al = (S1*S2*S3*S4*S5*S6)**(1./6.)
#            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_Al,pwat,deprelmar,baldcyp,fresh,inter,brack,deep,s,edge) 
# no deepwater
            HSI_Al = (S1*S2*S3*S4*S5)**(1./5.)
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,9999,%s,%s\n' %(gridID,HSI_Al,pwat,deprelmar,baldcyp,fresh,inter,brack,s,edge) 
            fAl.write(writestring)

# map alligator HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    del(pwat,deep,deprelmar,baldcyp,fresh,inter,brack,S1,S2,S3,S4,S5,S6)
    

######################################
##         SEASIDE SPARROW HSI      ##
######################################
# 2023 Update -- SES 6/30/20 added new HSI #    
    print ' Calculating Seaside Sparrow HSI'

    HSIcsv = r'%sSPARR.csv' % csv_outprefix
    HSIasc = r'%sSPARR.asc' % asc_outprefix

    with open(HSIcsv,'w') as fSp:
      
        headerstring = 'GridID,HSI,v1a,v1b,v1c,v2,elv\n'
        fSp.write(headerstring)

        for gridID in gridIDs:
 #           w = waterdict[gridID]    # SES 6/30/20 did not set other v1_s from from report because all v1_s multiplied by 0.0        
 #           v1a = frattdict[gridID]  # fresh attached marsh (V1e)
 #           v1b = frfltdict[gridID]  # fresh floating marsh (V1d)
            v1c = interdict[gridID]   # intermed marsh (V1c)
            v1b = brackdict[gridID]  # brackish marsh (V1b)
            v1a = salmardict[gridID] # saline marsh (V1a)
 #           v1f = swfordict[gridID]  # swamp forest (V1g) (same as bottomland - LULC reclass doesn't differentiate between the two)
 #           v1g = btfordict[gridID]  # bottomland forest (V1f) (same as swamp forest - LULC reclass doesn't differentiate between the two)

            elv = melevdict[gridID]   # Eric will have to create marsh elevation dictionanary and outputs for us

            S1 = 1.0*(v1a/100.) + 0.7*(v1b/100.) + 0.3*(v1c/100.)   # divde by 100 to go from percent to proportion veg type for equation
            
#   why is btfordict added to wetlndict?  the equation is from fish hsis for percent marsh (wetland)
            v2 =  max(0.0,min(wetlndict[gridID]+btfordict[gridID],100.0))   # percent of cell that is wetland (includes floating marsh as of 6/30/20)

            if v2 < 65.:                       
                S2 = 0.0154*v2                                        
            else:                       
                S2 = 1.0

            if elv <= 0.09:
                S3 = 0.0
            elif elv < 0.285
                S3 = 5.025*elv - 0.452
            else:
                S3 = 1.0

            HSI_Spar = (S1*S2*S3)**(1./3.)
            writestring = '%s,%s,%s,%s,%s,%s,%s\n' %(gridID,HSI_Spar,v1a,v1b,v1c,v2,elv) 
            fSp.write(writestring)

# map sparrow HSI to Ascii grid
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    del(v1a,v1b,v1c,elv,v2,S1,S2,S3)
    
    
############################
##    Nitrogen Uptake     ##
############################
    print ' Calculating Nitrogen Uptake'

    BLH_Dn = 491.1
    Swampforest_Dn = 288.7
    Freshwatermarsh_Dn = 385.4
    Brackishmarsh_Dn = 404.3
    Saltmarsh_Dn = 130.9
  
    HSIcsv = r'%sNITUP.csv' % csv_outprefix
    HSIasc = r'%sNITUP.asc' % asc_outprefix

    with open(HSIcsv,'w') as fNU:
        fNU.write('GridID, HSI, tmp, tmp_multiplier, sal_JanDec_ave, benthic_Dn, waterdict, frattdict, frfltdict, brackdict, salmardict, swfordict, btfordict, landdict, vegetation_value\n')    
    
        for gridID in gridIDs:
        
            if tmp_JanDec_ave[gridID] <= 0.0:
                tmp_multiplier = 0.0
            elif tmp_JanDec_ave[gridID] <= 10:
                tmp_multiplier = 0.01*tmp_JanDec_ave[gridID]
            elif tmp_JanDec_ave[gridID] <= 20:
                tmp_multiplier = 0.1+(0.9-0.1)*(tmp_JanDec_ave[gridID]-10.)/(20.-10.)
            elif tmp_JanDec_ave[gridID] <= 30:    
                tmp_multiplier = 0.9+(1-0.9)*(tmp_JanDec_ave[gridID]-20.)/(30.-20.)
            else:    
                tmp_multiplier = 1.0
            
            if sal_JanDec_ave[gridID] <= 5.0:
            	benthic_Dn = 134.0
            elif sal_JanDec_ave[gridID] <= 20.0:	
            	benthic_Dn = 56.2
            else:	
            	benthic_Dn = 64.4            
            if waterdict[gridID] >= 0.0: # if water percentage is negative (i.e. -9999), set benthic value to 0
                benthic_value = waterdict[gridID]/100.*500*500*benthic_Dn*tmp_multiplier
            else:
                benthic_value = 0.0
                
            fresh_marsh = (max(0,frattdict[gridID]) + max(0,frfltdict[gridID]) )/100.0
            brack_marsh = max(0,brackdict[gridID])/100.0
            salt_marsh = max(0,salmardict[gridID])/100.0
            swamp_for = max(0,swfordict[gridID])/100.0
            btm_for = max(0,btfordict[gridID])/100.0
            land_area = 500*500*max(0,landdict[gridID])/100.0
                       
            vegetation_value = (fresh_marsh*Freshwatermarsh_Dn + brack_marsh*Brackishmarsh_Dn + salt_marsh*Saltmarsh_Dn + swamp_for*Swampforest_Dn + btm_for*BLH_Dn)*land_area*tmp_multiplier
            
            N_Uptake = max(0.0,vegetation_value + benthic_value)
            
            fNU.write(str(gridID) + ',' + str(N_Uptake) + ',' + str(tmp_JanDec_ave[gridID]) + ',' + str(tmp_multiplier) + ',' + str(sal_JanDec_ave[gridID]) + ',' + str(benthic_Dn) + ',' + str(waterdict[gridID]) + ',' + str(frattdict[gridID]) + ',' + str(frfltdict[gridID]) + ',' + str(brackdict[gridID]) + ',' + str(salmardict[gridID]) + ',' + str(swfordict[gridID]) + ',' + str(btfordict[gridID]) + ',' + str(landdict[gridID]) + ',' + str(vegetation_value) + '\n')
        
    HSIascii_grid(HSIcsv,HSIasc,ascii_grid_lookup,n500cols,n500rows,ascii_header)
    

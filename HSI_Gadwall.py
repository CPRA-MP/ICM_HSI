import numpy as np
import os
import ICM_HelperFunctions as hf


def HSI_Gadwall(year,HSI_dir,csv_outprefix,asc_outprefix,gridIDs,saldict,waterdict,frattdict,frfltdict,interdict,brackdict,salmardict,swfordict,btfordict,watsavdict,map2grid,ascii_grid_lookup,ascii_header,ascii_header_nrows,n500rows,n500cols):
    #######################################
    ###           GADWALL HSI            ##
    #######################################
    
    #months for averaging salinity and temperature data
    months2use = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'octb', 'nov', 'dec']
    
    # This HSI assumes is coded based on the depth tables generated by WM.HSIreclass provide depth in centimeters
    GadwallCSV = os.path.normpath("%s/GadwallDepths_cm_%s.csv" % (HSI_dir,year))

    print( ' Calculating Gadwall HSI')
    Gadwall = np.genfromtxt(GadwallCSV,delimiter = ',',  skip_header = 1)
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

    sal2use = hf.monthly_avg_dict(saldict,months2use)

    with open(HSIcsv,'w') as fG:
        
        headerstring = 'GridID,HSI,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v3i,v3j,v3k,v3l\n'
        fG.write(headerstring)
        
        for gridID in gridIDs:
            s = sal2use[gridID]
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
            
    # 2023 Update -- SES 6/18/20
            S1 = 0.68*(v1a + w_fresh)/100.0 + 0.68*v1b/100.0 + (v1c + w_inter)/100.0 + 0.50*(v1d +  w_brack)/100.0 + 0.09*(v1e + w_saline)/100.0 + 0.05*v1f/100.0 + 0.05*v1g/100.0
            
            v2 = watsavdict[gridID]/100.0     #watsavdict is in percentages, therefore divide by 100 to get V2 in correct unit
            
            if v2 < 0.30:
                S2 = 0.08
            elif v2 >= 0.30 and v2 < 0.70:
                S2 = 2.3*v2 - 0.61
            else:  
                S2 = 1.0
                    
            area = 0.0
            for x in range(1,14):
                area = area + GadwallDepdict[gridID][x]  #determine area of cell (not exactly equal to 500x500 since the 30x30 m grid doesn't fit in the 500x500
                if area < 480*480:
                    area = 480*480
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
                    
                #TODO check if needed for alligator: save deep water from GadwallDepdict for use in Alligator HSI
                deepwat[gridID] = GadwallDepdict[gridID][13]/area        # portion of cell greater than 150 cm deep
                    
            S3 = 0.05*v3a + 0.15*v3b + 0.35*v3c + 0.6*v3d + 0.83*v3e + 1.0*v3f + 0.86*v3g + 0.61*v3h + 0.37*v3i + 0.2*v3j + 0.1*v3k + 0.05*v3l
                    
            HSI_Gadwall = (S1*S2*S3)**(1./3.)
            
            writestring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (gridID,HSI_Gadwall,w,s,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v3i,v3j,v3k,v3l)
            fG.write(writestring)

    # NOTE this was commented out; map2grid flag added (2029MP_Task06)
    hf.HSIascii_grid(HSIcsv,HSIasc,map2grid,ascii_grid_lookup,n500cols,n500rows,ascii_header)

    # delete any dictionaries that aren't used in any other HSIs - frees up memory
    del(GadwallDepdict,Gadwall,GDdict)

    # delete temporary variables so they do not accidentally get used in other HSIs
    del(area,v1a,v1b,v1c,v1d,v1e,v1f,v1g,v2,v3a,v3b,v3c,v3d,v3e,v3f,v3g,v3h,v3i,v3j,v3k,v3l,S1,S2,S3)
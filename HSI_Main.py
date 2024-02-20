#HSI imports
import HSI_PreProcessing as hsi_pp
import HSI_BlueCrab as bc
import HSI_AdultMenhaden as am
import HSI_Gadwall as gad

#ICM imports
import ICM_Settings as icm

#Python imports
import os
import sys

##############################################
##    HABITAT SUITABILITY INDICES ~ HSIs    ##
##############################################

def RunHSI(startyear, endyear):

    #defined variables
    HSI_dir = icm.HSI_dir
    gridIDs = icm.gridIDs
    map2grid = icm.map2grid
    ascii_grid_lookup = icm.ascii_grid_lookup
    ascii_header = icm.ascii_header
    ascii_header_nrows = icm.ascii_header_nrows
    n500rows = icm.n500rows
    n500cols = icm.n500cols

    # change working directory to veg folder
    os.chdir(HSI_dir)
    sys.path.append(HSI_dir)

    for year in range(startyear,endyear+1):

        print('\n--------------------------------------------------')
        print('  RUNNING HABITAT SUITABILITY INDICES - Year %s' % year)
        print('--------------------------------------------------\n')
        
        #unpack returned vars
        csv_outprefix,asc_outprefix,land_mult,fresh_for_mult,bare_mult,saldict,tmpdict,wetlndict,btfordict,swfordict,frattdict,interdict,brackdict,salmardict,watsavdict,frfltdict,baldcypdict,cultchdict,waterdict  = hsi_pp.HSIyearlyVars(year)
        
        bc.HSI_BlueCrab(HSI_dir,csv_outprefix,asc_outprefix,gridIDs,land_mult,fresh_for_mult,bare_mult,saldict,tmpdict,wetlndict,watsavdict,cultchdict,map2grid,ascii_grid_lookup,ascii_header,ascii_header_nrows,n500rows,n500cols )
        
        am.HSI_AdultMenhaden(HSI_dir,csv_outprefix,asc_outprefix,gridIDs,land_mult,fresh_for_mult,bare_mult,saldict,tmpdict,wetlndict,map2grid,ascii_grid_lookup,ascii_header,ascii_header_nrows,n500rows,n500cols )
        
        gad.HSI_Gadwall(year,HSI_dir,csv_outprefix,asc_outprefix,gridIDs,saldict,waterdict,frattdict,frfltdict,interdict,brackdict,salmardict,swfordict,btfordict,watsavdict,map2grid,ascii_grid_lookup,ascii_header,ascii_header_nrows,n500rows,n500cols )


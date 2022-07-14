ecoregions = ['MRP','LPO','LBO','CHS','CHSbi','UBR','LBR','LBRbi','BFD','UBA','MBA','LBAnw','LBAne','LBAnebi','LBAsw','LBAswbi','LBAse','LBAsebi','UVR','VRT','ETB','ETBbi','PEN','WTE','WTEbi','ATB','ATD','TVB','MEL','CHR','CAL','SAB','-9999']
spp = ['ALLIG','BLUCJ','BSHRL','BSHRS','CRAYF_3month','CRAYF_aug2nov','CRAYF_oct2dec','CRAYF_sep2dec','EAGLE','GADWA','GMENA','GMENJ','LMBAS','MOTDK','OYSTE','SPARR','SPSTA','SPSTJ','WSHRL','WSHRS']
years = range(1,53)
Ss = [7,8]
Gs = [500]#range(601,653)


# read in grid-to-compartmant lookup tables
print('Reading in grid-to-compartment lookup data')
grid2comp = {}
grid2comp_file = 'S07/G500/hydro/grid_lookup_500m.csv'
with open(grid2comp_file,mode='r') as g2cf:
    nl = 0
    for line in g2cf:
        if nl > 0:
            grid = int(line.split(',')[0])
            comp = int(line.split(',')[1])
            grid2comp[grid] = comp
        nl += 1


# read in compartmant-to-ecoregion lookup tables
print('Reading in compartment-to-ecoregion data')
comp2eco = {}
comp2eco_file = 'S07/G500/geomorph/input/compartment_ecoregion.csv'
with open(comp2eco_file,mode='r') as c2ef:
    nl = 0
    for line in c2ef:
        if nl > 0:
            comp = int(line.split(',')[0])
            er   = line.split(',')[4]
            comp2eco[comp] = er
        nl += 1
        
        
for s in Ss:
    for g in Gs:
        print('Processing S%02d G%03d:' % (s,g) )
        # build out initial zero dictionary for [species > ecoregion > year > cumulative habitat units]
        hu_er = {}
        for sp in spp:
            hu_er[sp] = {}
            for er in ecoregions:
                hu_er[sp][er] = {}
                for y in years:
                    hu_er[sp][er][y] = 0.0
        
        # read in HSI data from gridded output files and cumulate over each ecoregion
        for sp in spp:
            print('  - %s' % sp)
            for y in years:
                csv = 'S%02d/G%03d/hsi/MP2023_S%02d_G%03d_C000_U00_V00_SLA_O_%02d_%02d_X_%s.csv' % (s,g,s,g,y,y,sp)
                
                with open(csv,mode='r') as hsi:
                    nr = 0
                    for row in hsi:
                        if nr > 0:
                            grid =   int(row.split(',')[0])
                            val  = float(row.split(',')[1])
                            try:
                                er = comp2eco[grid2comp[grid]]
                                hu_er[sp][er][y] += max(val,0.0)        # use 0 filter to remove any NoData/-9999 values
                            except:
                                print('failed on: S%02d G%03d %s %s %s' % (s,g,sp,y,row) )
                        nr += 1
        
        # write output files for each ecoregion
        for er in ecoregions:
            outfile = 'S%02d/G%03d/hsi/MP2023_S%02d_G%03d_C000_U00_V00_%s_O_01_52_X_hsi.csv' % (s,g,s,g,er,sp)
            print('Writing: %s' % outfile)

            with open(outfile,mode='w') as outcsv:
                header = 'ModelGroup,Scenario,Year_FWOA,Year_ICM'
                
                for sp in spp:
                    header = '%s,%s' % (header,sp)
                outcsv.write('%s\n' % header)
                
                for y in years:
                    if y == 1:
                        y_fwoa = -2
                    elif y == 2:
                        y_fwoa = -1
                    else:
                        y_fwoa = y - 2
                        
                    line2write = '%d,%d,%d,%d' % (g,s,y_fwoa,y)
                    for sp in spp:
                        line2write = '%s,%0.4f' % (line2write,hu_er[sp][er][y])
                    
                    outcsv.write('%s\n' % line2write)
                    

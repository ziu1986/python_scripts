import os, glob, sys
import numpy as np
import xarray as xr

def write_grid(lat, lon, outname):
    '''
    Generate grid definition files needed for cdo regriding routine.
    '''
    outfile = open(outname, 'w')
    mesh = np.meshgrid(lat, lon)
    # Write file line by line
    outfile.write("# OsloCTM3 grid file\n")
    outfile.write("gridtype = curvilinear\n")
    outfile.write("gridsize = %s\n" % (lat.size*lon.size))
    outfile.write("xsize = %s\n" % lon.size)
    outfile.write("ysize = %s\n" % lat.size)
    outfile.write("\n")
    
    outfile.write("# Longitudes\n")
    outfile.write("xvals = \n")
    for ilon in np.ravel(mesh[1], order='F'):
        outfile.write("%s  " % ilon)
        
    outfile.write("\n\n")
    
    outfile.write("# Longitudes of cell corners\n")
    outfile.write("xbounds = \n")
    buffering = []
    for ilon in np.arange(lon.size-1):
        delta_lon = np.abs(lon[ilon+1]-lon[ilon])
        buffering.append("%s %s %s %s " % (delta_lon/2+lon[ilon], delta_lon/2+lon[ilon],
                                             -delta_lon/2+lon[ilon], -delta_lon/2+lon[ilon]))
    delta_lon = np.abs(lon[-1]-lon[-2])
    buffering.append("%s %s %s %s " % (delta_lon/2+lon[-1], delta_lon/2+lon[-1],
                                         -delta_lon/2+lon[-1], -delta_lon/2+lon[-1]))
    for ilat in np.arange(lat.size):
        outfile.write("%s \n" % (buffering))
                         
    outfile.write("\n\n")
    
    outfile.write("# Latitudes\n")
    outfile.write("yvals = \n")
    for ilat in np.ravel(mesh[0], order='F'):
        outfile.write("%s  " % ilat)
          
    outfile.write("\n\n")
    
    outfile.write("# Latitudes of cell corners\n")
    outfile.write("ybounds = \n")

    for ilat in np.arange(lat.size-1):
        delta_lat = np.abs(lat[ilat+1]-lat[ilat])
        for ilon in np.arange(lon.size):    
            outfile.write("%s %s %s %s " % (-delta_lat/2+lat[ilat], delta_lat/2+lat[ilat],
                                              delta_lat/2+lat[ilat], -delta_lat/2+lat[ilat]))
    delta_lat = np.abs(lat[-1]-lat[-2])
    for ilon in np.arange(lon.size):    
        outfile.write("%s %s %s %s " % (-delta_lat/2+lat[-1], delta_lat/2+lat[-1],
                                          delta_lat/2+lat[-1], -delta_lat/2+lat[-1]))
    outfile.write("\n\n")
    outfile.close()

nc_src = os.environ['DATA']+'/astra_data/ctm_results/C3RUN_emep_ppgs_2005/scavenging_daily/scavenging_daily_2d_20050101.nc'
# Read the data
try:
    input_resolution
except NameError:
    input_resolution = xr.open_dataset(nc_src)

#write_grid(input_resolution['dry_O3'].lat.data, input_resolution['dry_O3'].lon.data,"osloctm3_grid.txt")

#hardacre_lat = np.arange(-88.5,89,3, dtype=np.float64)
#hardacre_lon = np.arange(1.5,360,3, dtype=np.float64)

#write_grid(hardacre_lat, hardacre_lon,"hardacre_grid.txt")

pft_pct_lat = np.arange(-89.8,90,0.5, dtype=np.float64)
pft_pct_lon = np.arange(-179.8,180,0.5, dtype=np.float64)

write_grid(pft_pct_lat, pft_pct_lon,"pft_pct_grid.txt")


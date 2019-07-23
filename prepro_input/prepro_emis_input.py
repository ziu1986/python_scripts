import os, glob, sys
from scipy.constants import *     # Get physics constants
import datetime as dt
from cdo import *
from mytools.met_tools import *
from mytools.netcdf_tools import *


var = ('NOx_em_anthro', 
'CO_em_anthro', 
'SO2_em_anthro', 
'NH3_em_anthro', 
'VOC02_ethane_em_speciated_VOC_anthro', 
'VOC03_propane_em_speciated_VOC_anthro', 
'VOC04_butanes_em_speciated_VOC_anthro', 
'VOC05_pentanes_em_speciated_VOC_anthro', 
'VOC06_hexanes_pl_em_speciated_VOC_anthro', 
'VOC07_ethene_em_speciated_VOC_anthro', 
'VOC08_propene_em_speciated_VOC_anthro', 
'VOC13_benzene_em_speciated_VOC_anthro', 
'VOC14_toluene_em_speciated_VOC_anthro', 
'VOC15_xylene_em_speciated_VOC_anthro', 
'VOC16_trimethylb_em_speciated_VOC_anthro', 
'VOC17_other_arom_em_speciated_VOC_anthro', 
'VOC21_methanal_em_speciated_VOC_anthro', 
'VOC22_other_alka_em_speciated_VOC_anthro', 
'VOC23_ketones_em_speciated_VOC_anthro', 
'BC_em_FOSSIL_FUEL_anthro', 
'BC_em_SOLID_BIOFUEL_anthro', 
'OC_em_FOSSIL_FUEL_anthro', 
'OC_em_SOLID_BIOFUEL_anthro')

data_dir = os.environ['DATA']+'/CTM3_input_data/EMIS/CEDS_SCENARIOS/SSP585/'
data_tgt = os.environ['DATA']+'/astra_data/input_data/ctm_input/EMIS/CEDS_SCENARIOS/SSP585_CICERO/'
cicero_sectornames = {0:"AGR", 1:"ENE", 2:"IND", 3:"TRA", 4:"RBC", 5:"SOL", 6:"WAS", 7:"SHI"}

for ivar in var:
    file_name = ivar.replace('_','-')+"_input4MIPs_emissions_ScenarioMIP_IAMC-REMIND-MAGPIE-ssp585*.nc"
    file = glob.glob(data_dir+file_name)[0]
    # Read the data
    print("Reading file %s" % (file))
    tgt_file_name = os.path.basename(file[:-3]+'_CICERO.nc')
    data = xr.open_dataset(file)
    # Copy the original and drop sector info
    output = data.drop((ivar,'sector','sector_bnds')).copy()
    # Add grid cell area?
    #cdo = Cdo()
    #cdo.gridarea(input=data_dir+file_name, options = "-f nc",  returnArray = 'P')

    # Seperate the sectors
    for isector in data['sector'].values:
        output[ivar+'_'+cicero_sectornames[isector]] = data[ivar].sel(sector=isector)
      
    # Info
    output.attrs['CICERO_title'] = "CEDS sectors separated for CICERO usage" 
    output.attrs['CICERO_history'] = "python script by S.Falk, MetOS" 
    output.attrs['CICERO_date'] = dt.datetime.now().isoformat()
    output.attrs['title'] = "Annual Anthropogenic Emissions of BC prepared for CMIP6" 
    output.attrs['source_id'] = "CEDS-2017-05-18" 
    output.attrs['references'] = "MGidden, M. J. , Riahi, K., Smith, S. J. , Fujimori, S., Luderer, G., Kriegler, E., van Vuuren, D. P., van den Berg, M., Feng, L., Klein, D., Calvin, K., Doelman, J. C., Frank, S., Fricko, O., Harmsen, M., Hasegawa, T., Havlik, P., Hilaire, J., Hoesly, R., Horing, J., Popp, A., Stehfest, E., Takahashi, K.: Global emissions pathways under different socioeconomic scenarios for use in CMIP6: a dataset of harmonized emissions trajectories through the end of the century, Geosci. Model Dev., doi:10.5194/gmd-12-1443-2019, 12 April 2019." 
    output.attrs['dataset_version_number'] = dt.datetime.now().isoformat()[:10]
    output.attrs['reporting_unit'] = "Mass flux of %s" % (ivar[:ivar.find('_')])
    output.attrs['frequency']  = "mon"

    # Write it to file
    output.to_netcdf(data_tgt+tgt_file_name)
    print("Written data to %s" % (data_tgt+tgt_file_name))

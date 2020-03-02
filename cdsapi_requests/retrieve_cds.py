import cdsapi
import calendar
import os
import shutil
from optparse import OptionParser
def main():
    usage = "usage: %prog --start_year YYYY --end_year YYYY [--start_month MM] [--end_month MM] [--start_day DD] [--end_day DD] [--area lat_start/lon_sstart/lat_end/lon_end] [--path]"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--start_year", dest="start_year",
                      help="start year YYYY", metavar="start_year",type=int )
    parser.add_option("-e", "--end_year", dest="end_year",
                      help="end_year YYYY", metavar="end_year", type=int)
    parser.add_option("--start_month", dest="start_month",
                      help="start month MM", metavar="start_month", type=int)
    parser.add_option("--end_month", dest="end_month",
                      help="end month DD", metavar="end_month", type=int)
    parser.add_option("--start_day", dest="start_day",
                      help="start day DD", metavar="start_day", type=int)
    parser.add_option("--end_day", dest="end_day",
                      help="end day DD", metavar="end_day", type=int)
    parser.add_option("--area", dest="area",
                      help="'lat_start/lon_start/lat_end/lon_end'", metavar="area", type=str)
    parser.add_option("--path", dest="path",
                      help="'output directory'", metavar="path", type=str)
    (options, args) = parser.parse_args()
    if not options.start_year:
        parser.error("start year must be specified!")
    else:
        start_year=options.start_year
    if not options.end_year:
        end_year=start_year
    else:
        end_year=options.end_year
    if not options.start_month:
        start_month=1
    else:
        start_month=options.start_month
    if not options.end_month:
        end_month=12
    else:
        end_month=options.end_month
    if not options.area:
        area = ""
    else:
        area = options.area
    if not options.path:
        path = "."
    else:
        path = options.path

    server = cdsapi.Client()
    
    print(start_year)
    print(end_year)
    for year in range(start_year, end_year+1):
        print('YEAR ',year)
        for month in range(start_month,end_month+1):
            if not options.start_day:
                sdate="%s%02d01"%(year,month)
            else:
                sdate="%s%02d%02d"%(year,month,int(options.start_day))
            if not options.end_day:
                lastday=calendar.monthrange(year,month)[1]
                edate="%s%02d%s"%(year,month,lastday)
            else:
                edate="%s%02d%02d"%(year,month,int(options.end_day))
            print('get data from %s/to/%s (YYYYMMDD)' % (sdate,edate))
            print('saving data to %s' % path)
            print('saving data to file Tair_WFDE5_CRU_%s.nc' % sdate[:-2])
            if area=="":
                server.retrieve(
                    'derived-near-surface-meteorological-variables',
                    {
                        'format': 'nc',
                        'variable': 'near_surface_air_temperature',
                        'reference_dataset': 'cru',
                        'year': '%d' % year,
                        'month': '%02d' % month,
                    },
                    '%s/Tair_WFDE5_CRU_%s.nc' % (path, sdate[:-2]))
            else:
                print(area)
                server.retrieve(
                    'derived-near-surface-meteorological-variables',
                    {
                        'format': 'nc',
                        'variable': 'near_surface_air_temperature',
                        'reference_dataset': 'cru',
                        'year': '%d' % year,
                        'month': '%02d' % month,
                        'area': "%s" % (area),
                    },
                    '%s/Tair_WFDE5_CRU_%s.nc' % (path, sdate[:-2]))
            
if __name__ == "__main__":
    main()

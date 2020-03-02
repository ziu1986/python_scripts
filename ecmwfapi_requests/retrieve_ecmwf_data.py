#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import calendar
import os
import shutil
from optparse import OptionParser
def main():
    usage = "usage: %prog --start_year YYYY --end_year YYYY [--start_month MM] [--end_month MM] [--start_day DD] [--end_day DD]"
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
    server = ECMWFDataServer()
    print start_year
    print end_year
    for year in range(start_year, end_year+1):
        print 'YEAR ',year
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
            print 'get data from ', sdate,' to ',edate,' (YYYYMMDD)'
            server.retrieve({
                    'dataset'  : "interim",
                    'date'     : "%s/to/%s"%(sdate,edate),
                    'time'     : "00/06/12/18",
                    'step'     : "00",
                    'stream'   : "oper",
                    'levtype'  : "sfc",
                    'type'     : "an",
                    'class'    : "ei",
                    'grid'     : "0.5/0.5",
                    'param'    : "165/166/167/168/134/151/235/31/34/33/141/139/170/183/236/39/40/41/42",
                    #'area'     : "70./-9./46./37.",
                    'target'   : "an_sfc_%s_%s.grb"%(sdate,edate),
                    })
            server.retrieve({
                    'dataset'  : "interim",
                    'date'     : "%s/to/%s"%(sdate,edate),
                    'time'     : "00/06/12/18",
                    'step'     : "00",
                    'stream'   : "oper",
                    'levtype'  : "pl",
                    'levelist' : "all",
                    'type'     : "an",
                    'class'    : "ei",
                    'grid'     : "0.5/0.5",
                    'param'    : "129/130/131/132/133/157/156",
                    'frame'    : "OFF",
                    #'area'     : "70./-9./46./37.",
                    'target'   : "an_pl_%s_%s.grb"%(sdate,edate),
                    })
if __name__ == "__main__":
    main()

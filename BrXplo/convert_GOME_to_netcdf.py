import os, glob
execfile('read_GOME_data.py')

# Access environment variable for directory
#data_dir = os.environ['DATA']
#subd = '/BrO/'
src_dir = '../../gome_data/'
for i in range(1,13):
    if i < 10:
        src = 'gome_bro_000%s*_v1*' % (i)
    else:
        src = 'gome_bro_00%s*_v1*' % (i)
    BrO_data = []

    for file in sorted(glob.glob(src_dir+src)):
        print file
        if file.find("asp")>-1:
            print("GOME asp-file")
            BrO_data.append(read_gome_data(file, l_header=24))
        elif file.find("asc")>-1:
            print("GOME asc-file")
            BrO_data.append(read_gome_data(file, l_header=22))
    if len(BrO_data)!=0:
        BrO_con = xr.concat(BrO_data, dim='time')
        BrO_con.to_netcdf(src[:src.find('*')]+'.nc')
        print("Wrote", src[:src.find('*')]+'.nc')

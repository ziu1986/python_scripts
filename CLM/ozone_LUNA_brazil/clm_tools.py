def load_data(src):
    import glob, os, sys
    import xarray as xr

    data_list = []
    for file in sorted(glob.glob(src)):
        print("Loading... %s" % (file))
        data = xr.open_dataset(file)
        data_list.append(data[['GSSHA', 'GSSUN', 'JMX25T', 'Jmx25Z', 'PSNSHA', 'PSNSUN', 'RSSHA','RSSUN', 'VCMX25T', 'Vcmx25Z', 'TOTVEGC', 'TOTVEGN', 'NPP', 'GPP']])
    data = xr.concat(data_list, dim='time')
    return(data)

def execfile(filepath, globals=None, locals=None):
    '''
    Re-implementation of python2 interpreters execfile.

    - Uses binary reading to avoid encoding issues
    - Guaranteed to close the file (Python3.x warns about this)
    - Defines __main__, some scripts depend on this to check if they are loading as a module or not for eg. if __name__ == "__main__"
    - Setting __file__ is nicer for exception messages and some scripts use __file__ to get the paths of other files relative to them.

    - Takes optional globals & locals arguments, modifying them in-place as execfile does - so you can access any variables defined by reading back the variables after running.

    - Unlike Python2's execfile this does not modify the current namespace by default. For that you have to explicitly pass in globals() & locals().

    '''
    if globals is None:
        globals = {}
    globals.update({
        "__file__": filepath,
        "__name__": "__main__",
    })
    with open(filepath, 'rb') as file:
        exec(compile(file.read(), filepath, 'exec'), globals, locals)


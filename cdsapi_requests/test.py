import cdsapi

c = cdsapi.Client()

c.retrieve(
    'derived-near-surface-meteorological-variables',
    {
        'format': 'zip',
        'variable': 'near_surface_air_temperature',
        'reference_dataset': 'cru',
        'year': '2018',
        'month': [
            '06', '07', '08',
        ],
    },
    'download.zip')

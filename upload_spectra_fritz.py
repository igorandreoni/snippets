import requests
import glob
import re

from astropy.io import ascii

# Copy the API token from your Fritz account
token = 'copied-from-fritz'

# Default observers and reducers
# users IDs can be found at: https://fritz.science/api/user
default_observer = [37]
default_reducer = [14]


def api(method, endpoint, data=None):
    """API request"""

    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response


def upload_spectrum(spec, observers, reducers, group_ids=[], date=None,
                    inst_id=3, ztfid=None, meta=None):
    """
    Upload a spectrum to the Fritz marshal

    ----
    Parameters

    spec astropy table object
        table with wavelength, flux, fluxerr
    observers list of int
        list of integers corresponding to observer usernames
    reducers list of int
        list of integers corresponding to reducer usernames
    group_ids list of int
        list of IDs of the groups to share the spectum with
    date str
         date (UT) of the spectrum
    inst_id int
        ID of the instrument (default DBSP)
    ztfid str
        ID of the ZTF source, e.g. ZTF20aclnxgz
    """

    data = {
            "observed_by": observers,
            "group_ids": group_ids,
#            "assignment_id": 0,
            "altdata": meta,
            "observed_at": str(date),
            "fluxes": list(spec['flux']),
            "errors": list(spec['fluxerr']),
#            "followup_request_id": 0,
            "wavelengths": list(spec['wavelength']),
            "instrument_id": inst_id,
            "reduced_by": reducers,
            "origin": "",
            "obj_id": ztfid
            }

    response = api('POST',
                   'https://fritz.science/api/spectrum',
                   data=data)

    print(f'HTTP code: {response.status_code}, {response.reason}')
    if response.status_code in (200, 400):
        print(f'JSON response: {response.json()}')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Upload spectra to Fritz.')
    parser.add_argument('ID', nargs='+', help='Spectra filenames')
    parser.add_argument('-d', action='store_true',
                        help='Use default reducer ID and observer - \
please customize the code')
    parser.add_argument('--date', dest='date', type=str,
                        required=True, help='Date of the observations, \
e.g. 2020-11-10T00:00:00', default='2020-11-11T00:00:00')
    parser.add_argument('--inst', dest='inst_id', type=int,
                        required=False, help='Instrument ID, \
e.g. inst_id = 3 for DBSP. Instrument IDs can be found here: \
https://fritz.science/api/instrument', default=3)
    args = parser.parse_args()

    # Observers, reducers
    if args.d is True:
        #default observer
        observers = default_observer
        #default reducer
        reducers = default_reducer
    else:
        print("users IDs can be found at: https://fritz.science/api/user")
        observers = input("enter comma-separated IDs (int) \
of the observers:\n")
        observers = observers.split(",")
        reducers = input("enter comma-separated IDs (int) of the reducers:\n")
        reducers = reducers.split(",")

    # For each ID, check which files are available
    for source_filename in args.ID:
        files = glob.glob(f"{source_filename}")
        if len(files) == 0:
            print(f"No files named {source_filename} found")
            filename = input(f"Enter the correct spectrum file name \
or enter 'c' to skip this source and continue: \n")
        else:
            filename = files[0] 
        if filename == 'c':
            continue
        elif filename[-4:] == 'fits':
            print("Reading of FITS files is not yet implemented, \
please enter the name of an ascii file")
            continue
        # Read the file
        spec = ascii.read(filename)
        spec.rename_column("col1", "wavelength")
        spec.rename_column("col2", "flux")
        # Uncertainty for DBSP
        if len(spec.colnames) > 2 and args.inst_id == 3:
            spec.rename_column("col3", "fluxerr")
        elif len(spec.colnames) > 2 and args.inst_id == 7:
            spec.rename_column("col4", "fluxerr")
        elif len(spec.colnames) == 2:
            spec["fluxerr"] = np.zeros(len(spec))
        # Metadata
        meta = spec._meta['comments']
        # Extract the source filename
        if not "ZTF" in source_filename:
            source = input(f"No 'ZTF' found in the file name, please enter \
the name of the source for the spectrum {source_filename}:\n")
        else:
            span = re.search("ZTF", source_filename).span()
            source = source_filename[span[0]:span[0]+12]
        print(f"Uploading spectrum {source_filename} for source {source}")
        upload_spectrum(spec, observers, reducers, group_ids=[], date=args.date,
                    inst_id=args.inst_id, ztfid=source, meta=meta)

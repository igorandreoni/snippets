import requests

from astropy.io import ascii
from astropy.table import Table, unique
from astropy.time import Time
import numpy as np


from penquins import Kowalski

email = "<your email address>"
userpass = "<your password>"


def get_lightcurve_alerts_aux(token, list_names):
    """Query the light curve for a list of candidates"""

    #k = Kowalski(username=username, password=password, verbose=False)
    k = Kowalski(token=token, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts_aux",
                   "filter": {
                              '_id': {'$in': list(list_names)}
                              },
                   "projection": {}
                       },
         "kwargs": {"hint": "_id_"}
         }

    r = k.query(query=q)
    #for ii, (kk, vv) in enumerate(indexes.items()):
    #    print(f'index #{ii+1}: "{kk}"\n{vv["key"]}\n')
    if r['default']['data'] == []:
        print("No candidates to be checked?")
        return None
    out = []
    for l in r['default']['data']:
        with_det = list({'objectId': l['_id'], 'candidate': s} for s in \
l['prv_candidates'] if 'magpsf' in s.keys())
        out = out + with_det

    return out


def get_lightcurve_alerts(token, list_names):
    """Query the light curve for a list of candidates"""
    #k = Kowalski(username=username, password=password, verbose=False)
    k = Kowalski(token=token, verbose=False)
    q = {"query_type": "find",
         "query": {
                   "catalog": "ZTF_alerts",
                   "filter": {
                              'objectId': {'$in': list(list_names)}
                              },
                   "projection": {
                                  "objectId": 1,
                                  "candidate.jd": 1,
                                  "candidate.ra": 1,
                                  "candidate.dec": 1,
                                  "candidate.magpsf": 1,
                                  "candidate.fid": 1,
                                  "candidate.sigmapsf": 1,
                                  "candidate.programid": 1,
                                  "candidate.magzpsci": 1,
                                  "candidate.magzpsciunc": 1,
                                  "candidate.sgscore1": 1,
                                  "candidate.sgscore2": 1,
                                  "candidate.sgscore3": 1,
                                  "candidate.distpsnr1": 1,
                                  "candidate.distpsnr2": 1,
                                  "candidate.distpsnr3": 1,
                                  "candidate.field": 1,
                                  "candidate.rcid": 1,
                                  "candidate.pid": 1
                                  }
                       },
         "kwargs": {"hint": "objectId_1"}
         }

    r = k.query(query=q)
    if r['default']['data'] == []:
        print("No candidates to be checked?")
        return None

    return r['default']['data']


def create_tbl_lc(light_curves, outfile=None):
    """Create a table with the light curves
    and write a CSV output file"""

    # fid -> filter
    filters = {'1': 'g', '2': 'r', '3': 'i'}

    tbl = Table([[], [], [], [], [], [], [], [], [], [], [], [], [], [], [],
                 [], [], [], []],
                names=('name', 'ra', 'dec', 'jd', 'magpsf', 'sigmapsf',
                       'filter', 'magzpsci', 'magzpsciunc',
                       'programid', 'field', 'rcid', 'pid',
                       'sgscore1', 'sgscore2', 'sgscore3',
                       'distpsnr1', 'distpsnr2', 'distpsnr3'),
                dtype=('S12', 'double', 'double', 'double',
                       'f', 'f', 'S', 'f', 'f', 'i', 'i', 'i', 'int_',
                       'f', 'f', 'f', 'f', 'f', 'f'))

    for l in light_curves:
        magzpsci = l["candidate"].get("magzpsci")
        magzpsciunc = l["candidate"].get("magzpsciunc")
        try:
            row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"], l["candidate"]["field"],
               l["candidate"]["rcid"], l["candidate"]["pid"], 
               l["candidate"]["sgscore1"], l["candidate"]["sgscore2"],
               l["candidate"]["sgscore3"], l["candidate"]["distpsnr1"],
               l["candidate"]["distpsnr2"], l["candidate"]["distpsnr3"]]
        except KeyError:
            row = [l["objectId"], l["candidate"]["ra"], l["candidate"]["dec"],
               l["candidate"]["jd"], l["candidate"]["magpsf"],
               l["candidate"]["sigmapsf"], filters[str(l["candidate"]["fid"])],
               magzpsci, magzpsciunc,
               l["candidate"]["programid"], l["candidate"]["field"],  
               l["candidate"]["rcid"], l["candidate"]["pid"], np.nan,
               np.nan, np.nan, np.nan, np.nan, np.nan]
        tbl.add_row(row)
    # Remove exact duplicates
    tbl = unique(tbl)
    tbl.sort("jd")
    #tbl.write(outfile, format='csv', overwrite=True)

    return tbl


if __name__ == "__main__":
    """Trigger Frank Masci's forced photometry service at IPAC"""

    import argparse

    parser = argparse.ArgumentParser(description="Trigger Frank Masci's \
forced photometry service at IPAC")
    parser.add_argument('-n', dest='names', nargs='+', required=False,
                        help='Names of the ZTF candidates; if given, \
the coordinates will be queried from kowalski', default=None)
    parser.add_argument('--ra', dest='ra', type=float, required=False,
    help='Right ascension (array, in degrees)', default=None)
    parser.add_argument('--dec', dest='dec', type=float, required=False,
    help='Declination (array, in degrees)', default=None)

    args = parser.parse_args()

    if args.names is None and (args.ra is None and args.dec is None):
        print("No input candidates. Please use --n and provide ZTF names \
or --ra ra1 ra2] --dec [dec1 dec2] to trigger by coordinates]")
        exit()

    if args.ra is not None and args.dec is not None: 
        if type(args.ra) != list and type(args.dec) != list:
            ra, dec = [args.ra], [args.dec]
        else:
            ra, dec = args.ra, args.dec
        outname = f"script_ra{ra}_dec{dec}.txt"
    else:
        outname = f"script_{'_'.join(args.names)}.txt"

    # Read the secrets
    secrets = ascii.read('/Users/igor/.secrets.csv', format='csv')
    username_kowalski = secrets['kowalski_user'][0]
    password_kowalski = secrets['kowalski_pwd'][0]

    # Get the token
    data={
          "username": username_kowalski,
          "password": password_kowalski
          }
    response = requests.request("POST",
                                'https://kowalski.caltech.edu/api/auth',
                                json=data)
    token = response.json()['token']

    # Get the light curves
    #light_curves_alerts = get_lightcurve_alerts(username_kowalski, password_kowalski, args.names)
    light_curves_alerts = get_lightcurve_alerts(token, args.names)
    
    # Add prv_candidates photometry to the light curve
    #light_curves_aux = get_lightcurve_alerts_aux(username_kowalski, password_kowalski, args.names)

    #light_curves = light_curves_alerts + light_curves_aux

    # Create a table and output CSV file
    t = create_tbl_lc(light_curves_alerts)
    
    # Get coords
    ra, dec = np.median(np.array(t['ra'])), np.mean(np.array(t['dec']))

    # Command
    jdstart = Time("2018-03-17").jd
    jdend = Time.now().jd
    command = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt \
'https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={ra}&\
dec={dec}&jdstart={jdstart}&jdend={jdend}&email={email}&userpass={userpass}'"
    print(command)
 
    with open(outname, "w") as o:
        o.write(command)

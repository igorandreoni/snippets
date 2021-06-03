# Author: Igor Andreoni
# ZTF alert cone search using Kowalski, returns summary table

import glob
import requests
import datetime

import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time

from penquins import Kowalski

def api(method, endpoint, data=None):
    """API request"""

    headers = {'Authorization': f'{token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)

    return response

def query_kowalski_cone(ra, dec, search_rad, jd_start=2458194, jd_end=2459406,
                        programid_list=[1], drb_min=0.64, ndethist_min=2):
    """
    Cone search of ZTF alerts using Kowalski

    Parameters
    ----------
    ra float
        Right ascension (equatorial)
    dec float
        Declination (equatorial)
    jd_start float
        JD of the start of the search
    jd_end float
        JD of the end of the search
    programid_list list of int
        list of ZTF program IDs
    drb_min
        min deep real/bogus score
    ndethist_min int
        min number of detections >3 sigma

    Returns
    -------
    info dict
        query output
    """

    q = {
    "query_type": "cone_search",
    "query": {
              "object_coordinates": {
                                     "cone_search_radius": search_rad,
                                     "cone_search_unit": "arcmin",
                                     "radec": {
                                               "object1": [
                                                           ra,
                                                           dec
                                                           ]
                                               }
                                     },
              "catalogs": {
                           "ZTF_alerts": {
                                          "filter": {'candidate.jd': {'$gt': jd_start, '$lt': jd_end},
                                                     'candidate.drb': {'$gt': drb_min},
                                                     'candidate.programid': {'$in': programid_list},
                                                     'candidate.ndethist': {'$gte': ndethist_min},
                                                     },
                                          "projection": {
                                                         "objectId": 1,
                                                         "candidate.ra": 1,
                                                         "candidate.dec": 1,
                                                         "candidate.jd": 1,
                                                         "candidate.drb": 1,
                                                         "candidate.ndethist": 1,
                                                         "candidate.jdstarthist": 1,
                                                         "candidate.jdendhist": 1,
                                                         "candidate.programid": 1,
                                                         "candidate.isdiffpos": 1,
                                                         "candidate.distpsnr1": 1,
                                                         "candidate.sgscore1": 1,
                                                         "candidate.srmag1": 1,
                                                         "candidate.distpsnr2": 1,
                                                         "candidate.sgscore2": 1,
                                                         "candidate.srmag2": 1,
                                                         "candidate.distpsnr3": 1,
                                                         "candidate.sgscore3": 1,
                                                         "candidate.srmag3": 1
                                                         }
                                          }
                           }
             },
    "kwargs": {
#"hint": "jd_field_rb_drb_braai_ndethhist_magpsf_isdiffpos"
               }
       }

    response = api('POST',
                   f'https://kowalski.caltech.edu/api/queries',
                   data=q)

    print(f"Using a search radius of {search_rad} arcmin")
    if response.status_code == 200:
        info = response.json()['data']['ZTF_alerts']['object1']
        return info
    else:
        print("FAILED kowalski query!")
        print(f'HTTP code: {response.status_code}, {response.reason}')
        print(f'{response.content}')
        exit()


if __name__ == '__main__':
    import argparse 

    parser = argparse.ArgumentParser(description='ZTF alert cone search \
using Kowalski')
    parser.add_argument('radec', metavar='RA, Dec', type=str, nargs='+',
                        help='RA and Dec (degrees)')
    parser.add_argument('-r', dest='radius', type=float,
                        required=False, help='Search radius (arcmin)',
                        default=1)
    parser.add_argument('--date-start', dest='date_start', type=str,
                        required=False,
                        help="Start date of the query, in ISO format. \
                        Example: '2017-08-17 12:41:04.4'", default=None)
    parser.add_argument('--date-end', dest='date_end', type=str,
                        required=False,
                        help="End date of the query, in ISO format. \
                        Example: '2017-08-18 12:00:00.0'", default=None),
    parser.add_argument('--ndethist', dest='ndethist_min', type=int,
                        required=False,
                        help='Minimum number of detections', default=2)
    parser.add_argument('--pid', dest='programid_list', nargs='+',
                        required=False,
                        help='ZTF program IDs', default=[1,3])
    parser.add_argument('--drb', dest='drb_min', type=float,
                        required=False,
                        help='Minimum drb score', default=0.64)
    parser.add_argument('--out', dest='out', type=str,
                        required=False,
                        help='Output filename: if given, a CSV file will be\
                        created', default=None)
    args = parser.parse_args()

    # RA and Dec
    ra, dec = float(args.radec[0]), float(args.radec[1])

    # Radius
    search_rad = args.radius

    # Dates
    if args.date_start is None:
        date_start = Time.now() - datetime.timedelta(days=1)
    else:
        try:
            date_start = Time(args.date_start, format='iso')
        except ValueError:
            print("Invalid start date. It must be a string in ISO format.")
            print("Example: '2017-08-17 12:41:04.4'")
            exit()

    if args.date_end is None:
        date_end = Time.now()
    else:
        try:
            date_end = Time(args.date_end, format='iso')
        except ValueError:
            print("Invalid end date. It must be a string in ISO format.")
            print("Example: '2018-01-01 12:41:04.4'")
            exit()

    # Program IDs
    programid_list = [int(p) for p in args.programid_list]

    # Read the secrets
    secrets = ascii.read('secrets.csv', format='csv')
    username = secrets['kowalski_user'][0]
    password = secrets['kowalski_pwd'][0]

    # Get the token
    data={
          "username": username,
          "password": password
          }
    response = requests.request("POST",
                                'https://kowalski.caltech.edu/api/auth',
                                json=data)
    token = response.json()['token']

    # Query Kowalski
    info_clean = query_kowalski_cone(ra, dec, search_rad, jd_start=date_start.jd, jd_end=date_end.jd,
                        programid_list=programid_list, drb_min=args.drb_min, ndethist_min=args.ndethist_min)

    # Object IDs
    object_ids = set([c['objectId'] for c in info_clean])
    if len(object_ids) == 0:
        print("No candidates in the 'cleaned' alert stream!")
        exit()

    # Initialize dict for summary info
    summary_info = {}

    for object_id in object_ids:
        summary_info[object_id] = {}
        summary_info[object_id]["ra"] = np.mean([c['candidate']['ra'] for c in info_clean if c['objectId']==object_id])
        summary_info[object_id]["dec"] = np.mean([c['candidate']['dec'] for c in info_clean if c['objectId']==object_id])
        summary_info[object_id]["ndethist"] = np.max([c['candidate']['ndethist'] for c in info_clean if c['objectId']==object_id])
        summary_info[object_id]["ndethist"] = np.max([c['candidate']['ndethist'] for c in info_clean if c['objectId']==object_id])
        summary_info[object_id]["jdstarthist"] = np.min([c['candidate']['jdstarthist'] for c in info_clean if c['objectId']==object_id])
        summary_info[object_id]["jdendhist"] = np.max([c['candidate']['jdendhist'] for c in info_clean if c['objectId']==object_id])
        summary_info[object_id]["max_programid"] = np.max([c['candidate']['programid'] for c in info_clean if c['objectId']==object_id])
        diff_pos = [c['candidate']['isdiffpos'] for c in info_clean if c['objectId']==object_id]
        if "t" in diff_pos:
            summary_info[object_id]["includes_diff_pos"] = 1
        else:
            summary_info[object_id]["includes_diff_pos"] = 0
        if "f" in diff_pos:
            summary_info[object_id]["includes_diff_neg"] = 1
        else:
            summary_info[object_id]["includes_diff_neg"] = 0
        # Find one alert for crossmatches
        for alert in info_clean:
            if alert['objectId']==object_id:
                summary_info[object_id]["distpsnr1"] = alert['candidate']['distpsnr1']
                summary_info[object_id]["sgscore1"] = alert['candidate']['sgscore1']
                summary_info[object_id]["srmag1"] = alert['candidate']['srmag1']
                summary_info[object_id]["distpsnr2"] = alert['candidate']['distpsnr2']
                summary_info[object_id]["sgscore2"] = alert['candidate']['sgscore2']
                summary_info[object_id]["srmag2"] = alert['candidate']['srmag2']
                summary_info[object_id]["distpsnr3"] = alert['candidate']['distpsnr3']
                summary_info[object_id]["sgscore3"] = alert['candidate']['sgscore3']
                summary_info[object_id]["srmag3"] = alert['candidate']['srmag3']
                break

    # Create a nice table
    names = ['objectId'] + list(summary_info[list(summary_info.keys())[0]].keys())
    columns = []
    for k in summary_info[list(summary_info.keys())[0]].keys():
        columns.append([summary_info[object_id][k] for object_id in object_ids])
    columns = [list(object_ids)] + columns

    tab = Table(columns, names=names)
    print(f"Done!")
    print(tab)
    if args.out is not None:
        outname = args.out
        tab.write(outname, format='csv', overwrite=True)
        print(f"Created output file: {outname}")

    # Print and say goodbye
    print(f"Found {len(tab)} objects in ZTF data")

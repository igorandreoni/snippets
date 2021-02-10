#Author: Igor Andreoni
#email: andreoni@caltech.edu

import requests
import glob
import re
import numpy as np
from astropy.io import ascii

# Copy the API token from your Fritz account
token = 'copied-from-fritz'

def api(method, endpoint, data=None):
    """API request"""

    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1', 'Yes', 'True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', 'No', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def check_source_exists(source):
    """check if the source exists"""

    response = api('HEAD', f'https://fritz.science/api/sources/{source}')

    if response.status_code == 200:
        print(f"Source {source} was found on Fritz")
    else:
        print(f"Source {source} does not exist on Fritz!")

    return response


def get_groups(source):
    """Get the groups a source belongs to"""

    response = api('GET',
                   f'https://fritz.science/api/sources/{source}'
                   )

    if response.status_code == 200:
        groups = response.json()['data']['groups']
    else:
        print(f'HTTP code: {response.status_code}, {response.reason}')

    return groups


def get_candidate(name):
    """
    Get a candidate from the Fritz marshal

    ----
    Parameters

    name str
        source ZTF name
    """

    response = api('GET',
                   f'https://fritz.science/api/candidates/{name}')

    print(f'HTTP code: {response.status_code}, {response.reason}')
    if response.status_code in (200, 400):
        print(f'JSON response: {response.json()}')

    return response


def get_alerts(name):
    """
    Get the alerts from the Fritz marshal

    ----
    Parameters

    name str
        source ZTF name
    """

    response = api('GET',
                   f'https://fritz.science/api/alerts/ztf/{name}')

    print(f'HTTP code: {response.status_code}, {response.reason}')
    if response.status_code == 400:
        print(f'JSON response: {response.json()}')

    return response.json()


def add_source(name, group_ids):
    """
    Add a new ZTF source to the db

    ----
    Parameters

    name str
        source ZTF name
    group_ids list of int
        list of group IDs
    """

    # Compute precise coordinates
    alerts = get_alerts(name)
    ra = np.median([a['candidate']['ra'] for a in alerts['data']])
    dec = np.median([a['candidate']['dec'] for a in alerts['data']])
   
    data = {"ra": ra,
            "dec": dec,
            "id": name,
            "group_ids": group_ids
            }
   
    response = api('POST',
                   'https://fritz.science/api/sources',
                   data=data)
   
    print(f'HTTP code: {response.status_code}, {response.reason}')
    if response.status_code in (200, 400):
        print(f'JSON response: {response.json()}')


def upload_source(name, group_ids):
    """
    Upload a list of sources to the Fritz marshal

    ----
    Parameters

    name str
        source ZTF name
    group_ids list of int
        list of IDs of the groups to share the spectum with
    """

    data = {
            "objId": name,
            "inviteGroupIds": group_ids,
            "unsaveGroupIds": []
            }

    response = api('POST',
                   'https://fritz.science/api/source_groups',
                   data=data)

    print(f'HTTP code: {response.status_code}, {response.reason}')
    if response.status_code in (200, 400):
        print(f'JSON response: {response.json()}')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Upload sources to Fritz.')
    parser.add_argument('-n', dest='ID', nargs='+', help='Object names')
    parser.add_argument('-f', dest='filename', type=str,
                        required=False, default=None,
                        help="CSV file including a 'name' column"
                        )
    parser.add_argument('-g', dest='groups', nargs='+',
                        required=False, help='IDs of the group to save the\
sources to. Group IDs can be found here: \
https://fritz.science/api/groups\
e.g.: SNLensing = 236')

    args = parser.parse_args()

    if args.ID is not None:
        sources = args.ID
    elif args.filename is not None:
        sources = list(ascii.read(args.filename, format='csv')['name'])
    else:
        print("Provide source names or a CSV filename")
        exit()

    # From str to int
    my_groups = [int(g) for g in args.groups]

    # For each ID, check which files are available
    for source in sources:
        # Does the source exist?
        response = check_source_exists(source)

        # If the source does not exist
        if response.status_code != 200:
            # Add the source and continue
            add_source(source, my_groups)
            continue

        groups = get_groups(source)
        group_ids = [g['id'] for g in groups]
        # Which of my groups still needs the source to be saved?
        missing_groups = [g for g in my_groups
                          if not (int(g) in group_ids)]
        print("Missing groups:", missing_groups)
        if len(missing_groups) > 0:
            # Save the source to the desired groups
            upload_source(source, missing_groups)

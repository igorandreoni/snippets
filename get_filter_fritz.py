# Author: Igor Andreoni
# email: igor.andreoni@gmail.com

import os
import requests

# Copy the API token from your Fritz account
token = 'copied-from-fritz'
token = os.environ.get('FRITZ_TOKEN')

def api(method, endpoint, data=None):
    """API request"""

    headers = {'Authorization': f'token {token}'}
    response = requests.request(method, endpoint, json=data, headers=headers)
    return response

def getFilter(filter_id):
    """GET a filter"""
    response = api('GET',
                   f'https://fritz.science/api/filters/{filter_id}'
                   )
    if response.status_code == 200:
        filt = response.json()['data']
    else:
        print(f'HTTP code: {response.status_code}, {response.reason}')

    return filt


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Get filters from Fritz.')
    parser.add_argument('filter_id', nargs='+', help='Filter IDs, space \
separated')
    parser.add_argument('--group', '-g', dest='group_id', default=1544,
                        type=int, help='Group ID')

    args = parser.parse_args()

    # Fix the filter IDs
    filter_ids = [int(x) for x in args.filter_id]

    for filter_id in filter_ids:
        filt = getFilter(filter_id)
        print(filt)

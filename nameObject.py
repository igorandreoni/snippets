# Written by Igor Andreoni (igor.andreoni@gmail.com)

import numpy as np
from astropy.time import Time

def getNext(names, base, mjd):
    """
    Find the next transient name

    Parameters
    ----------
    names list
        list of names of the last transients with the same base
        e.g. ["AT2018gfo", "AT2022cmc"]
    base str
        base for the name, e.g. "AT"
    mjd float
        MJD of the discovery

    Returns
    ------- 
    name_new str
        new transient name
    """
    # find the year of discovery
    yy = Time(mjd, format='mjd').iso[2:4] 
    # select from the list sources with same base and same year
    others = [n for n in names if (n[0:len(base)] == base) and
              n[len(base): len(base)+2] == yy]
    # select the latest
    others = [n for n in others if len(n) == max([len(o) for o in others])]
    if len(others) == 0:
        name_new = f"{base}{yy}a"
    else:
        others.sort()
        name_last = others[-1]
        # get the suffix
        sfx = list(name_last[len(base)+2:])
        # iterate from the end
        for i in np.arange(len(sfx))+1:
            # check if the last letter is not a z
            if sfx[-1*i] != "z":
                for idx in np.arange(i)+1:
                    sfx[-1*idx] = next_alpha(sfx[-1*idx]).lower()
                name_new = f"{base}{yy}{''.join(sfx)}"
                break
            # continue the iteration until there are no more items
            elif i == len(sfx):
                name_new = f"{base}{yy}{'a'*(len(sfx)+1)}"

    return name_new


def next_alpha(s):
    return chr((ord(s.upper())+1 - 65) % 26 + 65)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Get the next transient name')
    parser.add_argument('-n','--name', nargs='+', dest='names',
                        help='<Required> List of transient names \
e.g. AT2018gfo AT2022cmc',
                        default=None, required=True)
    parser.add_argument('--mjd', dest='mjd', type=float,
                        required=False,
                        help='MJD of the new transient detection \
(default: todays MJD)',
                        default=Time.now().mjd)
    parser.add_argument('-b', '--base', dest='base', type=str,
                        required=False,
                        help='Base name for the transients, e.g. AT',
                        default="GWM")
    args = parser.parse_args()

    #names = ["GWM23a", "GWM22c", "GWM25zz", "GWM23b", "GWM23ac", "GWM23bc"]
    new = getNext(args.names, args.base, args.mjd)
    print(new)

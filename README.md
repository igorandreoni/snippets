# snippets
Collection of useful astronomy snippets

### Get galaxies
usage: get_galaxies.py [-h] --ra RA --dec DEC [--r RAD] [--dist-min DIST_MIN]<br>
                       [--dist-max DIST_MAX] [--sep-max SEP_MAX_KPC]<br>
                       [--out OUT]<br>
<br>
Query the GLADE catalog (Dalya et al., 2018) to find galaxies around given coordinates<br>
<br>
optional arguments:<br>
  -h, --help            show this help message and exit<br>
  --ra RA               Right Ascension of the center of the query (deg)<br>
  --dec DEC             Declination of the center of the query (deg)<br>
  --r RAD               Radius of the query (deg), default=1deg<br>
  --dist-min DIST_MIN   Minimum distance (Mpc)<br>
  --dist-max DIST_MAX   Maximum distance (Mpc)<br>
  --sep-max SEP_MAX_KPC<br>
                        Maximum projected separation (kpc)<br>
  --out OUT             Output file name (CSV) default: galaxies.csv<br>

Example:

```
python get_galaxies.py --ra 219.51950000 --dec -60.34800000 --r 3. --out my_galaxies.csv --dist-max 100
```

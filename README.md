Collection of useful astronomy snippets

## Crop and display DBSP spectra
usage: dbsp_crop_spec.py [-h] [--doPlot] [-n NAMES [NAMES ...]] <br>
                         [--suffix OUT_SUFFIX]<br>
<br>
Crop DBSP spectra, prepare them for upload on Fritz, and plot them up. If no<br>
option is given, all spectra in the format ./ZTF\*.fits will be processed<br>
<br>
optional arguments:<br>
  -h, --help            show this help message and exit<br>
  --doPlot              Plot up the spectra<br>
  -n NAMES [NAMES ...]<br>
                        Names of the spectra; if not provided, all spectra in <br>
                        the format ./ZTF\*.fits will be processed <br>
  --suffix OUT_SUFFIX   suffix for the output; default = '_crop.txt' <br>

**Example:**

```
python dbsp_crop_spec.py --n ZTF19cool.fits ZTF20lesscool.fits --doPlot
```

**Limitations:** Headers are not yet transferred


## Upload spectra to Fritz

usage: upload_spectra_fritz.py [-h] [-d] [--date DATE] [--inst INST_ID]
                               ID [ID ...] <br>

Upload spectra to the Fritz marshal, individually or in bulk. <br>

positional arguments:<br>
  ID              Spectra filenames <br>

optional arguments:<br>
  -h, --help      show this help message and exit<br>
  -d              Use default reducer ID and observer - please customize the
                  code<br>
  --date DATE     Date of the observations, e.g. 2020-11-10T00:00:00<br>
  --inst INST_ID  Instrument ID, e.g. inst_id = 3 for DBSP. Instrument IDs can<br>
                  be found here: https://fritz.science/api/instrument<br>

**Notes:**
* If the name of the ZTF source is in the filename, it will be recognized automatically; filenames without the ZTF source name are fine, too
* Make sure that 'ZTF' appears only once in the filename, if any
* The spectra are assumed to have wavelength in the first column (in angstrom), flux in the second column, possibly flux error in the third column
* Please update your API before using the code. You can find the API code in your Fritz account page

**Useful links:**
* Instrument IDs: https://fritz.science/api/instrument
* User IDs: https://fritz.science/api/user
* Group IDs: https://fritz.science/api/groups
* Skyportal API guide: https://skyportal.io/docs/api.html
* Fritz marshal user's guide https://fritz-marshal.org/doc/

**Example:**
```
python upload_spectra_fritz.py ZTF*txt --inst 3 --date 2020-10-28T07:00:00
```

## Get galaxies
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

**Example:**

```
python get_galaxies.py --ra 219.51950000 --dec -60.34800000 --r 3. --out my_galaxies.csv --dist-max 100
```

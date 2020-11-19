# snippets
Collection of useful astronomy snippets

### Crop and display DBSP spectra
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

Example:

```
python dbsp_crop_spec.py --n ZTF19cool.fits ZTF20lesscool.fits --suffix _cropped.txt --doPlot
```

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

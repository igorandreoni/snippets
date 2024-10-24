from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pdb
from astropy.time import Time
from astropy.table import Table
from astropy.stats import sigma_clip

def doClean(t):
    # Clean images with null forcedimmimflux
    idx = [i for i in np.arange(len(t)) if t["forcediffimflux"][i] != 'null']
    if len(idx) < len(t):
        print(f"WARNING: {len(t)-len(idx)} points removed for null forcediffimflux")
    t = t[idx]
    # Convert to floats
    t["forcediffimflux"] = [float(x) for x in t["forcediffimflux"]]
    t["forcediffimfluxunc"] = [float(x) for x in t["forcediffimfluxunc"]]

    # Find outliers in scisigpix distribution for non-photometric nights
    t['scisigpix']
    filtered = sigma_clip(np.array(t['scisigpix']), sigma=2, maxiters=5)
    print(f"{np.size([x for x in filtered.mask if x==False])}/{len(t)} \
points were selected from scisigpix the distribution cut")

    # Points to be kept in the table
    idx = [i for i in np.arange(len(filtered.mask)) if filtered.mask[i]==False]
    t = t[idx]

    return t


def stack_lc(tbl, days_stack=1., snt_det=3, snt_ul=5):
    """Given a dataframe with a maxlike light curve,
    stack the flux """

    if 'jdobs' in list(tbl.colnames):
        key_jd = 'jdobs'
    elif 'jd' in list(tbl.colnames):
        key_jd = 'jd'
    else:
        print("What is the column for the JD??")
        pdb.set_trace()
    t_out = Table([[],[],[],[],[],[],[],[],[],[],[]],
                  names=(key_jd, 'forcediffimflux', 'forcediffimfluxunc', 'zpdiff', 'ezp',
                         'mag', 'mag_unc', 'diffmaglim', 'filter', 'programid','exptime'),
                  dtype=('double', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'S', 'f','f'))
    # Bin separately by filter
    filters = list(set(tbl['filter']))
    for f in filters:

        t = tbl[tbl['filter'] == f]
        bins = np.arange(int(np.max(t[key_jd]) - np.min(t[key_jd]))+2)
        dt0 = np.min(t[key_jd]) - int(np.min(t[key_jd]))
        if dt0 <= 0.4:
            start = int(np.min(t[key_jd])) - 0.6
        else:
            start = int(np.min(t[key_jd])) + 0.4
        bins = bins + start
        for b in bins:
            temp = t[(t[key_jd] > b) & (t[key_jd] < b+1)]
            # Remove rows of null flux
            temp = temp[[i for i in np.arange(len(temp)) if temp['forcediffimflux'][i] != "null"]]
            if len(temp) == 0:
                continue
            new_jd = np.mean(np.array(temp[key_jd]))
            new_exptime = np.sum(np.array(temp['exptime']))

            if len(set(temp['zpdiff'])) == 1:
                zp = temp['zpdiff'][0]
                #new_flux = np.mean(np.array(temp['forcediffimflux']))
                flux = np.array(temp['forcediffimflux'])
                # Make sure we are working with floats
                flux = np.array([float(fl) for fl in flux])
                nan_index = [i for i in np.arange(len(flux)) if np.isnan(flux[i])]
                for i in nan_index:
                    flux[i] = 0
                flux_unc = np.array(temp['forcediffimfluxunc'])
                flux_unc = np.array([float(fl) for fl in flux_unc])
                # Use weights only if there are only detections
                if np.min(flux/flux_unc) >= snt_det:
                    weights = (flux/flux_unc)**2
                    new_flux = np.sum(np.array([float(x) for x in temp['forcediffimflux']])*weights)/np.sum(weights)
                else:
                    new_flux = np.mean(np.array([float(ff) for ff in temp['forcediffimflux']]))
                new_flux_unc = np.sqrt(np.sum(np.array([float(ff) for ff in temp['forcediffimfluxunc']])**2))/len(temp)
            else:
                zp = temp['zpdiff'][0]
                flux1 = np.array(temp['forcediffimflux'])
                flux1_aux = []
                for i in np.arange(len(flux1)):
                   try:
                       flux1_aux.append(float(flux1[i]))
                   except:
                       flux1_aux.append(0)
                flux1 = np.array(flux1_aux)
                nan_index = [i for i in np.arange(len(flux1)) if np.isnan(flux1[i])]
                for i in nan_index:
                    flux1[i] = 0
                flux1_unc = np.array(temp['forcediffimfluxunc'])
                flux1_unc_aux = []
                for i in np.arange(len(flux1_unc)):
                   try:
                       flux1_unc_aux.append(float(flux1_unc[i]))
                   except:
                       flux1_unc_aux.append(0)
                flux1_unc = np.array(flux1_unc_aux)
                zp1 = np.array(temp['zpdiff'])
                flux = 10**((2.5*np.log10(flux1) - zp1 + zp ) / 2.5)
                flux_unc = 10**((2.5*np.log10(flux1_unc) - zp1 + zp ) / 2.5)
                flux[np.isnan(flux)] = 0
                # Use weights only if there are only detections
                if np.min(flux/flux_unc) >= snt_det:
                    weights = (flux/flux_unc)**2
                    new_flux = np.sum(flux*weights)/np.sum(weights)
                else:
                    new_flux = np.mean(flux)
                new_flux_unc = np.sqrt(np.sum(flux_unc**2))/len(temp)
            if new_flux/new_flux_unc > snt_det:
                mag_stack = -2.5*np.log10(new_flux) + zp
                mag_unc_stack = np.abs(-2.5*np.log10(new_flux-new_flux_unc) + 2.5*np.log10(new_flux))
                maglim_stack = 99.
            else:
                mag_stack = 99.
                mag_unc_stack = 99.
                maglim_stack = -2.5 * np.log10(snt_ul * new_flux_unc) + zp
            ezp = 0
            #ezp = np.sum(temp['zpdiffunc']**2)/len(temp)
            t_out.add_row([new_jd, new_flux, new_flux_unc, zp, ezp, mag_stack, mag_unc_stack, maglim_stack, f, 3, new_exptime])
    t_out[t_out['mag'] < 50].write(f"photometry_binned_snt{int(snt_det)}.csv",
                                   format='csv', overwrite=True)
    t_out_det = t_out[t_out["mag"] < 60]
    t_out_det.sort("jd")
    for i in np.arange(len(t_out_det)):
        print(f"{t_out_det['jd'][i]:.5f} & {t_out_det['mag'][i]:.3f} & {t_out_det['mag_unc'][i]:.3f} & {t_out_det['filter'][i]} & {t_out_det['exptime'][i]} \\\\")
    return t_out



def main(args):
    try:
        t = ascii.read(args.input_filename, header_start=0, data_start=1)
    except:
        t = ascii.read(args.input_filename, header_start=0, data_start=200)
    #Remove commas from the column names
    for n in t.colnames:
        t[n].name = t[n].name[:-1]
    #Set t0
    t0 = float(Time.now().jd)
    #t0 = 2458710.3823958
    time_arr=[]

    # Clean the data using image quality metrics
    if args.doClean is True:
        t = doClean(t)

    # Select only data for given program IDs
    if args.pid is not None:
        pids = [int(p) for p in args.pid.split(",")]
        idx = [i for i in np.arange(len(t)) if t["programid"][i] in pids]
        t = t[idx]
        t.write(args.input_filename.replace(".txt", f"_pid{args.pid}.txt"), overwrite=True, format='csv')
    filter_color={'ZTF_g':'darkturquoise', 'ZTF_r':'r', 'ZTF_i':'y', 'c': 'cyan', 'o': 'orange'}
    t_small = Table([t["jd"],t['filter']],names=("JD", "filter"))


    # Stack the light curve
    if args.doStack is True:
        stacked = stack_lc(t, snt_det=5)
        t = stacked

    #for l in stacked:
    #    print(f"{Time(l['jd'], format='jd').iso} | {l['filter'][-1]} | >{l['limmag']} <br>")

    #import pdb
    #pdb.set_trace()

    plt.figure(figsize=(9,6))
    for j in t['jd']:
        time_arr.append(j-t0)
    t['jd_modified'] = time_arr

    for f in ['ZTF_g', 'ZTF_r', 'ZTF_i']:
        time=t['jd_modified'][t['filter']==f]
        
        flux_arr = t['forcediffimflux'][t['filter'] == f]	
        flux_err_arr = t['forcediffimfluxunc'][t['filter'] == f]
        zpdiff_arr=t['zpdiff'][t['filter'] == f]
        sci_maglim_arr = t['diffmaglim'][t['filter'] ==f]
        pid_arr = t['programid'][t['filter']==f]
        mag_arr=[]
        mag_err_arr=[]
        ul_arr=[]
        exptime_arr = t['exptime'][t['filter']==f]
        SNT = 5.
        SNU = 5.
        for flux, flux_err, zpdiff, sci_maglim in zip(flux_arr, flux_err_arr, zpdiff_arr, sci_maglim_arr):
#            if f == 'ZTF_i':
#                pdb.set_trace()
            try:
                flux, flux_err, zpdiff = float(flux), float(flux_err), float(zpdiff)
            except ValueError:
                flux, flux_err = 1., 1.
                try:
                    zpdiff = float(zpdiff)
                except ValueError:
                    pdb.se_trace()
            try:
                x = (flux, flux_err, flux / flux_err)
            except TypeError:
                pdb.set_trace()
            if (flux / flux_err) > SNT:
                # we have a 'confident' detection, compute and plot mag with error bar:
                mag = zpdiff - 2.5*np.log10(flux)
                mag_err = 1.0857* flux_err / flux 
                mag_arr.append(mag)
                mag_err_arr.append(mag_err)
                ul_arr.append(0)
            else:
            # compute flux upper limit and plot as arrow:mag = zpdiff - 2.5*log10[SNU*forcediffimfluxunc]
                #ul_arr.append(zpdiff - 2.5*np.log10(SNU*flux_err)) 
                ul_arr.append(sci_maglim)
                mag_arr.append(99)
                mag_err_arr.append(0) 

               			
        #PLOT!
        plt.errorbar(time, mag_arr, yerr=mag_err_arr, color=filter_color[f], marker='o', markeredgecolor='k', linestyle="", label=f"{f}")
        #plt.plot(time, mag_arr, filter_color[f]+'o',color='black')
        plt.plot(time, ul_arr, color=filter_color[f], marker="v", linestyle="")

        print("############")
        print(f)
        count = 0
        print("time | mag | programid | exptime")
        for tt, mm, mmerr, ull, pid, exptime in zip(time, mag_arr, mag_err_arr, ul_arr, pid_arr, exptime_arr):
            time_iso = Time(tt+t0, format='jd').iso
            if mm > 5. and mm < 30:
                print(f"{time_iso} | {f[-1]} = {mm:.2f} +- {mmerr:.2f} | {pid} | {exptime}")
                count += 1
            else:
                print(f"{time_iso} | {f[-1]} > {ull:.1f} | {pid}")
        print(f"Total of {count} detections in {f} filter with S/N>{SNT}")
        print("-----")
        #print(mag_arr, mag_err_arr)

    # Plot ATLAS
    if args.atlas is not None:
        atlas = ascii.read(args.atlas)
        atlas.rename_column("##MJD", "MJD")
        for f in set(atlas['F']):
            atlasf = atlas[atlas['F'] == f]
            atlasf_det = atlasf[atlasf['dm'] <= 0.3]
            time_det = np.array(Time(atlasf_det["MJD"], format='mjd').jd - t0)
            atlasf_ndet = atlasf[atlasf['dm'] > 0.3]
            time_ndet = np.array(Time(atlasf_ndet["MJD"], format='mjd').jd - t0)
            plt.errorbar(time_det, atlasf_det['m'], yerr=atlasf_det['dm'], color=filter_color[f], marker='o', markeredgecolor='k', linestyle="", label=f"{f} ATLAS")
            plt.plot(time_ndet, atlasf_ndet['mag5sig'], color=filter_color[f], marker="v", linestyle="")

    plt.xlabel('Days from '+str(Time.now().iso), fontsize=19)
    plt.ylabel('magnitude [AB]', fontsize=19)
    #plt.axis([xmin,xmax,ymin,ymax])
    plt.tick_params(labelsize=18)
    plt.gca().invert_yaxis()
    plt.legend(fontsize = 18)

    print(f"Total observations: {len(t)}")
    print(f"Total observations g: {len(t[t['filter'] == 'ZTF_g'])}")
    print(f"Total observations r: {len(t[t['filter'] == 'ZTF_r'])}")
    print(f"Total observations i: {len(t[t['filter'] == 'ZTF_i'])}")
    print(f"Total nights: {len(set([np.round(jd) for jd in t['jd']]))}")

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot light curves \
            obtained with ZTF forced photometry service (Masci et al., 2019')

    parser.add_argument('-i', dest='input_filename', type=str, required=True,
    help='Light curve filename')
    parser.add_argument('-a', dest='atlas', type=str, required=False,
    help='ATLAS light curve')
    parser.add_argument('-pid', dest='pid', type=str, required=False,
    help='program IDs, example 1,2', default=None)
    parser.add_argument('--doClean', dest='doClean', action='store_true', required=False,
    help='Clean the data using image quality metrics', default=False)
    parser.add_argument('--doStack', dest='doStack', action='store_true', required=False,
    help='Stack nightly', default=False)
#    parser.add_argument('-s', dest='saveornot', type=str2bool, required=False, \
#    help='Save the plot (PNG format)')
#    parser.add_argument('-o', dest='output_filename', type=str, required=False, default=' ', \
#    help='Output light curve file name if saved (-s) ')
#    parser.add_argument('-t', dest='plot_title', type=str, required=False, default=' ', \
#    help='Plot title')
#    parser.add_argument('-yl', dest='y_label', type=str, required=False, default='g mag', \
#    help='Y axis label')
	    
    args = parser.parse_args()
    main(args)

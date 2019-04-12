'''
Some misc functions.
'''
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def load_data_pd(fn, tp=''):
    '''Load data file.'''
    print('>> Loading data: {}'.format(fn))
    tb = pd.read_table(fn, delim_whitespace=True, comment='#', header=None)
    tb = tb.values
    if tp == 'knife':
        # RA, DEC, weight
        return np.column_stack((tb[:, 0], tb[:, 1], tb[:, 3]))
    else:
        return tb


def save_jk_map(jk_map, fn):
    '''Save jackknife map to fits file.'''
    hp.write_map(fn, jk_map, overwrite=True)
    print(':: Jackknife map saved to file: {}'.format(fn))


def plot_jk_map(jk_map, shuffle=False, njr=0):
    '''Plot jackknife map.'''
    print(':: Plotting jackknife regions in Healpix map')
    if shuffle:
        print('-- shuffle the labels, looks better, demo only')
        assert (njr != 0), 'zero njr'
        lb_max = np.int(np.amax(jk_map))
        lb_min = lb_max - njr + 1
        arr = np.array([i for i in range(lb_min, lb_max+1, 1)])
        np.random.shuffle(arr)
        for i, p in enumerate(jk_map):
            if p != hp.UNSEEN:
                jk_map[i] = arr[np.int(p)-lb_min]

    hp.mollview(jk_map, coord='GC', title='Jackknife regions')
    plt.show()


def save_jk_bounds(jk_bounds, fn):
    '''Save jackknife bounds to txt file.'''
    header = 'Number of jackknife regions: {0:d}\n'.format(len(jk_bounds))
    header += 'RA_min   RA_max   DEC_min   DEC_max'
    np.savetxt(fn, jk_bounds, header=header)
    print(':: Jackknife bounds saved to file: {}'.format(fn))


def combine_bounds(bdfs, fn):
    '''Combine a few bounds files into one.'''
    for i, f in enumerate(bdfs):
        if i == 0:
            data = np.loadtxt(f)
            continue
        tmp = np.loadtxt(f)
        data = np.row_stack((data, tmp))
    header = 'Number of jackknife regions: {0:d}\n'.format(len(data))
    header += 'RA_min   RA_max   DEC_min   DEC_max'
    np.savetxt(fn, data, header=header)
    print(':: Combined jackknife bounds saved to file: {}'.format(fn))


def rand2map(rand, nside):
    '''Make Healpix map with labeled random points.'''
    npix = hp.nside2npix(nside)
    jk_map = np.full(npix, hp.UNSEEN)

    rot = hp.Rotator(coord=['C', 'G'])
    # convert to theta, phi used by default in Healpy
    theta_equ, phi_equ = np.deg2rad(90.-rand[:, 1]), np.deg2rad(rand[:, 0])
    theta_gal, phi_gal = rot(theta_equ, phi_equ)
    ipix = hp.ang2pix(nside, theta_gal, phi_gal)

    jk_map[ipix] = rand[:, 4]

    return jk_map


def analyze_rand(rand, sf):
    '''Analyze the labeled random points.
    Columns in [0:RA, 1:DEC, 2:Z, 3:weight, 4:jackknife label].'''

    jk_min = np.int(np.amin(rand[:, 4]))
    jk_max = np.int(np.amax(rand[:, 4]))
    njr = np.int(jk_max - jk_min + 1)
    print('>> Number of jackknife regions: {0:d}'.format(njr))
    print('>> Labeled from {0:d} to {1:d}'.format(jk_min, jk_max))

    w_tot = np.sum(rand[:, 3])
    print('>> Total weight: {0:f}'.format(w_tot))
    w_ave = w_tot / njr
    print('>> Average weight per region: {0:f}'.format(w_ave))

    print('>> Analyzing weights for each region')
    w_jk = np.zeros(njr)
    for p in rand:
        w_jk[np.int(p[4])-jk_min] += p[3]

    print('>> Percent deviation from average for each region')
    pcdev = 100. * (w_jk - w_ave) / w_ave
    print(pcdev)

    jk_map = rand2map(rand, 256)
    plot_jk_map(jk_map, shuffle=sf, njr=njr)


def merge_masks(fmasks, nside, fo):
    '''Merge masks each with jk labels start from 0.'''
    npix = hp.nside2npix(nside)
    mask_tot = np.full(npix, hp.UNSEEN)
    njr = 0
    for fn in fmasks:
        print('>> Loading mask: {}'.format(fn))
        mask = hp.read_map(fn)
        mask_tot = np.where(mask != hp.UNSEEN, mask + njr, mask_tot)
        njr = np.amax(mask) + 1

    hp.write_map(fo, mask_tot)
    print('>> Merged mask written to file: {}'.format(fo))


def get_ra_dec(theta, phi):
    '''Get RA, DEC [degree] from theta, phi [radians].'''
    rot = hp.Rotator(coord=['G', 'C'])
    theta_equ, phi_equ = rot(theta, phi)
    dec, ra = 90. - np.rad2deg(theta_equ), np.rad2deg(phi_equ)
    # move RA in [-180,0) to [180,360)
    ra = np.where(ra < 0., ra + 360., ra)

    return ra, dec


def sep_mask_ra(fmask, RAs, fo):
    '''Seperate one mask into two in RA direction. [0,360]'''
    print('>> Loading mask: {}'.format(fmask))
    mask = hp.read_map(fmask)
    nside = hp.get_nside(mask)
    npix = hp.nside2npix(nside)
    ipix = np.array([i for i in range(npix)])

    us = np.min(mask)

    theta, phi = hp.pix2ang(nside, ipix)
    ra, dec = get_ra_dec(theta, phi)
    del dec

    mask1 = np.where((RAs[0] < ra) & (ra < RAs[1]), mask, us)
    mask2 = np.where((RAs[0] > ra) | (ra > RAs[1]), mask, us)

    hp.write_map(fo[0], mask1, overwrite=True)
    print('>> Mask 1 written to: {}'.format(fo[0]))
    print('>> Mask 1 total weights: {0:f}'.format(np.sum(mask1)))

    hp.write_map(fo[1], mask2, overwrite=True)
    print('>> Mask 2 written to: {}'.format(fo[1]))
    print('>> Mask 2 total weights: {0:f}'.format(np.sum(mask2)))

'''
Some useful functions.
'''
import pandas as pd
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt


#--- General ---#


def load_data_pd(fn, tp='', verbose=True):
    '''Load data file.'''
    if verbose:
        print('>> Loading data: {}'.format(fn))
    df = pd.read_csv(fn, delim_whitespace=True, comment='#', header=None)
    df = df.to_numpy()
    if tp == 'knife':
        # RA, DEC, weight
        return np.column_stack((df[:, 0], df[:, 1], df[:, 3]))
    else:
        return df


def get_ra_dec(theta, phi):
    '''Get RA, DEC [degree] from theta, phi [radians] used in Healpy.'''
    rot = hp.Rotator(coord=['G', 'C'])
    theta_equ, phi_equ = rot(theta, phi)
    dec, ra = 90. - np.rad2deg(theta_equ), np.rad2deg(phi_equ)
    # move RA in [-180,0) to [180,360)
    ra = np.where(ra < 0., ra + 360., ra)

    return ra, dec


def get_theta_phi(ra, dec):
    '''Get theta, phi [radians] used in Healpy from RA, DEC [degree].'''
    rot = hp.Rotator(coord=['C', 'G'])
    theta_equ, phi_equ = np.deg2rad(90.-dec), np.deg2rad(ra)
    theta_gal, phi_gal = rot(theta_equ, phi_equ)
    # move phi in [-pi,0) to [pi,2pi)
    phi_gal = np.where(phi_gal < 0., phi_gal + 2*np.pi, phi_gal)

    return theta_gal, phi_gal


#--- Map ---#


def save_jk_map(jk_map, fn):
    '''Save jackknife map to fits file.'''
    hp.write_map(fn, jk_map, overwrite=True)
    print(':: Jackknife map saved to file: {}'.format(fn))


def plot_jk_map(jk_map, shuffle=False, njr=0, cmap=None):
    '''Plot jackknife map.'''
    print(':: Plotting jackknife regions in Healpix map')
    if shuffle:
        print('-- shuffle the labels, looks better, demo only')
        lb_max = np.int(np.amax(jk_map))
        lb_min = lb_max - njr + 1
        arr = np.array([i for i in range(lb_min, lb_max+1, 1)])
        np.random.shuffle(arr)
        for i, p in enumerate(jk_map):
            if p != hp.UNSEEN:
                jk_map[i] = arr[np.int(p)-lb_min]
    if cmap is None:
        hp.mollview(jk_map, coord='GC', title='jackknife regions')
    else:
        hp.mollview(jk_map, coord='GC', cmap=cmap, title='jackknife regions')
    plt.show()


def merge_jk_maps(maps, fo=None):
    '''Merge jk maps each with label start from 0.'''
    npix = len(maps[0])
    map_tot = np.full(npix, hp.UNSEEN)
    njr = 0
    for i in range(len(maps)):
        map_ = maps[i]
        map_tot = np.where(map_ != hp.UNSEEN, map_ + njr, map_tot)
        njr += np.amax(map_) + 1

    if fo is not None:
        hp.write_map(fo, map_tot)
        print('>> Merged jackknife map written to file: {}'.format(fo))

    return map_tot


#--- Bounds ---#


def save_bounds(jk_bounds, fn):
    '''Save jackknife bounds to txt file.'''
    header = 'Number of jackknife regions: {0:d}\n'.format(len(jk_bounds))
    header += 'RA_min   RA_max   DEC_min   DEC_max'
    np.savetxt(fn, jk_bounds, header=header)
    print(':: Jackknife bounds saved to file: {}'.format(fn))


def combine_bounds(bds):
    '''Combine a few bounds into one.'''
    for i, bd in enumerate(bds):
        if i == 0:
            data = bd
            continue
        data = np.row_stack((data, bd))
    return data


def in_bound(ang, bd):
    '''Check if array ang(ra,dec) are in the bound.'''
    res = np.full(ang.shape[0], False, dtype=np.bool)
    # check DEC
    idx = np.where((bd[2] <= ang[:, 1]) & (ang[:, 1] <= bd[3]))[0]
    ang = ang[idx]
    # check RA
    if bd[0] < bd[1]:
        idx_ = np.where((bd[0] <= ang[:, 0]) & (ang[:, 0] <= bd[1]))[0]
    else:  # region crossing zero
        idx_ = np.where((bd[0] <= ang[:, 0]) | (ang[:, 0] <= bd[1]))[0]
    res[idx[idx_]] = True

    return res

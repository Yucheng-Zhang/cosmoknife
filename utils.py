'''
Some useful functions.
'''
import pandas as pd
import numpy as np
import healpy as hp


#--- General ---#


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

    return theta_gal, phi_gal


#--- Map ---#


def save_jk_map(jk_map, fn):
    '''Save jackknife map to fits file.'''
    hp.write_map(fn, jk_map, overwrite=True)
    print(':: Jackknife map saved to file: {}'.format(fn))


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

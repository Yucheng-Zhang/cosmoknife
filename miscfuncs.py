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


def plot_jk_map(jk_map):
    '''Plot jackknife map.'''
    print(':: Plotting jackknife regions in Healpix map')
    print(':::: note: this plot is labeled with [0, 1, ..., njr-1]')
    hp.mollview(jk_map, coord='GC')
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

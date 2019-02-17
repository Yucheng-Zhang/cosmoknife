'''
Make jackknife regions based on the weights of the random points.
'''
import healpy as hp
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make jackknife regions.')
    parser.add_argument('-njr', type=int, default=50,
                        help='Number of jackknife regions.')
    parser.add_argument('-nra', type=int, default=10,
                        help='Number of slices in RA.')
    parser.add_argument('-rand', type=str, default='',
                        help='Random file, with columns: RA, DEC, Z, weight')
    parser.add_argument('-nside', type=int, default=256,
                        help='Jackknife regions will be a Healpix map.')
    args = parser.parse_args()

    njr = args.njr
    n_ra = args.nra
    nside = args.nside


def load_data_pd(fn):
    '''Load random data file.'''
    print('>> Loading data: {}'.format(fn))
    tb = pd.read_table(fn, delim_whitespace=True, comment='#', header=None)
    rand = tb.values
    # RA, DEC, weight
    return np.column_stack((rand[:, 0], rand[:, 1], rand[:, 3]))


def cut_in_ra(rand, w_ra):
    '''Cut in the RA direction.'''
    print('>> Cutting in the RA direction')
    d_ra = {}
    rand = rand[rand[:, 0].argsort()]  # sort along RA
    tmp = 0.
    i0 = 0
    j = 0
    nps = len(rand[:, 0])
    for i1, p in enumerate(rand, 1):
        tmp += p[2]
        if tmp > w_ra or i1 == nps:  # or reach the end
            d_ra[j] = rand[i0:i1, :]
            j += 1
            tmp = 0.
            i0 = i1

    del rand
    return d_ra


def cut_in_dec(d_ra, w_dec):
    '''Cut in the DEC direction.'''
    print('>> Cutting in the DEC direction')
    d_dec = {}
    j = 0
    for i in range(len(d_ra)):
        rand = d_ra[i]  # points in the RA piece
        del d_ra[i]
        rand = rand[rand[:, 1].argsort()]  # sort each RA piece along DEC
        tmp = 0
        i0 = 0
        nps = len(rand[:, 0])
        for i1, p in enumerate(rand, 1):
            tmp += p[2]
            if tmp > w_dec or i1 == nps:
                d_dec[j] = rand[i0:i1, :]
                j += 1
                tmp = 0.
                i0 = i1

    del d_ra
    return d_dec


def make_jk_map(d_dec, nside):
    '''Make healpix map for the jackknife regions.'''
    print('>> Making healpix map for jackknife regions')
    npix = hp.nside2npix(nside)
    jk_map = np.zeros(npix)
    rot = hp.Rotator(coord=['C', 'G'])
    for i in range(len(d_dec)):
        rand = d_dec[i]
        del d_dec[i]
        # convert to theta, phi used by default in Healpy
        theta_equ, phi_equ = np.deg2rad(90.-rand[:, 1]), np.deg2rad(rand[:, 0])
        theta_gal, phi_gal = rot(theta_equ, phi_equ)
        ipix = hp.ang2pix(nside, theta_gal, phi_gal)
        jk_map[ipix] = i + 1.

    del d_dec
    return jk_map


def knife(rand, njr, n_ra, nside, plot=False):
    '''Main function.'''
    n_dec = int(njr/n_ra)

    w_total = np.sum(rand[:, 2])
    w_ra = w_total / n_ra
    w_dec = w_ra / n_dec

    d_ra = cut_in_ra(rand, w_ra)
    d_dec = cut_in_dec(d_ra, w_dec)
#    print('len(d_dec) == njr :', len(d_dec) == njr)

    jk_map = make_jk_map(d_dec, nside)
    if plot:
        hp.mollview(jk_map, coord='GC')
        plt.show()

    return jk_map


if __name__ == "__main__":

    rand = load_data_pd(args.rand)
    jk_map = knife(rand, njr, n_ra, nside, plot=True)

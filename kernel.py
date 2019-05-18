'''
Kernel functions.
'''
import numpy as np
import healpy as hp
from . import utils


def cut_in_ra(rand, w_ra, rra):
    '''Cut in the RA direction. RA in [0, 360].'''
    print('>> Cutting in the RA direction')
    d_ra = {}

    if rra != 0.:  # rotate rra if cross 0
        print('++ Rotate RA for {0:f} degrees'.format(rra))
        tmp = (rand[:, 0] + rra) % 360.
    else:
        tmp = rand[:, 0]

    rand = np.column_stack((rand, tmp))  # helping column for RA

    rand = rand[rand[:, 3].argsort()]  # sort along RA
    tmp, i0, j = 0., 0, 0
    nps = len(rand[:, 0])

    for i1, p in enumerate(rand, 1):
        tmp += p[2]
        if tmp > w_ra[j] or i1 == nps:  # or reach the end
            d_ra[j] = rand[i0:i1, :]
            j += 1
            tmp = 0.
            i0 = i1

    del rand
    return d_ra


def cut_in_dec(d_ra, w_dec):
    '''Cut in the DEC direction. DEC in [-90, 90].'''
    print('>> Cutting in the DEC direction')
    d_dec = {}
    j = 0
    for i in range(len(d_ra)):
        rand = d_ra[i]  # points in the RA piece
        rand = rand[rand[:, 1].argsort()]  # sort each RA piece along DEC
        tmp, i0 = 0, 0
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


def knife(data, njr, nra, rra):
    '''Knife function. data includes 3 columns [ra, dec, weight].'''
    w_total = np.sum(data[:, 2])
    w_dec = w_total / njr  # weight for final jackknife regions

    # make weights for RA pieces
    if njr % nra == 0:
        w_ra = np.full(nra, w_total/nra)
    else:
        res = njr % (nra-1)
        ndec = (njr - res) / (nra - 1)
        w_ra = np.full(nra - 1, ndec * w_dec)
        w_ra = np.append(w_ra, res * w_dec)

    d_ra = cut_in_ra(data, w_ra, rra)
    d_dec = cut_in_dec(d_ra, w_dec)

    return d_dec


def make_jk_bounds(d_dec):
    '''Make bounds [ra_min, ra_max, dec_min, dec_max] for jackknife regions.'''
    print('>> Making RA, DEC bounds for jackknife regions')
    jk_bounds = np.zeros((len(d_dec), 4))

    for i in range(len(d_dec)):
        rand = d_dec[i]
        k = np.argmin(rand[:, 3])
        jk_bounds[i, 0] = rand[k, 0]  # RA min
        k = np.argmax(rand[:, 3])
        jk_bounds[i, 1] = rand[k, 0]  # RA max
        jk_bounds[i, 2] = np.amin(rand[:, 1])  # DEC min
        jk_bounds[i, 3] = np.amax(rand[:, 1])  # DEC max

    return jk_bounds


def make_jk_map(d_dec, nside):
    '''Make healpix map for the jackknife regions.'''
    print('>> Making healpix map for jackknife regions')
    npix = hp.nside2npix(nside)
    jk_map = np.full(npix, hp.UNSEEN)

    for i in range(len(d_dec)):
        rand = d_dec[i]
        # convert to theta, phi used by default in Healpy
        theta_gal, phi_gal = utils.get_theta_phi(rand[:, 0], rand[:1])
        ipix = hp.ang2pix(nside, theta_gal, phi_gal)
        jk_map[ipix] = i

    return jk_map

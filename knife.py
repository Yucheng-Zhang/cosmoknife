'''
Make jackknife regions based on the weights of the random points.
'''
import healpy as hp
import numpy as np


def cut_in_ra(rand, w_ra, rra):
    '''Cut in the RA direction., RA in [0, 360].'''
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
        if tmp > w_ra or i1 == nps:  # or reach the end
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


def make_jk_map(d_dec, nside):
    '''Make healpix map for the jackknife regions.'''
    print('>> Making healpix map for jackknife regions')
    npix = hp.nside2npix(nside)
    jk_map = np.full(npix, hp.UNSEEN)
    rot = hp.Rotator(coord=['C', 'G'])
    for i in range(len(d_dec)):
        rand = d_dec[i]
        # convert to theta, phi used by default in Healpy
        theta_equ, phi_equ = np.deg2rad(90.-rand[:, 1]), np.deg2rad(rand[:, 0])
        theta_gal, phi_gal = rot(theta_equ, phi_equ)
        ipix = hp.ang2pix(nside, theta_gal, phi_gal)
        jk_map[ipix] = i

    return jk_map


def make_jk_bounds(d_dec, rra):
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


def knife(rand, njr, n_ra, nside, rra):
    '''Main function.'''
    n_dec = int(njr/n_ra)

    w_total = np.sum(rand[:, 2])
    w_ra = w_total / n_ra
    w_dec = w_ra / n_dec

    d_ra = cut_in_ra(rand, w_ra, rra)
    d_dec = cut_in_dec(d_ra, w_dec)

    jk_map = make_jk_map(d_dec, nside)
    jk_bounds = make_jk_bounds(d_dec, rra)

    return jk_map, jk_bounds

'''
Make jackknife regions based on the weights of the random points.
'''
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


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


def save_jk_map(jk_map, fn):
    '''Save jackknife map to fits file.'''
    hp.write_map(fn, jk_map, overwrite=True)
    print(':: Jackknife map saved to file: {}'.format(fn))


def make_jk_map(d_dec, nside):
    '''Make healpix map for the jackknife regions.'''
    print('>> Making healpix map for jackknife regions')
    npix = hp.nside2npix(nside)
    jk_map = np.zeros(npix)
    rot = hp.Rotator(coord=['C', 'G'])
    for i in range(len(d_dec)):
        rand = d_dec[i]
        # convert to theta, phi used by default in Healpy
        theta_equ, phi_equ = np.deg2rad(90.-rand[:, 1]), np.deg2rad(rand[:, 0])
        theta_gal, phi_gal = rot(theta_equ, phi_equ)
        ipix = hp.ang2pix(nside, theta_gal, phi_gal)
        jk_map[ipix] = i + 1.

    return jk_map


def save_jk_bounds(jk_bounds, fn):
    '''Save jackknife bounds to txt file.'''
    header = 'Number of jackknife regions: {0:d}\n'.format(len(jk_bounds))
    header += 'RA_min   RA_max   DEC_min   DEC_max'
    np.savetxt(fn, jk_bounds, header=header)
    print(':: Jackknife bounds saved to file: {}'.format(fn))


def make_jk_bounds(d_dec):
    '''Make bounds [ra_min, ra_max, dec_min, dec_max] for jackknife regions.'''
    print('>> Making RA, DEC bounds for jackknife regions')
    jk_bounds = np.zeros((len(d_dec), 4))
    for i in range(len(d_dec)):
        rand = d_dec[i]
        jk_bounds[i, 0] = np.amin(rand[:, 0])  # RA min
        jk_bounds[i, 1] = np.amax(rand[:, 0])  # RA max
        jk_bounds[i, 2] = np.amin(rand[:, 1])  # DEC min
        jk_bounds[i, 3] = np.amax(rand[:, 1])  # DEC max

    return jk_bounds


def knife(rand, njr, n_ra, nside, fmap='', fbounds='', plot=False):
    '''Main function.'''
    n_dec = int(njr/n_ra)

    w_total = np.sum(rand[:, 2])
    w_ra = w_total / n_ra
    w_dec = w_ra / n_dec

    d_ra = cut_in_ra(rand, w_ra)
    d_dec = cut_in_dec(d_ra, w_dec)

    jk_map = make_jk_map(d_dec, nside)
    jk_bounds = make_jk_bounds(d_dec)

    if fmap != '':
        save_jk_map(jk_map, fmap)

    if fbounds != '':
        save_jk_bounds(jk_bounds, fbounds)

    if plot:
        print(':: Plotting jackknife regions in Healpix map')
        hp.mollview(jk_map, coord='GC')
        plt.show()

    return jk_map, jk_bounds

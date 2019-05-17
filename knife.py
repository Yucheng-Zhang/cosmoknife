'''
Make jackknife regions.
'''
import utils
import kernel
import numpy as np
import healpy as hp


def knife_rand(frand, njr, nra, rra, nside=None):
    '''Make jackknife regions with randoms.
       Columns: [RA, DEC, redshift, weight].'''
    rand = utils.load_data_pd(frand, tp='knife')

    d_dec = kernel.knife(rand, njr, nra, rra)

    jk_bounds = kernel.make_jk_bounds(d_dec)

    if nside is None:
        return jk_bounds
    else:
        jk_map = kernel.make_jk_map(d_dec, nside)
        return jk_bounds, jk_map


def knife_mask(fmask, njr, nra, rra, nside):
    '''Make jackknife regions with mask.
       Value on each pixel should be in range [0, 1] or hp.UNSEEN+(0,1].'''
    print('>> Loading mask: {}'.format(fmask))
    mask = hp.read_map(fmask)
    nside = hp.get_nside(mask)
    npix = hp.nside2npix(nside)
    ipix = np.array([i for i in range(npix)])

    # cut off the pixels with value 0 or UNSEEN
    idx = np.where((mask != 0.) & (mask != hp.UNSEEN))
    mask, ipix = mask[idx], ipix[idx]

    theta, phi = hp.pix2ang(nside, ipix)
    ra, dec = utils.get_ra_dec(theta, phi)
    data = np.column_stack((ra, dec, mask))

    d_dec = kernel.knife(data, njr, nra, rra)

    jk_map = kernel.make_jk_map(d_dec, nside)

    return jk_map

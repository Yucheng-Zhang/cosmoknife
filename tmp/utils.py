'''
Script to use some useful functions.
'''
import argparse
import miscfuncs as mf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Useful funcs for jackknife analysis.')

    parser.add_argument('-func', type=str, default='',
                        help='''Function to use.''')
    #--- 0: merge_masks ---#
    parser.add_argument('-fmasks0', type=str, nargs='+',
                        help='Masks to merge.')
    parser.add_argument('-nside0', type=int, help='nside of the masks.')
    parser.add_argument('-fo0', type=str, default='merged-mask.fits')
    #--- 1: sep_mask_ra ---#
    parser.add_argument('-fmask1', type=str, help='Mask to seperate.')
    parser.add_argument('-RAs1', type=float, nargs='+',
                        help='Two seperations in RA.')
    parser.add_argument('-fo1', type=str, nargs='+',
                        default=['sep-mask1.fits', 'sep-mask2.fits'])
    #--- 2: make_jk_masks ---#
    parser.add_argument('-fmask2', type=str)
    parser.add_argument('-fjkmask2', type=str)
    parser.add_argument('-froot2', type=str)

    args = parser.parse_args()


if __name__ == "__main__":

    if args.func == 'merge_masks':
        mf.merge_masks(args.fmasks0, args.nside0, args.fo0)
    elif args.func == 'sep_mask_ra':
        mf.sep_mask_ra(args.fmask1, args.RAs1, args.fo1)
    elif args.func == 'make_jk_masks':
        mf.make_jk_masks(args.fmask2, args.fjkmask2, args.froot2)
    else:
        pass

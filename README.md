# Jackknife

- Make jackknife regions based on the weights of the random points.
  - First sort and cut into pieces in the RA direction.
  - Then for each RA piece, sort and cut into sub-pieces in the DEC direction.

- Mark the jackknife regions.
  - Bounds in RA and DEC for each region: [RA_min, RA_max, DEC_min, DEC_max].
  - Healpix map may also be used to mark the angular jackknife regions.

- For now, input random and data file should have four columns: `RA, DEC, redshift, weight`.
- Output file have one more columns with jackknife regions label, starting from `0`.

## Usage
Example
```bash
#!/bin/bash

frand=rand_DR14_QSO_S_radec_z.dat # random file
fdata=data_DR14_QSO_S_radec_z.dat # data file

n_jk=20 # number of jackknife regions
n_ra=5 # number of slices in the RA direction
n_side=256 # nside for Healpix map of jk regions

f_map=out_jk_map.fits # map for jk regions
f_bounds=out_jk_bounds.dat # bounds for jk regions

fo_rand=rand_DR14_QSO_S_rdzwj.dat # output random file
fo_data=data_DR14_QSO_S_rdzwj.dat # output data file

mark=bounds # should be better than map

jk0=0 # Initial label of jackknife regions. All the labels will be [jk0, jk0+1, ..., jk0+njr-1].

python jack.py -rand $frand -data $fdata -forand $fo_rand -fodata $fo_data -njr $n_jk -nra $n_ra -fmap $f_map -fbounds $f_bounds -nside $n_side -tp $mark -jk0 $jk0

```

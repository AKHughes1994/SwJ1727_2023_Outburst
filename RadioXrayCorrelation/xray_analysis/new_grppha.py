#!/usr/bin/python3

import os, optparse, sys
import numpy as np
from astropy.io import fits

def parse_args( argv):

    desc="""%prog creates a grouped pha spectral file suitable for Xspec over a specified energy range. User must specify input and ouput file. User can specifiy minimum binning based on spectral resolution and minimum number of counts per bin."""
    
    parser = optparse.OptionParser(description=desc)

    parser.add_option('-i', '--input', \
                    help='input pi file', \
                    dest='input_pi')
    parser.add_option('-r', '--rmf', \
                        help='rmf file', \
                        dest='rmf')
    parser.add_option('-a', '--arf', \
                        help='arf file', \
                        dest='arf')
    parser.add_option('-b', '--background', \
                        help='background pi file', \
                        dest='background')
    parser.add_option('-o', '--output', \
                        help='output (grouped) pi file', \
                        dest='output_pi')
    parser.add_option('-e', '--energy_bin', \
                        help='minimum energy change per bin in eV or fraction [%default]', \
                        dest='bin_min_energy', \
                        type='float', \
                        default=100.0)
    parser.add_option('-f', '--fractional_energy_binning', \
                        help='switch to set fractional energy binning [%default]', \
                        dest='fraction_switch', \
                        default=False, \
                        action='store_true')
    parser.add_option('-l', '--lower_energy', \
                        help='minimum energy of first bin in eV [%default]', \
                        dest='lower_energy', \
                        type='float', \
                        default=300.0)
    parser.add_option('-u', '--upper_energy', \
                        help='maximum energy of last bin in eV [%default]', \
                        dest='upper_energy', \
                        type='float', \
                        default=10000.0)
    parser.add_option('-c', '--count_bin', \
                        help='minimum counts per bin [%default]', \
                        dest='bin_min_counts', \
                        type='long', \
                        default=25)

    if (argv == []):
        parser.print_help()
        exit(-1)

    return parser

def valid_args_check(opts, parser):

    if opts.input_pi is None:
        print("The input pi file is missing\n")
        parser.print_help()
        exit(-1)
    if os.path.isfile(opts.input_pi) == 0:
        print("The specified input pi file ("+opts.input_pi+") does not exist\n")
        exit(-1)

    input_fits   = fits.open(opts.input_pi)
    input_pi     = input_fits[1]
    #  input_pi_hdr = input_pi.header.ascardlist()
    input_pi_hdr = input_pi.header

    if opts.output_pi is None:
        print("The output pi file is missing\n")
        parser.print_help()
        exit(-1)

    if opts.rmf is None:
        if 'RESPFILE' in input_pi_hdr:
            opts.rmf = input_pi_hdr['RESPFILE']
            print("Found RESPFILE: {}".format(opts.rmf))
        if opts.rmf is None:
            print("No rmf file is specified in either input pi file or command line\n")
            parser.print_help()
            exit(-1)

    if opts.background is None:
        if 'BACKFILE' in input_pi_hdr:
            opts.background = input_pi_hdr['BACKFILE']
            print("Found BACKFILE: {}".format(opts.background))
        if opts.background is None:
            print("*WARNING*: he grouped spectra has no background file\n")

    if opts.arf is None:
        if 'ANCRFILE' in input_pi_hdr:
            opts.arf = input_pi_hdr['ANCRFILE']
            print("Found ANCRFILE: {}".format(opts.arf))
        if opts.arf is None:
            print("*WARNING*: The grouped spectra has no arf file\n")

    input_fits.close

def group_pha(opts, pi_data, ebounds):
    pi_rows = pi_data.shape[0]

    count_index = np.where(np.array(pi_data.columns.names).astype(str) == 'COUNTS')[0][0] # AKH added this to find the index that corresponded to counts
    quality             = np.zeros(pi_rows)+5
    grouping            = np.zeros(pi_rows)+1

    channels            = ebounds.data.field('CHANNEL')
    e_min               = ebounds.data.field('E_MIN')*1000
    e_max               = ebounds.data.field('E_MAX')*1000
    energy              = np.sqrt(e_min*e_max)
    energy_width        = e_max-e_min
    min_width_at_energy = np.zeros(pi_rows)+opts.bin_min_energy
    
    if (opts.fraction_switch):
        min_width_at_energy = opts.bin_min_energy * energy
        
    tot_counts          = 0
    row_counter         = 0
    last_new_row        = 0
    counts_in_bin       = opts.bin_min_counts
    energy_width_of_bin = np.sum(energy_width)

    for row in (pi_data):
        pi_row = row[0] # This is also just the channels so I'm not sure what the difference is
        pi_counts = row[count_index] # This was originally an index of 2 potentially a difference between chandra v. swift
        matching_index = (np.where(channels==pi_row))[0][0]

        if ((e_min[matching_index] > opts.lower_energy) and (opts.upper_energy > e_max[matching_index])):
            
            quality[row_counter] = 0
            tot_counts += pi_counts

            if((counts_in_bin >= opts.bin_min_counts) and (energy_width_of_bin >= min_width_at_energy[row_counter])):
                counts_in_bin         = pi_counts
                energy_width_of_bin   = energy_width[row_counter]
                last_new_row          = row_counter
                grouping[row_counter] =1
            else:
                counts_in_bin         += pi_counts
                energy_width_of_bin   += energy_width[row_counter]
                grouping[row_counter] =-1

        row_counter += 1
    

    if counts_in_bin < opts.bin_min_counts:
        grouping[last_new_row] =-1
    
    return (quality, grouping, tot_counts)

def fitsio_grppha(opts):
    input_fits   = fits.open(opts.input_pi)

    ebounds = fits.open(opts.rmf)['EBOUNDS']

    (quality, grouping, tot_counts) = group_pha(opts, input_fits[1].data, ebounds)
    if tot_counts < 200: 
        opts.bin_min_counts = 1
        (quality, grouping, tot_counts) = group_pha(opts, input_fits[1].data, ebounds)
    
    quality = fits.Column(name='QUALITY', format='I', array=quality)
    grouping = fits.Column(name='GROUPING', format='I', array=grouping)

    if 'QUALITY' in input_fits[1].columns.names:
        print('Deleting QUALITY column')
        input_fits[1].columns.del_col('QUALITY')

    if 'GROUPING' in input_fits[1].columns.names:
        print('Deleting GROUPING column')
        input_fits[1].columns.del_col('GROUPING')

    print(f'Final Count Grouping: {opts.bin_min_counts}')
    print(f'Total Counts: {tot_counts}')

    output_fits = input_fits.copy() # This is a reference not a unique object despite it claiming to be a unique object
    output_fits[1] = fits.BinTableHDU.from_columns(input_fits[1].columns + quality + grouping, header=input_fits[1].header)
    output_fits[1].header['BACKFILE'] = opts.background
    output_fits[1].header['ANCRFILE'] = opts.arf
    output_fits[1].header['COUNTGRP'] = opts.bin_min_counts
    output_fits[1].header['COUNTTOT'] = tot_counts

    output_fits.writeto(opts.output_pi, overwrite=True)
    output_fits.close()
    input_fits.close()

def main(argv):
    parser = parse_args(argv)
    (opts, args) = parser.parse_args()
    valid_args_check(opts, parser)
    fitsio_grppha(opts)

if __name__ == "__main__":
    main(sys.argv[1:])

#!/usr/bin/env python3

import sys
if not sys.version_info.major == 3: # Paul's being a dick
    sys.exit("Please use Python 3 for this script")

import os
import re
import argparse

class MapFile(object):
    """ A class for the PLINK MAP File"""

    def __init__(self, harvest):
        self._genetic_data = set()
        with open(harvest, 'r') as h:
            for index, line in enumerate(h):
                if index == 0:
                    continue
                split_line = line.strip().split('\t')
                #   The columns we want are
                #   1: SNP ID
                #   14: chromosome (2011 map)
                #   15: cM (2011 map)
                #   18: physical position (new reference genome)
                snpid = split_line[0]
                chrom = split_line[13]
                cm = split_line[14]
                pp = split_line[17]
                if snpid == '' or snpid == 'NA':
                    snpid = '0'
                if chrom == '' or chrom == 'NA':
                    chrom = '0'
                if cm == '' or cm == 'NA':
                    cm = '0'
                if pp == '' or pp == 'NA':
                    pp = '0'
                snp_tuple = (snpid, chrom, cm, pp)
                self._genetic_data.add(snp_tuple)

    def sort_map(self):
        """Sort the genetic map by chromosome, then genetic map position"""
        sorted_snps = sorted(
            self._genetic_data,
            key=lambda tup : (tup[1], float(tup[2]))
        )
        return sorted_snps


class Pedigree(object):
    """A class to hold the information of a .ped for PLINK"""

    _GET_FAMILY = re.compile(r'([A-Z]+)', re.IGNORECASE)
    _GET_GENOTYPES = re.compile(r'[A-Z0-9]')

    def __init__(self, sample):
        family = self._GET_FAMILY.search(sample)
        self._family = family.groups()[0]
        self._individual = sample
        self._paternal = 0
        self._maternal = 0
        self._sex = 0
        self._phenotype = -9
        self._genotypes = dict()
        self._selected_genotypes = []

    def assign_genotype(self, snp, genotype, probability, cutoff=0.8):
        if float(probability) < cutoff:
            gen = ['0', '0']
        else:
            gen = [genotype[0], genotype[1]]
        self._genotypes[snp] = gen

    def select_genotypes(self, snp):
        try:
            self._selected_genotypes += self._genotypes[snp]
        except KeyError:
            raise

    def format_pedigree(self):
        # geno_values = str(self._genotypes.values())
        # genotypes = self._GET_GENOTYPES.findall(geno_values)
        ped_line = [
            str(self._family),
            str(self._individual),
            str(self._paternal),
            str(self._maternal),
            str(self._sex),
            str(self._phenotype)
        ] + self._selected_genotypes
        return '\t'.join(ped_line)


def make_argument_parser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        '-a',
        '--alchemy',
        dest='alchemy',
        type=str,
        required=True,
        default=None,
        metavar='ALCHEMY FILE',
        help='Path to Alchemy file'
    )
    parser.add_argument(
        '-v',
        '--harvest',
        dest='harvest',
        type=str,
        required=True,
        default=None,
        metavar='HARVEST FILE',
        help='Path to HarvEST file with physical positions added to the end of it'
    )
    parser.add_argument(
        '-o',
        '--outname',
        dest='output',
        type=str,
        required=False,
        default='output',
        metavar='OUTPUT BASENAME',
        help="Basename for output files; will place in same directory as alchemy file. Defaults to 'output'"
    )
    return parser


def overwrite(filename):
    message = filename + " exists! Would you like to overwrite (y/n)? "
    error = "Please select (y)es or (n)o"
    isvalid = lambda response : re.search(r'([yn])', response, re.IGNORECASE)
    res = None
    while res == None:
        res = input(str(message))
        if not isvalid(res):
            print(error, file=sys.stderr)
            res = None
    return res.upper()


def main():
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args())
    #   Read in the input data
    try:
        #   HarvEST data
        print("Reading in", args['harvest'], file=sys.stderr)
        mapfile = MapFile(args['harvest'])
        sorted_map = mapfile.sort_map()
        #   Alchemy data
        print("Reading in", args['alchemy'], file=sys.stderr)
        pedigree = dict() # PED FILE
        with open(args['alchemy'], 'r') as a:
            for line in a:
                if line.startswith('#'):
                    continue
                split_line = line.strip().split()
                snp = split_line[0]
                sample = split_line[1]
                genotypes = split_line[3]
                probability = split_line[4]
                if sample not in pedigree.keys():
                    p = Pedigree(sample)
                    pedigree[sample] = p # Add our single pedigree to our pedigree file (DICT)
                pedigree[sample].assign_genotype(snp, genotypes, probability)
    except FileNotFoundError as error:
        sys.exit("Failed to find " + error.filename)
    #   Select the genotypes for our SNPs
    unfound_snps = []
    print("Selecting SNPs from", args['harvest'], file=sys.stderr)
    for sample in pedigree.keys():
        for snp in sorted_map:
            (snpid, chrom, cm, pp) = snp
            try:
                pedigree[sample].select_genotypes(snpid)
            except KeyError:
                unfound_snps.append(sorted_map.pop(snp))
    #   Write the outputs
    print("Writing output files", file=sys.stderr)
    if os.path.dirname(args['alchemy']) == '':
        outdir = os.getcwd()
    else:
        outdir = os.path.dirname(args['alchemy'])
    ped_output = outdir + '/' + args['output'] + '.ped'
    map_output = outdir + '/' + args['output'] + '.map'
    unfound_out = outdir + '/' + 'unfound_snps.log'
    if os.path.isfile(ped_output):
        poo = overwrite(ped_output)
    else:
        poo = 'Y'
    if poo == 'Y':
        with open(ped_output, 'w') as po:
            for sample in sorted(pedigree.keys()):
                po.write(pedigree[sample].format_pedigree())
                po.write('\n')
            po.flush()
    if os.path.isfile(map_output):
        moo = overwrite(map_output)
    else:
        moo = 'Y'
    if moo == 'Y':
        with open(map_output, 'w') as mo:
            for snp in sorted_map:
                (snpid, chrom, cm, pp) = snp
                outline = [chrom, snpid, cm, pp]
                mo.write('\t'.join(outline))
                mo.write('\n')
            mo.flush()
    if len(unfound_snps) > 0:
        print("Failed to find", str(len(unfound_snps)), "snps in", args['alchemy'], file=sys.stderr)
        with open(unfound_out, 'w') as uo:
            for snp in unfound_snps:
                (snpid, chrom, cm, pp) = snp
                uo.write(snpid)
                uo.write('\n')


if __name__ == '__main__':
    main()

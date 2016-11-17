#!/usr/bin/env python3

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")

import os
import argparse

try:
    from overload import overload
except ImportError as error:
    sys.exit("Please install " + error.name)


#   A custom error
class MismatchedSamplesError(Exception):
    """Your samples do not match"""


#   A class for every line in an 012 file
class ZeroOneTwo(object):
    """This is a class for a line of an 012 file from VCFTools"""

    _HIERF_TRANS = {
            0 : '11',
            1 : '12',
            2 : '22',
            -1 : 'NA'
        }

    @overload
    def __init__(self, ind_num, snp_dict, name, levels=None):
        try:
            assert isinstance(ind_num, int)
            assert isinstance(snp_dict, dict)
            for s in snp_dict.values():
                assert s is -1 or s is 0 or s is 1 or s is 2
            assert isinstance(levels, list) or levels is None
            if levels:
                for l in levels:
                    assert isinstance(l, int)
            assert isinstance(name, str)
        except AssertionError:
            raise
        self._num = ind_num
        self._snps = snp_dict
        self._levels = levels
        self._name = name
        self._hier = {}

    @__init__.add
    def __init__(self, arg_dict):
        try:
            assert isinstance(arg_dict, dict)
            assert len(arg_dict) >= 2 and len(arg_dict) <= 4
        except AssertionError:
            raise
        args = ['ind_num', 'snp_dict', 'name', 'levels']
        for a in args:
            if a not in arg_dict:
                arg_dict[a] = None
        self.__init__(
            ind_num=arg_dict['ind_num'],
            snp_dict=arg_dict['snp_dict'],
            name=arg_dict['name'],
            levels=arg_dict['levels']
        )

    def __repr__(self):
        return self._name + ':' + str(len(self._snps))

    def __eq__(self, other):
        if isinstance(other, ZeroOneTwo):
            equal_snps = len(self._snps) == len(other._snps)
            if self._levels:
                equal_levels = len(self._levels) == len(other._levels)
            else:
                equal_levels = True
            return equal_snps and equal_levels
        else:
            raise NotImplementedError("Cannot equate a ZeroOneTwo to " + type(other))

    def get_ind_num(self):
        """Get the number of this individual"""
        return self._num

    def get_num_snps(self):
        """Get the number of snps"""
        return len(self._snps)

    def get_levels(self):
        """Get the levels of this individual"""
        return self._levels

    def get_num_levels(self):
        """Get the number of levels for this individual"""
        return len(self._levels)

    def get_name(self):
        """Get the name of this individual"""
        return self._name

    def convert_to_hierfstat(self):
        """Convert the 012 format to hierfstat format"""
        for snp, geno in self._snps.items():
            self._hier[snp] = self._HIERF_TRANS[geno]

    def format_hierfstat(self, snp_order=None):
        """Format as hierfstat"""
        try:
            assert isinstance(snp_order, list) or snp_order is None
        except AssertionError:
            raise
        if not self._hier:
            self.convert_to_hierfstat()
        if not snp_order:
            snp_order = sorted(self._hier)
        hierf = []
        hierf.append(self._name)
        if self._levels:
            hierf += [str(l) for l in self._levels]
        for snp in snp_order:
            try:
                hierf.append(self._hier[snp])
            except KeyError:
                raise
        return '\t'.join(hierf)


#   Parse the levels
def parse_levels(level_files):
    """Parse any level files"""
    try:
        assert isinstance(level_files, list) or level_files is None
    except AssertionError:
        raise
    if level_files == None:
        return None
    levels_by_pop = {}
    for index, file in enumerate(level_files):
        try:
            print("Level ", index + 1, ": ", file, sep='', end='\n', file=sys.stderr, flush=True)
            with open(file, 'r') as f:
                for line in f:
                    sp = line.strip().split()
                    ind = sp[0]
                    lev = int(sp[1])
                    if ind not in levels_by_pop.keys():
                        levels_by_pop[ind] = [lev]
                    else:
                        levels_by_pop[ind].append(lev)
        except FileNotFoundError:
            raise
    try:
        num_levels = len(level_files)
    except TypeError:
        num_levels = 0
    return(levels_by_pop, num_levels)


#   Parse the SNPs
def parse_snps(geno_list, pos_dict):
    """Parse the SNPs into proper positions"""
    try:
        assert isinstance(geno_list, list)
        for g in geno_list:
            assert g == '0' or g == '1' or g == '2' or g == '-1'
        assert isinstance(pos_dict, dict)
    except AssertionError:
        raise
    snp_dict = {}
    for index, geno in enumerate(geno_list):
        snp = pos_dict[index]
        snp_dict[snp] = int(geno)
    return snp_dict


#   Make an argument parser
def make_argument_parser():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        '-012',
        dest='012',
        type=str,
        required=True,
        default=None,
        metavar='012 FILE',
        help="'.012' file from VCFTools"
    )
    parser.add_argument(
        '-indv',
        dest='indv',
        type=str,
        required=True,
        default=None,
        metavar='INDV FILE',
        help="'.indv' file from VCFTools"
    )
    parser.add_argument(
        '-pos',
        dest='pos',
        type=str,
        required=True,
        default=None,
        metavar='POS FILE',
        help="'.pos' file from VCFTools"
    )
    parser.add_argument(
        '-pop',
        dest='pops',
        type=str,
        required=False,
        nargs='*',
        default=None,
        metavar='POPULATION FILE(S)',
        help="One or more files with samples and their population ('1' or '2') listed"
    )
    return parser


#   Run the program here
def main():
    """Run the program"""
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args())
    try:
        inds = {}
        poss = {}
        genos = []
        no_match = []
        levels_by_pop, num_levels = parse_levels(args['pops'])
        with open(args['indv'], 'r') as i:
            for index, line in enumerate(i):
                inds[index] = line.strip()
        with open(args['pos'], 'r') as p:
            for index, line in enumerate(p):
                sp = line.strip().split()
                p = sp[0] + ':' + sp[1]
                poss[index] = p
        with open(args['012'], 'r') as f:
            for index, line in enumerate(f):
                try:
                    sp = line.strip().split()
                    num = sp[0]
                    # snps = [int(s) for s in sp[1:]]
                    snps = parse_snps(sp[1:], poss)
                    name = inds[index]
                    if levels_by_pop:
                        levels = levels_by_pop[name]
                    else:
                        levels = None
                    genos_dict = {
                        'ind_num' : int(num),
                        'snp_dict' : snps,
                        'name' : name,
                        'levels' : levels
                    }
                    z = ZeroOneTwo(arg_dict=genos_dict)
                    genos.append(z)
                except KeyError:
                    print("No match for", name, file=sys.stderr)
                    no_match.append(name)
                    continue
    except FileNotFoundError as error:
        sys.exit("Failed to find " + error.filename)
    outbase = os.path.basename(args['012'])
    outfile = os.getcwd() + '/' + os.path.splitext(outbase)[0] + '.hierfstat'
    with open(outfile, 'w') as o:
        print("Writing output to", outfile, file=sys.stderr)
        snp_order = sorted(poss.values())
        o.write('\t')
        for num in range(num_levels):
            o.write('Level ' + str(num + 1))
            o.write('\t')
        o.write('\t'.join(snp_order))
        o.write('\n')
        for ind in genos:
            o.write(ind.format_hierfstat(snp_order))
            o.write('\n')


if __name__ == '__main__':
    main()

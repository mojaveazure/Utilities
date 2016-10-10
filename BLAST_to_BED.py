#!/usr/bin/env python3
"""Run a BLAST search and create a BED file from the resulting hits"""

import sys
if sys.version_info.major is not 3:
    sys.exit("Please use Python 3 for this script")


import os
import re
import argparse

try:
    from bs4 import BeautifulSoup
    from bs4 import element
    from bs4 import FeatureNotFound
    from Bio.Blast.Applications import NcbiblastnCommandline
except ImportError as error:
    sys.exit("Please install " + error.name)


DEFAULT_OUTPUT = os.getcwd() + '/output.bed'

#   An error for no reference provided
class NoReferenceError(Exception):
    """A reference was not provided"""


#   An error I probably overuse...
class NoHitError(Exception):
    """A SNP has not been found"""


#   A custom error
class BLASTFailedError(Exception):
    """BLAST seems to have failed..."""


#   A class definition for a BLAST hit
class Hit(object):
    """This is a class for a BLAST hit
    It continas the following information:
        Chromosome name
        Query name
        Hit e-value
        Query sequence
        Subject sequence
        Hit start relative to subject
        Hit end relative to subject
        Subject strand (forward or reverse)
        SNP Position relative to subject (once calcualted with Hit.get_snp_position())
        """
    def __init__(self, chrom, name, evalue, hstart, hend, hstrand):
        try:
            assert isinstance(chrom, str)
            assert isinstance(name, str)
            assert isinstance(evalue, float)
            assert isinstance(hstart, int)
            assert isinstance(hend, int)
            assert isinstance(hstrand, int)
        except AssertionError:
            raise TypeError
        try:
            assert hstrand == 1 or hstrand == -1
        except AssertionError:
            raise ValueError
        self._chrom = chrom
        self._name = name
        self._evalue = evalue
        self._start = hstart
        self._end = hend
        self._hstrand = hstrand

    def __repr__(self):
        return self._name + ":" + str(self._evalue)

    def get_chrom(self):
        """Get the chromosome that the hit matched to"""
        return self._chrom

    def get_name(self):
        """Get the query name of the hit"""
        return self._name

    def get_rc(self):
        """Did the query align to the forward (False) or reverse (True) strand"""
        return self._hstrand == -1

    def format_bed(self):
        """Format a BED file"""
        bed = [
            self._chrom
        ]
        if self.get_rc():
            bed += [str(self._end - 1), str(self._start)]
        else:
            bed += [str(self._start - 1), str(self._end)]
        bed.append(self._name)
        return '\t'.join(bed)


#   A class definition for holding Iterations
class SNPIteration(object):
    """This is a class for holding BLAST iterations
    It holds the SNP name and a list of hits"""
    @staticmethod
    def GET_VALUE(tag, value):
        try:
            assert isinstance(tag, element.Tag)
        except AssertionError:
            raise TypeError
        return tag.findChild(value).text

    _VALS = ['Hsp_evalue', 'Hsp_hit-from', 'Hsp_hit-to', 'Hsp_hit-frame']
    def __init__(self, iteration):
        try:
            assert isinstance(iteration, element.Tag)
        except AssertionError:
            raise TypeError
        self._snpid = SNPIteration.GET_VALUE(iteration, 'Iteration_query-def')
        self._hits = []
        self._fail = False
        #   Start parsing hits
        for hit in iteration.findAll('Hit'):
            self._parse_hit(hit)
        # If we don't have any hits, set 'self._fail' to True
        if len(self._hits) < 1:
            self._fail = True

    def __repr__(self):
        return self._snpid + '(' + str(len(self._hits)) + ' hit(s))'

    def _parse_hsp(self, hsp):
        """Parse the hsp section of a BLAST XML"""
        try:
            assert isinstance(hsp, element.Tag)
        except AssertionError:
            raise TypeError
        #   Collect all values
        hsp_vals = []
        for val in self._VALS:
            hsp_vals.append(SNPIteration.GET_VALUE(hsp, val))
        # hsp_vals = map(SNPIteration.GET_VALUE, repeat(hsp, len(self._VALS)), self._VALS)
        #   Return as a tuple
        return tuple(hsp_vals)

    def _parse_hit(self, hit):
        """Parse the hit section of a BLAST XML"""
        try:
            assert isinstance(hit, element.Tag)
        except AssertionError:
            raise TypeError
        chrom = SNPIteration.GET_VALUE(hit, 'Hit_def') # Get the chromosome information
        hsps = [] # A list to hold hsps
        for hsp in hit.findAll('Hsp'):
            try:
                #   Try to parse the hsp
                hsp_vals = self._parse_hsp(hsp)
                hsps.append(hsp_vals)
            except NoHitError: # If there's no gap, skip
                continue
        for hsp in hsps:
            #   Unpack our tuple
            (evalue, hit_start, hit_end, strand) = hsp
            #   Make a Hit
            hit = Hit(
                chrom=chrom,
                name=self._snpid,
                evalue=float(evalue),
                hstart=int(hit_start),
                hend=int(hit_end),
                hstrand=int(strand)
            )
            #   Add our Hit to the list of Hits
            self._hits.append(hit)

    def get_snpid(self):
        """Get the SNP ID for this iteration"""
        return self._snpid

    def check_fail(self):
        """See if we are lacking any hits"""
        if self._fail:
            raise NoHitError

    def format_bed(self):
        try:
            self.check_fail()
        except NoHitError:
            raise
        bed_lines = []
        for hit in self._hits:
            bed_lines.append(hit.format_bed())
        return bed_lines


#   Ensure we have our BLAST database
def validate_db(db_path):
    """Find a BLAST nucleotide database"""
    try:
        assert isinstance(db_path, str)
    except AssertionError:
        raise TypeError
    db_name = os.path.basename(db_path)
    nhr = re.compile(r'(%s\.[0-9\.]*nhr)' % db_name).search
    nin = re.compile(r'(%s\.[0-9\.]*nin)' % db_name).search
    nsq = re.compile(r'(%s\.[0-9\.]*nsq)' % db_name).search
    nal = re.compile(r'(%s\.*nal)' % db_name).search
    db_directory = os.path.abspath(os.path.realpath(os.path.dirname(db_path)))
    if not db_directory:
        db_directory = os.getcwd()
    print('Searching for proper database files for', db_name, 'in', db_directory, file=sys.stderr)
    db_contents = '\n'.join(os.listdir(db_directory))
    if not nhr(db_contents) and not nin(db_contents) and not nsq(db_contents) and not nal(db_contents):
        raise FileNotFoundError("Failed to find the BLAST nucleotide database")


#   A function to run BLASTn
def run_blastn(query, database, evalue, max_seqs, max_hsps):
    """Run BLASTn"""
    try:
        assert isinstance(query, str)
        assert isinstance(database, str)
        assert isinstance(evalue, float)
        assert isinstance(max_seqs, int)
        assert isinstance(max_hsps, int)
    except AssertionError:
        raise TypeError
    #   Create an output name
    query_base = os.path.basename(os.path.splitext(query)[0])
    database_base = os.path.basename(os.path.splitext(database)[0])
    blast_out = os.getcwd() + '/' + query_base + '_' + database_base + '_BLAST.xml'
    validate_db(database)
    #   Setup BLASTn
    blastn = NcbiblastnCommandline(
        query=query,
        db=database,
        evalue=evalue,
        outfmt=5,
        max_target_seqs=max_seqs,
        max_hsps=max_hsps,
        out=blast_out
    )
    #   Run BLASTn
    print(blastn, file=sys.stderr)
    blastn()
    if not os.path.exists(blast_out):
        raise BLASTFailedError
    return blast_out


#   Make an argument parser
def make_argument_parser():
    parser = argparse.ArgumentParser(add_help=True)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument( # Take FASTA as input
        '-f',
        '--fasta',
        dest='fasta',
        type=str,
        default=None,
        metavar='FASTA FILE',
        help="Input FASTA file to run BLAST on, incompatible with '-x | --xml'"
    )
    inputs.add_argument( # Take XML as input
        '-x',
        '--xml',
        dest='xml',
        type=str,
        default=None,
        metavar='XML FILE',
        help="Input BLAST XML file to turn into BED file, incompatible with '-f | --fasta'"
    )
    outputs = parser.add_mutually_exclusive_group(required=False)
    outputs.add_argument( # Append to preexisting BED file
        '-b',
        '--bed',
        dest='bed',
        type=str,
        default=None,
        metavar='BED FILE',
        help="BED file to append results to, incompatible with '-o | --outfile'"
    )
    outputs.add_argument( # Write to new BED file
        '-o',
        '--outfile',
        dest='outfile',
        type=str,
        default=DEFAULT_OUTPUT,
        metavar='OUTPUT FILE',
        help="Name of output file to write to, defaults to '" + DEFAULT_OUTPUT + "', incompatible with '-b | --bed'"
    )
    parser.add_argument( # BLAST database to BLAST against
        '-d',
        '--database',
        dest='database',
        type=str,
        required=False,
        default=None,
        metavar='REFERENCE BLAST DATABASE',
        help="Reference BLAST database in nucleotide format, used only when running BLAST"
    )
    parser.add_argument( # E-value threshold
        '-e',
        '--evalue',
        dest='evalue',
        type=float,
        required=False,
        default=1e-1,
        metavar='E-VALUE THRESHOLD',
        help="Evalue threshold for BLAST, defaults to '1e-1', used only when running BLAST"
    )
    parser.add_argument( # Maximum number of hits
        '-s',
        '--max-hits',
        dest='max_hits',
        required=False,
        type=int,
        default=1,
        metavar='MAX HITS',
        help="Maximum hits per query, defaults to '1', used only when running BLAST"
    )
    parser.add_argument( # Maximum number of HSPs
        '-m',
        '--max-hsps',
        dest='max_hsps',
        required=False,
        type=int,
        default=1,
        metavar='MAX HSPS',
        help="Maximum HSPs per hit, defaults to '1', used only when running BLAST"
    )
    parser.add_argument( # Do we keep the XML file from BLAST?
        '-k',
        '--keep-xml',
        dest='keep_xml',
        required=False,
        action='store_const',
        const=True,
        default=False,
        metavar='KEEP XML',
        help="Do we keep the XML results? pas '-k | --keep-xml' to say 'yes', used only when running BLAST"
    )
    return parser


#   Run the program
def main():
    """Run BLAST_to_BED.py"""
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    # args = vars(parser.parse_args())
    args = {key : value for key, value in vars(parser.parse_args()).items() if value is not None}
    try:
        #   Are we running a BLAST search given a FASTA?
        if 'fasta' in args.keys():
            if not args['database']:
                raise NoReferenceError
            print("BLASTing", args['fasta'], "against", args['database'], file=sys.stderr)
            blast_xml = run_blastn(
                query=args['fasta'],
                database=args['database'],
                evalue=args['evalue'],
                max_seqs=args['max_hits'],
                max_hsps=args['max_hsps']
            )
        #   Or are we given a BLAST XML file?
        elif 'xml' in args.keys():
            print("Using", args['xml'], "as input", file=sys.stderr)
            blast_xml = args['xml']
        else:
            raise NotImplementedError("Whatever you're trying to do, we don't do yet...")
        #   Read in the XML as soup
        blast_soup = BeautifulSoup(open(blast_xml, 'r'), 'xml')
        #   Make a dictionary to hold iterations and find all iterations
        # iterations = []
        if args['bed']:
            print("Appending to", args['bed'], file=sys.stderr)
            no_hit_name = os.path.basename(args['bed']) + '_failed.log'
            outhandle = open(args['bed'], 'a')
        else:
            print("Writing to", args['outfile'], file=sys.stderr)
            no_hit_name = os.path.basename(args['outfile']) + '_failed.log'
            outhandle = open(args['outfile'], 'w')
    except FeatureNotFound:
        sys.exit("Pleast install 'lxml' to properly parse the BLAST results")
    except NoReferenceError:
        sys.exit(NoReferenceError)
    except FileNotFoundError as error:
        sys.exit("Cannot find " + error.filename)
    no_hit = []
    for query in blast_soup.findAll('Iteration'):
        try:
            iteration = SNPIteration(query)
            # iterations.append(iteration)
            bed_lines = iteration.format_bed()
            outhandle.write('\n'.join(bed_lines))
            outhandle.write('\n')
        except NoHitError:
            print("No hit found for", iteration.get_snpid(), file=sys.stderr)
            no_hit.append(iteration.get_snpid())
            continue
    #   If we make our own BLAST XML and we aren't told to keep it
    outhandle.close()
    if len(no_hit) > 0:
        print("Writing", len(no_hit), "failed searches to", no_hit_name, file=sys.stderr)
        with open(no_hit_name, 'w') as n:
            for fail in no_hit:
                n.write(fail)
                n.write('\n')
    if args['fasta'] and not args['keep_xml']:
        print("Removing intermediate files", file=sys.stderr)
        os.remove(blast_xml) # Remove


if __name__ == '__main__':
    main()



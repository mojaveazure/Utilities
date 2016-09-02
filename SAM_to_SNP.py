#!/usr/bin/env python3
"""Script to find SNP positions based off of mapped alignment (SAM File) and capture design"""

#   Borrowed heavily from Thomas Kono

import sys
if not sys.version_info.major == 3:
    sys.exit("Please use Python 3 for this script")

import re
import argparse

try:
    from Bio import Seq
    from Bio import SeqIO
except ImportError:
    sys.exit("Please install BioPython for this script")


class NotABaseError(Exception):
    """You have not provided a single-character string"""


#   A class definition for a SNP
class SNP(object):
    """A class to hold VCF information about an individual SNP. Stores the
    following information:
        SNP ID
        Contig
        Reference Position
        Reference Base
        Alternate Base
    """
    @staticmethod
    def reverse_complement(base):
        """Get the reverse complement"""
        try:
            assert isinstance(base, str)
            assert len(base) is 1
            rc = str.maketrans('ACGT', 'TGCA') # Traslation table for reverse complentary sequences
            return base.translate(rc)
        except AssertionError:
            raise NotABaseError

    def __init__(self, lookup, alignment, reference):
        try:
            #   Ensure we've been given a lookup and alignment object
            assert isinstance(lookup, Lookup)
            assert isinstance(alignment, Alignment)
            assert isinstance(reference, dict)
        except AssertionError:
            raise
        self._snpid = lookup.get_snpid() # The SNP ID is the same as the one in the lookup
        self._contig = alignment.get_contig() # The contig is found in the alignment
        #   Everything else gets made in a bit
        self._position = None
        self._reference = None
        self._alternate = None
        #   Get the rest of the information
        self._calculate_position(lookup, alignment) # True SNP position
        self._find_states(lookup, alignment, reference) # Reference and alternate states

    def __repr__(self):
        return self._snpid

    def __eq__(self, other):
        if isinstance(other, SNP):
            return self._snpid == other._snpid
        elif isinstance(other, str):
            return self._snpid == other
        else:
            return False

    def _calculate_position(self, lookup, alignment):
        """Calculate the position of the SNP in the reference sequence"""
        index = 0 # Index of our split CIGAR string
        if alignment.get_rc() or lookup.get_rc(): # If we're reverse complementing
            qpos = lookup.get_reverse_position() - 1 # Start with the reverse position of the SNP, must subtract one
        else: # Otherwise
            qpos = lookup.get_forward_position() # Start with the forward posittion
        while True: # Endless loop to do weird things...
            try: # While we have a CIGAR string to parse
                old = qpos # Store our previously calculated SNP position
                #   Seach the CIGAR string as a list, starting with index 0, for indels
                if re.search('M', alignment.get_cigar()[index]): # If we have a perfect match
                    if qpos < int(''.join(re.findall(r'\d+', alignment.get_cigar()[index]))): # If our SNP is in the perfect match
                        break # Exit the loop, we have our position
                if re.search('D', alignment.get_cigar()[index]): # If we have a deletion relative to reference
                    qpos += int(''.join(re.findall(r'\d+', alignment.get_cigar()[index]))) # Add the deletion to our SNP position
                if re.search('[IS]', alignment.get_cigar()[index]): # If we have an insertion relative to reference
                    qpos -= int(''.join(re.findall(r'\d+', alignment.get_cigar()[index]))) # Subtract the insertion from our SNP postion
                index += 1 # Increase the index
                if qpos <= 0 or qpos >= lookup.get_length(): # If we've gone beyond the scope of our lookup: 0 is before the sequence, lookup.get_length() is after
                    qpos = old # Go back to our previously calculated SNP postion
                    break # Exit the loop, we have our position
            except IndexError: # If we run out of CIGAR string codes
                break # Exit the loop, we have our position
        self._position = alignment.get_position() + qpos # Our SNP position is at the mapping position plus the SNP position

    def _find_states(self, lookup, alignment, reference):
        """Get the reference and alternate alleles"""
        #   Get the reference allele, given our contig and position found above
        self._reference = reference[self._contig][self._position - 1] # Subtract one as FASTA is 1-based and Python is 0-based
        if alignment.get_rc(): # If we're reverse complement
            alt, do_rc = lookup.get_alternate(self.reverse_complement(self._reference))
            self._alternate = self.reverse_complement(alt)
        else:
            self._alternate, do_rc = lookup.get_alternate(self._reference) # An 'N' will be returned if the reference allele doesn't match with our IUPAC code
        if do_rc:
            self._reference = self.reverse_complement(self._reference)

    def check_masked(self):
        """Check to see if our alternate allele is masked"""
        if self._alternate == 'N': # If our alternate allele is masked, or an 'N'
            return True # Return True
        else: # Otherwise
            return False # Return False

    def format_vcf(self):
        """Format the information in VCF style"""
        #   Create a list of information for a VCF file
        vcf_line = [
            self._contig,
            str(self._position),
            self._snpid,
            self._reference,
            self._alternate,
            '.',
            '.',
            's'
            ]
        return '\t'.join(vcf_line) # Join everything together with a tab


#   A class definition for a SAM alignment
class Alignment(object):
    """A class to hold a SAM Alignment
    It contains the following information:
        Query Name
        Bitwise Flag
        Reference sequence name
        1-based leftmost mapping position
        Cigar string
    """
    _CIGAR = re.compile(u'([0-9]+[A-Z]+)') # Regex to break the CIGAR string into component codes using re.findall()
    def __init__(self, line):
        #   There's only some information that we need for our alignment, everything else is forgotten
        split_line = line.strip().split() # Remove leading and trailing whitespace, then split the line by column
        self._qname = split_line[0] # First column in a SAM file
        self._flag = int(split_line[1]) # Second column
        self._rname = split_line[2] # Third column
        self._pos = int(split_line[3]) # Fourth column, should be an int
        self._cigar = self._CIGAR.findall(split_line[5]) # Sixth column, after breaking up the CIGAR string into a list of component codes

    def __repr__(self):
        return self._rname + ':' + self._qname

    def get_rc(self):
        """Check to see if we're reverse complementing our sequence"""
        #   If the 16th bit is set, it's reverse complement
        return self._flag is 16

    def get_name(self):
        """Get the alignment name"""
        return self._qname

    def get_position(self):
        """Get the alignment position"""
        return self._pos

    def get_contig(self):
        """Get the reference sequence name"""
        return self._rname

    def get_cigar(self):
        """Get the CIGAR string as a list"""
        return self._cigar

    def check_flag(self):
        """Make sure we don't have extraneous alignments"""
        return self._flag is 0 or self._flag is 16


#   A class definition for a lookup sequence in Illumina format
class Lookup(object):
    """This is a class for a SNP lookup sequence in Illumina format
    It contains the following information:
        SNP ID
        Sequence with Illumina syntax
        Sequence in IUPAC codes
        Sequence Length
        SNP Position from forward
        SNP Position from reverse
    """
    # A dictionary of IUPAC codes for SNPs
    _IUPAC_CODES = {
        'R' : 'AG',
        'Y' : 'CT',
        'S' : 'CG',
        'W' : 'AT',
        'K' : 'GT',
        'M' : 'AC'
    }
    def __init__(self, snpid, sequence):
        #   We're given the SNP ID and sequence when making the object, everything else 
        #   can be made with _capture_snp() and _find_iupac()
        self._snpid = snpid
        self._sequence = sequence
        self._forward_position = None
        self._reverse_position = None
        self._snp = None
        self._code = None
        self._iupac = None
        self._length = None
        self._rc = False
        #   Get the rest of the information we need for our lookup
        self._capture_snp()
        self._find_iupac()

    def __repr__(self):
        return self._snpid + ':' + self._code

    def _capture_snp(self):
        """Capture the SNP and it's position from the start and end of the sequence"""
        #   Get the forward position
        self._forward_position = self._sequence.find('[')
        #   Get the reverse position
        self._reverse_position = len(self._sequence) - self._sequence.find(']')
        #   Get the SNP
        self._snp = self._sequence[self._forward_position:self._sequence.find(']') + 1]

    def _find_iupac(self):
        """Create an IUPAC version of the sequence and calculate it's length"""
        #   Create a string of the two states of the SNP in alphabetical order
        ordered_snp = ''.join(sorted(re.findall('[ACGT]', self._snp)))
        #   Find the IUPAC code for the SNP
        self._code = ''.join([c for c, o in self._IUPAC_CODES.items() if ordered_snp == o])
        #   Create the IUPAC version of the sequence
        self._iupac = re.sub(r'\[%s\]' % self._snp[1:-1], self._code, self._sequence)
        #   Calculate the length of the sequence
        self._length = len(self._iupac)

    def set_rc(self):
        """Set the lookup sequence to reverse complement"""
        self._rc = True

    def get_snpid(self):
        """Get the SNP ID"""
        return self._snpid

    def get_forward_position(self):
        """Get the SNP position in the forward direction"""
        return self._forward_position

    def get_reverse_position(self):
        """Get the SNP position in the reverse direction"""
        return self._reverse_position

    def get_length(self):
        """Get the length of the IUPAC sequence"""
        return self._length

    def get_rc(self):
        """See if the lookup is reverse complement"""
        return self._rc

    #   Search the IUPAC codes for an alternate allele of a SNP
    def get_alternate(self, reference):
        """Get the alternate allele given an IUPAC code and reference allele"""
        ref = re.compile(u'(%s)' % reference) # Regex to ensure that our found reference allele is covered by the IUPAC code
        rc = re.compile(u'(%s)' % SNP.reverse_complement(reference)) # Regex to see if our reference is actually the reverse complement
        alt = re.compile(u'([^%s])' % reference) # Regex to find the alternate allele
        alt_rc = re.compile(u'([^%s])' % SNP.reverse_complement(reference)) # Regex to find the alternate allele to our reverse complementary sequence
        if ref.search(self._IUPAC_CODES[self._code]): # If our reference allele is plausible given our IUPCA code
            alternate = alt.search(self._IUPAC_CODES[self._code]).group() # Get the alternate
            return alternate, False
        elif rc.search(self._IUPAC_CODES[self._code]):
            alternate = alt_rc.search(self._IUPAC_CODES[self._code]).group()
            return alternate, True
        else: # Otherwise, give an 'N'
            return 'N', False

    def format_fasta(self):
        fasta = [
            '>' + self._snpid,
            self._iupac
        ]
        return '\n'.join(fasta)


#   Make an argument parser
def make_argument_parser():
    """Set and parse the arguments"""
    arguments = argparse.ArgumentParser(
        add_help=True,
        usage='%(prog)s -s <sam_file> -l <lookup_table> -r <reference_fasta> [-o <outname>]',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Provide a SAM file of SNP contextual sequences mapped\nto a reference genome, an Illumina lookup table,\nand a reference FASTA file\n\nThe lookup table should have two columns in this order:\n\tSNP ID\n\tIllumina-formatted sequence with SNP\n\nNo headers are allowed in the lookup table\n'
    )
    arguments.add_argument(
        '-s',
        '--sam-file',
        dest='samfile',
        type=str,
        default=None,
        required=True,
        metavar='SAM FILE',
        help="Name of SAM file"
    )
    arguments.add_argument(
        '-l',
        '--lookup-table',
        dest='lookup',
        type=str,
        default=None,
        required=True,
        metavar='LOOKUP TABLE',
        help="Name of lookup table"
    )
    arguments.add_argument(
        '-r',
        '--reference',
        dest='reference',
        type=str,
        default=None,
        required=True,
        metavar='REFERENCE FASTA',
        help="Path to reference fasta file"
    )
    arguments.add_argument(
        '-o',
        '--outname',
        dest='outname',
        type=str,
        default='output',
        required=False,
        metavar='OUTPUT NAME',
        help="Name of output file, without suffix. Defaults to 'output'"
    )
    return arguments


#   Run the program
def main():
    #   Make an argument parser
    parser = make_argument_parser()
    if not sys.argv[1:]: # If we're missing arguments
        sys.exit(parser.print_help()) # Print the help message and exit
    args = vars(parser.parse_args()) # Create a dictionary out of our arguments
    lookup_dict = {} # Create a dictionary to hold lookups
    alignment_dict = {} # Create a dictionary to hold alignments
    SNPs = [] # Create a list of SNPs
    unmapped = [] # Create a list of unmapped designs
    masked = [] # Create a list of SNPs that have a masked alternate allele
    #   Create a header for our VCF
    header = '##fileformat=VCFv4.2\n##INFO<ID=s,Number=1,Type=Flag,Description="Variant is calculated from SAM">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    try: # Attempt to read in our data
        print("Reading in reference FASTA " + args['reference'], file=sys.stderr)
        reference = SeqIO.to_dict(SeqIO.parse(args['reference'], 'fasta')) # Read in our reference 
        print("Reading lookup table " + args['lookup'], file=sys.stderr)
        with open(args['lookup'], 'r') as f: # Read in our lookup table
            for line in f: # For every line in our lookup table
                split = line.strip().split() # Remove leading and trailing whitespace and split the columns
                l = Lookup(split[0], split[1]) # Column 0 should be the SNP ID and column 1 should be the sequence
                lookup_dict[l.get_snpid()] = l # Add the lookup informaiton to our dictionary, using the SNP ID as our key
        print("Reading SAM file " + args['samfile'], file=sys.stderr)
        with open(args['samfile'], 'r') as s: # Read in our SAM file
            for line in s: # For every line in our SAM file
                if line.startswith('@'): # If this is a header line
                    continue # Skip it
                else: # Otherwise
                    a = Alignment(line) # Create an alignment from the line
                    if a.check_flag():
                        alignment_dict[a.get_name()] = a # Add the alignment to our alignment dictionary
    except FileNotFoundError as error: # If any files were not found
        sys.exit("Failed to find " + error.filename) # Exit with error
    for l in lookup_dict: # For every lookup we have
        if l in alignment_dict.keys() and alignment_dict[l].get_contig() is not '*': # If we have an alignment for it and the alignment mapped back
            try:
                s = SNP(lookup_dict[l], alignment_dict[l], reference) # Find information about the SNP
            except AssertionError:
                sys.exit("Something happened with creating a SNP object...")
            if s.check_masked(): # Check to see if the alternate allele was masked
                lookup_dict[l].set_rc()
                new_s = SNP(lookup_dict[l], alignment_dict[l], reference)
                if new_s.check_masked():
                    masked.append(s) # If so, add SNP to list of masked SNPs
                    continue # Don't do anything else
                else:
                    s = new_s
            SNPs.append(s) # Append the SNP to our list of SNPs
        else: # Otherwise
            unmapped.append(l) # Append lookup SNP ID to our list of unmapped
    if len(SNPs) < 1:
        print("Failed to find any SNPs in " + args['samfile'] + " that were designed in " + args['lookup'], file=sys.stderr)
    else:
        #   Write our output file
        outname = args['outname'] + '.vcf'
        out = open(outname, 'w')
        print("Writing SNPs to " + outname, file=sys.stderr)
        out.write(header)
        for snp in SNPs:
            out.write(snp.format_vcf())
            out.write('\n')
        out.close()
    if len(masked) > 0:
        # Write any masked SNPs to a VCF
        maskedname = args['outname'] + '_masked.vcf'
        maskedfile = open(maskedname, 'w')
        print("Writing " + str(len(masked)) + " masked SNPs to " + maskedname, file=sys.stderr)
        maskedfile.write(header)
        for snp in masked:
            maskedfile.write(snp.format_vcf())
            maskedfile.write('\n')
        maskedfile.close()
    if len(unmapped) > 0:
        #   Write any unmapped SNPs to a log file
        print("Failed to map " + str(len(unmapped)) + " SNPs back to the reference", file=sys.stderr)
        print("Writing unmapped SNP IDs to unmapped.log", file=sys.stderr)
        un = open('unmapped.log', 'w')
        for fail in unmapped:
            un.write(fail)
            un.write('\n')
        un.close()


if __name__ == '__main__':
    main()

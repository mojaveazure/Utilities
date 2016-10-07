#!/usr/bin/env python3

import sys
if not sys.version_info.major == 3:
    sys.exit("Please use Python 3 for this script")

import re
import os
import argparse

try:
    import regex
except ImportError as error:
    sys.exit("Failed to load " + error.name)


PROTEIN_DICTIONARY = {
    'F': ['TTT', 'TTC'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'],
    'M': ['ATG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'],
    'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['TGT', 'TGC'],
    'W': ['TGG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'STOP': ['TAA', 'TAG', 'TGA']
}

#   A class for mutants from an amino acid lookup table
class AminoAcidMutant(object):
    """This is a class for mutant amino acids
    It holds the following information:
        GeneID
        Protein
        Position of amino acid in protein
        Reference amino acid
        Mutant amino acid
    """

    def __init__(self, protein_name, aa_position, reference, mutant):
        try:
            assert isinstance(protein_name, str)
            assert isinstance(aa_position, int)
            assert isinstance(reference, str)
            assert isinstance(mutant, str)
            assert len(reference) == 1
            assert len(mutant) == 1
        except AssertionError:
            raise
        self._protein_name = protein_name
        self._aa_position = aa_position
        self._reference_aa = reference
        self._mutant_aa = mutant

    def __repr__(self):
        return self._protein_name + ' (' + self._reference_aa + '/' + self._mutant_aa + ')'

    def get_protein(self):
        """Return the name of the protein"""
        return self._protein_name

    def get_position(self):
        """Return the position of the amino acid"""
        return self._aa_position

    def get_reference(self):
        """Get the reference amino acid"""
        return self._reference_aa.upper()

    def get_mutant(self):
        """Get the mutant amino acid"""
        return self._mutant_aa.upper()

    def format_fail(self):
        """Format for the failed log"""
        fail_list = [
            self._protein_name,
            self._reference_aa,
            self._mutant_aa,
            str(self._aa_position)
        ]
        return '\t'.join(fail_list)


#   A class for lines of a GFF file
class Annotation(object):
    """This is a class for annotation information
    It holds the following information:
        Contig name
        Start position
        End position
        Gene ID
    """

    def __init__(self, start, end, geneid):
        try:
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert isinstance(geneid, str)
        except AssertionError:
            raise
        self._start = start
        self._end = end
        self._geneid = geneid

    def __repr__(self):
        return self._geneid + ':' + str(self._start) + '-' + str(self._end)

    def get_start(self):
        """Return the start of the region"""
        return self._start

    def get_end(self):
        """Return the end of the region"""
        return self._end

    def get_geneid(self):
        """Return the gene id"""
        return self._geneid


#   An object for genes from a gene model
class Gene(object):
    """
    This is a class for a gene model
    It holds the following information:
        Gene ID
        Contig/Chromosome
        Start position of contig
        End position of contig
        Sequence"""

    def __init__(self, geneid, contig, start, end, sequence):
        try:
            assert isinstance(geneid, str)
            assert isinstance(contig, str)
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert isinstance(sequence, str)
        except AssertionError:
            raise
        self._geneid = geneid
        self._contig = contig
        self._start = start
        self._end = end
        self._sequence = sequence

    def __repr__(self):
        return self._contig + ':' + self._geneid

    def get_geneid(self):
        """Get the gene name"""
        return self._geneid

    def get_contig(self):
        """Get the contig"""
        return self._contig

    def get_sequence(self):
        """Get the sequence"""
        return self._sequence


#   A class for holding SNP information
class SNP(object):
    """A class to hold VCF information about an individual SNP. Stores the
    following information:
        SNP ID
        Contig
        Reference Position
        Reference Base
        Alternate Base
    """

    _REVERSE_COMPLEMENT = str.maketrans('ACGT', 'TGCA') # Translation table for reverse complentary sequences

    def __init__(self, snpid, contig, position, reference, alternate):
        try:
            assert isinstance(snpid, str)
            assert isinstance(contig, str)
            assert isinstance(position, int)
            assert isinstance(reference, str)
            assert isinstance(alternate, str)
            assert len(reference) == 1
            assert len(alternate) == 1
        except AssertionError:
            raise
        self._snpid = snpid
        self._contig = contig
        self._position = position
        self._reference = reference
        self._alternate = alternate

    def __repr__(self):
        return self._snpid

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


#   A function to create a gene object from a FASTA record
def parse_gene(gene_model):
    """Parse a line of a FASTA gene model to create a gene object"""
    try:
        #   Two regexes to find stuff
        find_region = re.compile(r'([A-Z0-9]+:\d+-\d+)', re.IGNORECASE)
        find_sequence = re.compile(r'([AGCT]+)$')
        #   The Gene ID is the first thing
        geneid = gene_model.split()[0]
        #   Create a list with the region section parsed
        region_list = find_region.search(gene_model).groups()[0].split(':')
        contig = region_list[0] # The contig is first in the list
        start = int(region_list[1].split('-')[0]) # Start is in the second bit
        end = int(region_list[1].split('-')[1]) # End is in the third bit
        #   Get the sequence using our regex
        sequence = find_sequence.search(gene_model).groups()[0]
        #   Create a gene object
        g = Gene(geneid, contig, start, end, sequence)
        return(g.get_geneid(), g) # And return, with the geneid
    except:
        raise


#   A function to create an annotation from a line of a GFF file
def parse_gff(gff_line):
    """Parse a line of a GFF file"""
    #   Regex to find Gene ID in a GFF
    find_geneid = re.compile(r'Parent=([A-Z0-9\.]+)', re.IGNORECASE)
    gff_list = gff_line.split() # Split the line
    start = int(gff_list[3]) # Start is fourth
    end = int(gff_list[4]) # End is fifth
    geneid = find_geneid.search(gff_list[-1]).groups()[0] # Gene ID is somewhere in the back
    a = Annotation(start, end, geneid) # Make an annotation
    return a # Return the annotation


#   Find a SNP
def find_snp_position(codon, alternate):
    """Find SNP positions given a codon and alternate codon"""
    #   Use regex, not re, as regex allows for one mismatch
    if codon == alternate:
        return -8
    if regex.search(r"(%s){s<=1}"%alternate, codon): # If we have one or fewer mistmatches
        for i in range(3): # For every base in our codon
            if alternate[i] is not codon[i]:
                position = i
        return position # Return the mismatch
    else: # Otherwise
        return -8 # Return bogus


#   A function to create SNP objects
def get_snps(mutant, model, annotation_list):
    #   Figure out where in the model our mutation happens
    model_pos = (mutant.get_position() * 3) # Get the end of the codon
    codon = model.get_sequence()[model_pos - 3:model_pos] # Get the codon
    #   Check to see if the codon is valid
    if codon in PROTEIN_DICTIONARY[mutant.get_reference()]:
        valid = True
    else:
        valid = False
    #   Find a list of all potential SNP positions
    snp_pos = [find_snp_position(codon, alt) for alt in PROTEIN_DICTIONARY[mutant.get_mutant()]]
    #   Calculate the genomic position of the mutation
    dna_pos = model_pos + annotation_list[0].get_start() # Start out with our first CDS
    try:
        for index, cds in enumerate(annotation_list): # For every CDS defined
            if dna_pos > cds.get_end(): # If our DNA position is outside of this CDS
                dna_pos -= cds.get_end() # Subtract the CDS from it
                dna_pos += annotation_list[index + 1].get_start() # Add the remainder to the next CDS
            else: # If not
                break # We have a position! Exit out of the loop
    except IndexError:
        print(mutant, "is beyond the CDS for", model, file=sys.stderr)
        valid = False
    # for cds in annotation_list: # For every CDS defined
    #     dna_pos += cds.get_start()
    #     region = cds.get_end() - cds.get_start()
    #     if dna_pos > cds.get_end():
    #         dna_pos -= region
    #     else:
    #         break
    #     region = cds.get_end() - cds.get_start() # Define a region for the current CDS
    #     if (dna_pos - region) > 0: # If our SNP lies outside of the region defined
    #         dna_pos -= region # Remove the region from our calculate
    #     else: # Otherwise, our SNP is still in this region
    #         dna_pos += cds.get_start() # Calculate the genomic position
    #         break # Break our loop
    #   Create our SNPs
    snp_list = []
    for index in range(len(snp_pos)):
        if snp_pos[index] == -8:
            continue
        snpid = mutant.get_protein() # Use the protein as the SNP ID
        contig = model.get_contig() # Get the contig informaiton from our gene model
        position = dna_pos + snp_pos[index] - 3 # Add the base position, then subtract 3, and avoid weird negation stuff
        if valid:
            reference = codon[snp_pos[index]] # Reference is found in our codon
            alternate = PROTEIN_DICTIONARY[mutant.get_mutant()][index][snp_pos[index]] # Alternate is from the protein dictionary
        else:
            reference = 'N'
            alternate = 'N'
        s = SNP(snpid, contig, position, reference, alternate) # Assemble into SNP
        snp_list.append(s) # Add to our list
    return snp_list


#   Make an argument parser
def make_argument_parser():
    """Set and parse the arugments"""
    parser = argparse.ArgumentParser(
        add_help=True
    )
    parser.add_argument(
        '-g',
        '--gff',
        dest='gff',
        type=str,
        required=True,
        default=None,
        metavar='GFF FILE',
        help="Path to GFF file"
    )
    parser.add_argument(
        '-m',
        '--gene-model',
        dest='model',
        type=str,
        required=True,
        default=None,
        metavar='GENE MODEL',
        help="Path to gene model in FASTA format"
    )
    parser.add_argument(
        '-l',
        '--lookup',
        dest='lookup',
        type=str,
        required=True,
        default=None,
        metavar='SWISSPROT LOOKUP',
        help="Path to list of amino acids. This should have four columns: Protein, Wild type, Mutant, and Position"
    )
    parser.add_argument(
        '-o',
        '--outfile',
        dest='outfile',
        type=str,
        required=False,
        default='output.vcf',
        metavar='OUTPUT NAME',
        help="Name of output file, defaults to 'output.vcf'"
    )
    return parser


#   Driver function
def main():
    header = '##fileformat=VCFv4.2\n##INFO<ID=s,Number=1,Type=Flag,Description="Variant is predicted from Gene Model">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args())
    try:
        #   Read in the gene model as a dictionary
        print("Reading " + args['model'], file=sys.stderr)
        gene_model = open(args['model']).read() # Read in the entire file as a string
        model_list = gene_model.split('>') # Split by the FASTA header
        model_dict = dict(map(parse_gene, model_list[1:])) # Parse the gene model, ignoring the first entry which is ''
        #   Read in the GFF file as a dictionary
        ann_dict = {}
        with open(args['gff'], 'r') as g:
            print("Reading " + args['gff'], file=sys.stderr)
            for line in g:
                if not 'CDS' in line: # If this line isn't a CDS
                    continue # Skip it
                a = parse_gff(line.strip()) # Parse the GFF line, get an annotation object
                if a.get_geneid() in ann_dict.keys(): # If we already have this protein
                    ann_dict[a.get_geneid()].append(a) # Add the annotation to its entry
                else: # Otherwise
                    ann_dict[a.get_geneid()] = [a] # Create an entry for it
        for geneid in ann_dict: # For every protein defined
            ann_dict[geneid].sort(key = lambda ann : ann.get_start())  # Sort by start position
        #   Read in the amino acid lookup table
        mut_list = []
        with open(args['lookup'], 'r') as l:
            print("Reading " + args['lookup'], file=sys.stderr)
            for line in l: # For every line in our lookup table
                split_line = line.strip().split() # Split the line into a list
                #   Create a mutation object
                try:
                    mut = AminoAcidMutant(split_line[0], int(split_line[3]), split_line[1], split_line[2])
                    mut_list.append(mut) # Append the mutation to our list
                except IndexError:
                    print("Skipping line '" + line.strip() + "'", file=sys.stderr)
    except FileNotFoundError as error:
        sys.exit("Failed to find " + error.filename) # Exit with error
    #   Start finding SNPs
    snp_list = [] # Create a holder list
    masked_list = [] # Create a list for masked
    failed = []
    print("Finding SNPs", file=sys.stderr)
    for mutant in mut_list:
        try:
            snps = get_snps(mutant, model_dict[mutant.get_protein()], ann_dict[mutant.get_protein()])
            for snp in snps:
                if snp.check_masked():
                    masked_list.append(snp)
                else:
                    snp_list.append(snp)
        except KeyError:
            failed.append(mutant)
    #   Write out our output file
    with open(args['outfile'], 'w') as out:
        print("Creating VCF file at " + args['outfile'], file=sys.stderr)
        out.write(header)
        for snp in snp_list:
            out.write(snp.format_vcf())
            out.write('\n')
    #   Write out our masked SNPs
    if len(masked_list) > 0:
        masked_name = args['outfile'][::-1].split('.')[1][::-1] + '_masked.vcf'
        print("Writing " + str(len(masked_list)) + " masked SNPs to " + masked_name, file=sys.stderr)
        with open(masked_name, 'w') as masked:
            masked.write(header)
            for snp in masked_list:
                masked.write(snp.format_vcf())
                masked.write('\n')
    #   Write out our failed-to-find mutants
    if len(failed) > 0:
        if os.path.dirname(args['outfile']) == '':
            outdir = os.getcwd()
        else:
            outdir = os.path.dirname(args['outfile'])
        failed_name = outdir + '/failed_mutants.log'
        print("Writing " + str(len(failed)) + " variants that could not be found to " + failed_name, file=sys.stderr)
        with open(failed_name, 'w') as f:
            f.write("Protein\tReference\tAlternate\tPosition\n")
            for mutant in failed:
                f.write(mutant.format_fail())
                f.write('\n')


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
"""
Week 1 & 2 Gene Finder Code
Computes average length of ORF segments in randomly generated DNA
and returns animo acids for given DNA's ORFs segments if they
are long enough
@author: Siena Okuno (sokuno222)
"""

import random
import math
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq



def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide=='A':
        return 'T'
    if nucleotide=='T':
        return 'A'
    if nucleotide=='C':
        return 'G'
    if nucleotide=='G':
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    i=0
    reverse=''
    while i<len(dna):
        letter=dna[i:i+1]
        reverse= get_complement(letter) + reverse
        i=i+1
    return reverse


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    #stop=TAG, TAA, TGA
    a=0
    stoppedDNA=''
    while a<=len(dna)/3:
        if dna[3*a:3+3*a]=='TAG' or dna[3*a:3+3*a]=='TAA' or dna[3*a:3+3*a]=='TGA':
            stoppedDNA=dna[:3*(a)]
            return stoppedDNA
        a=a+1
    if stoppedDNA=='':
        return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    #start=ATG
    #stop=TAG, TAA, TGA
    b=0
    readMe = []
    while b<=len(dna)/3:
        if dna[3*b:3+3*b]=='ATG':
            readMe.append(rest_of_ORF(dna[3*b:]))
            addMe = len(rest_of_ORF(dna[3*b:]))/3
            b=b+(int)(addMe)
        else:
            b=b+1
    return readMe


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    i=0
    allORFs = []
    while i<3:
        allORFs = allORFs+(find_all_ORFs_oneframe(dna[i:]))
        i=i+1
    return allORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    part1=find_all_ORFs(dna)
    part2=find_all_ORFs(get_reverse_complement(dna))
    finalList=part1+part2
    return finalList



#week 1 finished



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATTTAATTCTCCAGGGAATGGGAGAC")
    'ATGGGAGAC'
    """
    # TODO: implement this
    allCombos=find_all_ORFs_both_strands(dna)
    i=0
    c=0
    allORFcombos = sorted(allCombos, key=len, reverse=True)
    while i<len(allORFcombos):
        if c+1==len(allORFcombos):
            pass
        elif len(allORFcombos[c]) > len(allORFcombos[c+1]):
            del allORFcombos[c+1]
        elif len(allORFcombos[c]) < len(allORFcombos[c+1]):
            del allORFcombos[c]
        elif len(allORFcombos[c]) == len(allORFcombos[c+1]):
            c=c+1
        else:
            return "Something went wrong..."
        i=i+1
    return allORFcombos[0]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    dnaList=[]
    endingList=[]
    ORFnoncodingList=[]
    for i in range(num_trials):
        dnaList=list(dna)
        random.shuffle(dnaList)
        dnaScrambled = ''.join(dnaList)
        endingList.append(longest_ORF(dnaScrambled))
    ORFnoncodingMax=max(endingList, key=len)
    if len(ORFnoncodingMax)>0:
        return len(ORFnoncodingMax)
    else:
        return 0


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    from amino_acids import aa_table
    codingStrand=''
    aminoLetter=''
    for i in range((math.floor(len(dna)/3))):
        aminoLetter=aminoLetter+aa_table[dna[i*3:i*3+3]]
    return aminoLetter


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFsList = find_all_ORFs_both_strands(dna)
    FinalORFList=[]
    FinalList=[]
    acidString=''
    for i in range(len(ORFsList)):
        if (len(ORFsList[i])>=threshold):
            FinalORFList.append(ORFsList[i])
    for i2 in range(len(FinalORFList)):
        acidString=''
        for i3 in range(math.floor(len(FinalORFList[i2])/3)):
            acidString+=coding_strand_to_AA(FinalORFList[i2][i3*3:i3*3+3])
        FinalList.append(acidString)
    return FinalList

from load import load_seq
dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))
print('done')

#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

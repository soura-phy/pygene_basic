"""
Python module pygene created by Souradeep Das (github.com/souradeepdas-iisc/)

Initially the gene_primitive.py module created in 2020. This is a modified version designed using OOPs.

This is a DNA:
    
    5' - A-T-G-C-A-T-G-C-..... - 3'
    3' - T-A-C-G-T-A-C-G-..... - 5'
    
            ^ Y
            |
            |
            |______ X
           /
          /
         / Z
  
"""

dict_AA = {'phe': ['Phe', 'UUU', 'UUC'], 'leu': ['Leu', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
           'ile': ['Ile', 'AUU', 'AUC', 'AUA'], 'met': ['Met', 'AUG'], 'val': ['Val', 'GUU', 'GUC', 'GUA', 'GUG'],
           'ser': ['Ser', 'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'pro': ['Pro', 'CCU', 'CCC', 'CCA', 'CCG'],
           'thr': ['Thr', 'ACU', 'ACC', 'ACA', 'ACG'], 'ala': ['Ala', 'GCU', 'GCC', 'GCA', 'GCG'],
           'tyr': ['Tyr', 'UAU', 'UAC'], 'his': ['His', 'CAU', 'CAC'], 'gln': ['Gln', 'CAA', 'CAG'],
           'asn': ['Asn', 'AAU', 'AAC'], 'lys': ['Lys', 'AAA', 'AAG'], 'asp': ['Asp', 'GAU', 'GAC'],
           'glu': ['Glu', 'GAA', 'GAG'], 'cys': ['Cys', 'UGU', 'UGC'], 'trp': ['Trp', 'UGG'],
           'arg': ['Arg', 'CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 
           'gly': ['Gly', 'GGU', 'GGC', 'GGA', 'GGG'], 'stop': ['UAA', 'UAG', 'UGA']}

dict_DNA={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
dict_RNA={'A':'U', 'U':'A', 'C':'G', 'G':'C'}

DNA_symbols=['A', 'T', 'C', 'G', '-', ' ']

dict_seq={'D': dict_DNA, 'R':dict_RNA}

dict_pol_53={True: ['5','3'], False: ['3','5']}

def separate_str(seq, sep=''):
    str1=seq[0]
    for i in seq[1:]:
        str1+=sep+i
    return str1


class Seq:
    # Default data type is string
    def __init__(self, seq, form='D'):
        self.seq = seq.upper()
        self.form = form.upper()
        pass

    def convert_list(self):
        '''

        Returns
        -------
        List
            Converts any string into a list. Analogous to list().

        '''
        seq1 = self.seq
        return list(seq1)

    def reverse(self):
        '''
        
        Returns
        -------
        String
            Converts any string to its reverse string.

        '''
        seq1 = self.seq
        seq_rev = ''
        for i in seq1:
            seq_rev = i + seq_rev
        return Seq(seq_rev, self.form)

    def print_seq(self):
        print(separate_str(self.seq))
        pass

    def get_form(self):
        return self.form

    def get_seq(self):
        return self.seq
    
    def complement(self):
        '''
        

        Returns
        -------
        Seq object
            Returns the complement sequence of the input sequence.
            The polarity of output sequence is 3'-5' for input sequence of polarity 5'-3'.

        '''
        dict1=dict_seq[self.form]
        seq1=self.seq
        seq2=''
        for char in seq1:
            seq2+=dict1[char]
        return Seq(seq2, self.form)
    
    def comp_reversed(self):
        return self.complement().reverse()


class DNAseq(Seq):
    def __init__(self, seq, pol_53=True):
        #pol_53 denotes whether the DNA is in 5'-3' polarity or not
        super().__init__(seq, 'D')
        self.compseq = self.complement().seq
        self.pol_53=pol_53
        pass
    
    def is_valid(self):
        for char in self.seq:
            if char not in DNA_symbols:
                return 0
        return 1
    
    def flip_x(self):
        return DNAseq(self.compseq, pol_53= not self.pol_53)
    
    def flip_y(self):
        Seq1=Seq(self.seq)
        seqn1=Seq1.reverse().seq
        return DNAseq(seqn1, pol_53= not self.pol_53)

    def flip_z(self):
        Seq1=Seq(self.compseq)
        seqn1=Seq1.reverse().seq
        return DNAseq(seqn1, pol_53= self.pol_53)
    
    def print_seq(self, sep=''):
        str1=separate_str(self.seq, sep)
        str2=separate_str(self.compseq, sep)
        TL=dict_pol_53[self.pol_53][0]#number on top-left
        TR=dict_pol_53[self.pol_53][1]#number on top-right
        
        print(TL+"' - "+str1+" - "+TR+"'")
        print(TR+"' - "+str2+" - "+TL+"'")
        pass
        
    def analyze(self):
        pass

#%%
class RNAseq(Seq):
    def __init__(self, seq):
        self.seq = seq
        self.form = 'R'


def compDNA(strand):
    """Returns the complementary DNA strand corresponding to any DNA or RNA strand.

    """
    z = ''
    for i in strand:
        if i == 'A':
            z = z + 'T'
        elif i == 'T' or i == 'U':
            z = z + 'A'
        elif i == 'G':
            z = z + 'C'
        elif i == 'C':
            z = z + 'G'
        else:
            print("\t\t\t\t\t\t! Error in code")
            break
    return z


def compRNA(strand):
    """Returns the complementary RNA strand corresponding to any DNA or RNA strand.

    """
    z = ''
    for i in strand:
        if i == 'A':
            z = z + 'U'
        elif i == 'T' or i == 'U':
            z = z + 'A'
        elif i == 'G':
            z = z + 'C'
        elif i == 'C':
            z = z + 'G'
        else:
            print("\t\t\t\t\t\t! Error in code")
            break
    return z


def transcribe(DNAstrand):
    """Transcribes any template strand or coding strand of a transcription unit into its corresponding hnRNA strand.

    Any of the two strands in any direction can be given as the input.
    The module verifies which strand you have inserted and in which polarity, and automatically prepares the 3 strands:
        1. The Coding Strand in 5'->3' direction,
        2. The Template Strand in 3'->5' direction, and
        3. The hnRNA Strand in 5'->3' direction."""
    global temp, code, anticode
    print('ENTRY PARAMETERS')
    sure = 'N'

    def fp1():
        global p1, p10, p11
        p1 = input("Please enter\n\tT if you enter Template Strand\n\tC if you enter Coding Strand\n\t\t")
        if p1 == 'T':
            p10 = 'Template Strand'
            p11 = 'Coding Strand'
        elif p1 == 'C':
            p11 = 'Template Strand'
            p10 = 'Coding Strand'
        else:
            print("This input is not allowed")
            fp1()

    def fp2():
        global p2, p20, p21
        p2 = input("Which is the direction in your strand?\n\tEnter 53 if 5\'->3\'\n\tEnter 35 if 3\'->5\'\n\t\t")
        if p2 == '53':
            p20 = '5\'->3\''
            p21 = '3\'->5\''
        elif p2 == '35':
            p21 = '5\'->3\''
            p20 = '3\'->5\''
        else:
            print("This input is not allowed")
            fp2()

    while sure == 'N':
        fp1()
        fp2()
        print("Are you sure to continue your input as ", p10, " with direction ", p20, " ?")
        sure = input("Please enter Y if yes, N if no.\t")
    rawcode = compDNA(DNAstrand) if p1 == 'T' else DNAstrand
    if p1 == 'T' and p2 == '53':
        code = rev(rawcode)
    elif p1 == 'T' and p2 == '35':
        code = rawcode
    elif p1 == 'C' and p2 == '53':
        code = rawcode
    else:
        code = rev(rawcode)
    temp = compDNA(code)
    anticode = compRNA(temp)
    print("Your coding DNA strand is :        5\'-", code, "-3\'\n")
    print("Your template DNA strand is :     3\'-", temp, "-5\'\n")
    print("Your transcribed RNA strand is :  5\'-", anticode, "-3\'\n")


def translate(RNAstrand):
    """Translates a valid mRNA sequence into its polypeptide sequence, dictated by the genetic code.

    Inputs a string of A,U,G,C and displays the list of codons as interpreted by ribosome.
    Codons list always starts from the first encountered AUG(Start codon).
    All codons are converted into corresponding amino acids.
    Translation terminates upon reaching a stop codon, any one of UAA, UAG, UGA"""
    global peptide, codon
    peptide, codon = [], []
    import string
    a = RNAstrand.find('AUG')
    while a <= len(RNAstrand) - 3:
        if a == -1:
            print("Absence of start codon AUG.\Sorry, the code you entered cannot be translated.")
            break
        else:
            codon.append(RNAstrand[a:a + 3])
            a += 3
    peptide = runcode(codon)
    peptidestring = '-'.join(peptide)
    print("\tThe codons encountered are :\n\t\t", codon, '\n')
    print("\tThe polypeptide chain translated is :\n\t\t", peptidestring, '\n')


def runcode(list):
    """Converts a list of valid codons into a list of their corresponding amino acids.

    Returns the list of amino acids in proper order."""
    global y
    y = []
    for v in list:
        c = 0
        for i in dict_AA.keys():
            if v in dict_AA['stop']:
                c = 1
                break
            elif v in i:
                y.append(i[0])
        if c == 1:
            break
    return y

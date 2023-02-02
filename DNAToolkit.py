from collections import Counter

Nucleotides=['A','C','G','T']
DNA_ReverseComplement = {"A":"T","T":"A","C":"G","G":"C"}

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

#Function to check the sequence to make sure it is a DNA string
def validateSeq(dna_seq):
    tmpseq=dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq
#Function to count nucleotides frequency
def countNucFrequency(seq):
    tmpFreqDict = {"A":0,"C":0,"G":0,"T":0}
    for nuc in seq:
        tmpFreqDict[nuc] +=1
    return tmpFreqDict
#Function to calculate GC content
def GC_content(seq):
    return round((seq.count('C')+ seq.count('G'))/len(seq)*100)
#Function to calculate GC content in subsections
def GC_Content_Subset(seq, k=20):
    res=[]
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(GC_content(subseq))
    return res
#Transcription: creatimg RNA from DNA
def transcription(seq):
    return seq.replace("T","U")
#Create reverse complement
def reverse_complement(seq):
    mapping =str.maketrans('ATCG','TAGC')
    return seq.translate(mapping)[::-1]

#Translates a DNA sequence into an aminoacid sequence
def translate_seq(seq, init_pos=0):
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]

#Provides the frequency of each codon encoding a given aminoacid in a DNA sequence
def codon_usage(seq, aminoacid):
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict

#Generate the six reading frames of a DNA sequence, including 3 reverse complement
def gen_reading_frames(seq):
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq,1))
    frames.append(translate_seq(seq,2))
    frames.append(translate_seq(reverse_complement(seq),0))
    frames.append(translate_seq(reverse_complement(seq),1))
    frames.append(translate_seq(reverse_complement(seq),2))
    return frames
def proteins_from_rf(aa_seq):
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
                # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins
def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)
        res = []
        for rf in rfs:
            prots = proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res
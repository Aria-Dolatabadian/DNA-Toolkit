from DNAToolkit import *
# randDNAStr='ATGACAGATGAAGCACGA'
# randDNAStr= input('Insert your sequence:')
# print(validateSeq(rndDNAStr))
import random
# #Create a random DNA seq for testing
randDNAStr=''.join([random.choice(Nucleotides)
                    for nuc in range(50)])
# print(validateSeq(randDNAStr))
# print(countNucFrequency(randDNAStr))
# print(transcription(randDNAStr))
from Utilities import colored
DNAStr=validateSeq(randDNAStr)

print(f'\nSequence: {colored(DNAStr)}\n')
print(f'[1] Sequence length: {len(DNAStr)}\n')
print(colored(f'[2] Nucleotides Frequency: {countNucFrequency(DNAStr)}\n'))
print(colored(f'[3] GC content (%): {GC_content(randDNAStr)}\n'))
print(colored(f'[4] GC content in subsections k=5: {GC_Content_Subset(randDNAStr, k=5)}\n'))
print(f'[5] DNA-->RNA Transcription: {colored(transcription(DNAStr))}\n')
print(f"[6] DNA String + Complement + Reverse complement:\n 5'...{colored(DNAStr)}...3'")
print(f"      {''.join(['|' for c in range(len(DNAStr))])}")
print(f" 3'...{colored(reverse_complement(DNAStr)[::-1])}...5' [Complement]")
print(f" 5'...{colored(reverse_complement(DNAStr))}...3' [Reverse Complement]\n")
print(f"[7] Amino acids sequence from DNA: {translate_seq(DNAStr,0)}\n")
print(f"[8] Codon frequency sequence (L): {codon_usage(DNAStr,'L')}\n")
print(f"[8] Codon frequency sequence (R): {codon_usage(DNAStr,'R')}\n")
print('[9] Reading frames:\n')
for frame in gen_reading_frames(DNAStr):
    print(frame)

print('\n[10] Allproteins in 6 open reading frames:\n')
for prot in all_proteins_from_orfs(DNAStr,0,0, True):
    print(f'{prot}')
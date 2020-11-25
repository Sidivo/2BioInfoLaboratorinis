



# https://en.wikipedia.org/wiki/FASTQ_format - info apie formatus.

from collections import Counter
from bioinfokit.analys import fastq
import math


# atrandame, kokia koduote parasyta
def check_encoding(FastqEncoding):
    for c in ASCIIsymbols:
        if c not in FastqEncoding:
            return False
    return True



# atrandame C ir G simboliu  santyki FASTQ failo biologineje sekoje.
def C_G_ratio_finder(FASTQ_sequence):
    FASTQ_sequence_length = float(len(FASTQ_sequence))
    C_NucleotideCount = FASTQ_sequence.count('C')
    G_NucleotideCount = FASTQ_sequence.count('G')
    ratio = (C_NucleotideCount + G_NucleotideCount) / FASTQ_sequence_length
    ratio = round_ceil_function(ratio)
    return ratio



# per paskaita sakete, kad reikia iki 2 skaiciu po kableliu apvalint, tai funckija ta ir daro.
def round_ceil_function(ratio):
    return(math.ceil(ratio * 100) / 100)



# konvertuojame is listo i stringa, kad apskaiciuoti ivairias reiksmes, nes lengviau dirbt negu su list.
def list_to_string(list):
    return "".join(list)


if __name__ == '__main__':




    # https://reneshbedre.github.io/blog/filereaders.html -  buvo pasinaudota siuo pavyzdziu ir biblioteka nuskaitymui FASTQ failo.
    fastq_iter = fastq.fastq_reader(file='C:/Users/MSI/PycharmProjects/pythonProject3/reads_for_analysis.fastq')

    qualitySequenceList = []
    C_G_SymbolsRatioList = []


    for record in fastq_iter:
        #header_1 - @SEQ_ID
        #header_2 - +
        header_1, sequence, header_2, qualitySequence = record

        # Sujungiame visas įverčio sekas.
        qualitySequenceList.append(qualitySequence)

        # Iskart tame paciame cikle surandame ir C G simboliu santykius.
        C_G_SymbolsRatioList.append(C_G_ratio_finder(sequence))
        # print(C_G_SymbolsRatioList)

    # print(qualitySequenceList)
    # nurodomi santykiai ir ju kiekis.
    allRatio = Counter(C_G_SymbolsRatioList)
    print(allRatio)


    ASCIIsymbols = list(set(list_to_string(qualitySequenceList)))
    print(ASCIIsymbols)


    #ASCII simboliai koduotei atrasti
    Sanger = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"""
    Solexa = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
    Illumina13 = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
    Illumina15 = "CDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi"
    Illumina18 = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"""


    if check_encoding(Sanger):
        print("FASTQ file given information uses Sanger encoding")
    if check_encoding(Solexa):
        print("FASTQ file given information uses Solexa encoding")
    if check_encoding(Illumina13):
        print("FASTQ file given information uses Illumina 1.3+ encoding")
    if check_encoding(Illumina15):
        print("FASTQ file given information uses Illumina 1.5+ encoding")
    if check_encoding(Illumina18):
        print("FASTQ file given information uses Illumina 1.8+ encoding")
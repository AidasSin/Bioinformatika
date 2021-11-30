import re
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

kodavimai = { 
    'Sanger Phred+33': (0, 40),
    'Solexa Solexa+64': (-5, 40),
    'Illumina 1.3+ Phred+64': (0, 40),
    'Illumina 1.5+ Phred+64': (3, 41),
    'Illumina 1.8+ Phred+33': (0, 41)
}

def nuskaityti_fastq_koduotes_skaicius(fastq_seka):

    skaiciai = []
    for record in SeqIO.parse(fastq_seka, "fastq"):
        skaiciai.append(record.letter_annotations["phred_quality"])
    
    return skaiciai

def gauti_skaiciu_intervala(skaiciai):

    minimalus = min(skaiciai)
    maximalus = max(skaiciai)

    return(minimalus, maximalus)

def gauti_skaiciu_intervalus(visi_skaiciai):

    intervalu_dictionary_masyvas = []
    for i in range(len(visi_skaiciai)):
        intervalu_dictionary_masyvas.append(gauti_skaiciu_intervala(visi_skaiciai[i]))
    
    return intervalu_dictionary_masyvas

def gauti_visos_sekos_intervala(dictionary):

    minimalus = min(min(dictionary))
    maximalus = max(max(dictionary))

    return(minimalus, maximalus)

def nustatyti_kokybės_kodavimą(rangeas):

    if(rangeas[0] >= 0 and rangeas[1] <= 40):
        return "Sanger Phred+33 arba Illumina 1.3+ Phred+64"
    if(rangeas[0] < 0):
        return "Solexa Solexa+64"
    if(rangeas[0] >= 3 and rangeas[1] <= (41)):
        return "Illumina 1.5+ Phred+64"
    if(rangeas[1] == 41):
        return "Illumina 1.8+ Phred+33"

def gauti_sekas(failas):

    sekos = []
    for record in SeqIO.parse(failas, "fastq"):
     sekos.append(str(record.seq))
    return sekos

def gauti_seku_id(failas):

    sekuId = []
    for record in SeqIO.parse(failas, "fastq"):
     sekuId.append(str(record.id))
    return sekuId


def ketvirto_a(failas):

    skaiciai = nuskaityti_fastq_koduotes_skaicius(fastq_seka=failas)
    dictionary = gauti_skaiciu_intervalus(skaiciai)
    rangeas = gauti_visos_sekos_intervala(dictionary)

    return(nustatyti_kokybės_kodavimą(rangeas))
    
def gauti_GC_dazni_sekoje(seka):

    GCskaicius = 0
    sekosIlgis = len(seka)
    for i in range(len(seka)):
        if seka[i] == 'G' or seka[i] == 'C':
            GCskaicius += 1
    return (GCskaicius * 100) / sekosIlgis

def gauti_GC_daznius(sekos):
    
    dazniai = []
    for i in range(len(sekos)):
        dazniai.append(gauti_GC_dazni_sekoje(sekos[i]))
    return dazniai


# Iš pradžių čia neteisingai mąsčiau - ėmiau kiekvieno read'o dažnį ir žiūrėjau, kuriuose didžiausias dažnis pavieniui.
# def gauti_grafika(sekos):

#     dazniai = gauti_GC_daznius(sekos)

#     plt.plot(dazniai, list(range(1, len(sekos)+1)), 'ro')
#     plt.xlabel("C/G Nukleotidų dalis (%)")
#     plt.ylabel("Read'ų skaičius")
#   # plt.show()

def gauti_grafika(sekos):
    dazniai = gauti_GC_daznius(sekos)

    GCkiekiai = []
    intervalai = []
    i = 0
    while i <= 100:
        intervalai.append(i)
        i += 4

    for i in range(len(intervalai)-1):
        # print(intervalai[i], intervalai[i+1])
        GCkiekiai.append(sum(intervalai[i] < x < intervalai[i+1] for x in dazniai))

    plt.plot(intervalai[1:], GCkiekiai)
    plt.xlabel("C/G Nukleotidų dalis (%)")
    plt.ylabel("Read'ų skaičius")
    #plt.show()

def top5_36(dazniai):
    procentai36 = []
    top5Procentai36 = []
    for i in range(len(dazniai)):
        if(dazniai[i] >= 36 and dazniai[i] <= 37):
            procentai36.append(i)

    reverseProcentai36 = procentai36[::-1]
    i=0
    while len(top5Procentai36) != 5:
        top5Procentai36.append(reverseProcentai36[i])
        i += 1
    return top5Procentai36

def top5_56(dazniai):
    procentai56 = []
    top5Procentai56 = []
    for i in range(len(dazniai)):
        if(dazniai[i] >= 56 and dazniai[i] <= 57):
            procentai56.append(i)

    reverseProcentai36 = procentai56[::-1]
    i=0
    while len(top5Procentai56) != 5:
        top5Procentai56.append(reverseProcentai36[i])
        i += 1
    return top5Procentai56

def top5_72(dazniai):
    procentai72 = []
    top5Procentai72 = []
    for i in range(len(dazniai)):
        if(dazniai[i] >= 72 and dazniai[i] <= 73):
            procentai72.append(i)

    reverseProcentai36 = procentai72[::-1]
    i=0
    while len(top5Procentai72) != 5:
        top5Procentai72.append(reverseProcentai36[i])
        i += 1
    return top5Procentai72

print(ketvirto_a("reads_for_analysis.fastq"))

gauti_grafika(sekos = gauti_sekas("reads_for_analysis.fastq"))

sekos = gauti_sekas("reads_for_analysis.fastq")
dazniai = gauti_GC_daznius(sekos)

sekuId = gauti_seku_id("reads_for_analysis.fastq")

penki36 = top5_36(dazniai)
penki56 = top5_56(dazniai)
penki72 = top5_72(dazniai)

print("\nTop 5 36%")
for i in range(len(penki36)):
    print(dazniai[penki36[i]])
    print(sekuId[penki36[i]])
    print(sekos[penki36[i]])

print("\nTop 5 56%")
for i in range(len(penki56)):
    print(dazniai[penki56[i]])
    print(sekuId[penki56[i]])
    print(sekos[penki56[i]])

print("\nTop 5 72%")
for i in range(len(penki72)):
    print(dazniai[penki72[i]])
    print(sekuId[penki72[i]])
    print(sekos[penki72[i]])


# print(dazniai[18520])
# print(dazniai[20095])
#print(dazniai[18520])
#print(sekos[18520])
#print(dazniai[6047])
#print(sekos[6047])
#print(dazniai[10558])
#print(sekos[10558])

# # top 4 vietos: [18520 20095 11479  6047]
# sekos = gauti_sekas("reads_for_analysis.fastq")
# dazniai = gauti_GC_daznius(sekos)

# npSekos = np.array(dazniai)
# # #top4 = np.argpartition(npSekos, -4)[-4:]

# top10 = np.argpartition(npSekos, -20)[-20:]
# print(top10)




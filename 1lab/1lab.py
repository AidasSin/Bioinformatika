import itertools
from re import S
import re
from Bio import SeqIO
from itertools import product
import pprint


# AUG - start codon -> ATG
# There are 3 STOP codons in the genetic code - UAG, UAA, and UGA. -> TAG, TAA, TGA


nukleorugstys = ["A", "T", "G", "C"]

visi_kodonai = [''.join(i) for i in itertools.product(nukleorugstys, repeat = 3)]
visi_dikodonai = [''.join(i) for i in itertools.product(visi_kodonai, repeat = 2)]


def nuskaityti_genoma(genomas):

    for mamalian1_seq_record in SeqIO.parse(genomas, "fasta"):
        (seka) = (mamalian1_seq_record.seq)
    
    return seka




def gauti_remeli(seka, skaicius):
    return [str(seka)[i:i+3] for i in range(skaicius, len(seka), 3)]

def gauti_remelius(seka): # prideti forus cia
    
    remeliai = []
    reverse_komplementari_seka = seka.reverse_complement()
    #pirmas_remelis  
    remeliai.append(gauti_remeli(seka=seka, skaicius=0))
    #antras_remelis
    remeliai.append(gauti_remeli(seka=seka, skaicius=1))
    #trecias_remelis = 
    remeliai.append(gauti_remeli(seka=seka, skaicius=2))
  #  pirmas_reverse_komplementarios_remelis = 
    remeliai.append(gauti_remeli(seka=reverse_komplementari_seka, skaicius=0))
  #  antras_reverse_komplementarios_remelis = 
    remeliai.append(gauti_remeli(seka=reverse_komplementari_seka, skaicius=1))
  # trecias_reverse_komplementarios_remelis = 
    remeliai.append(gauti_remeli(seka=reverse_komplementari_seka, skaicius=2))

    return remeliai


def gauti_start_stop_poras(remelis):
    
    remelio_poros = []
    remelio_ilgis = len(remelis)
    ilgio_check = 0

    while ilgio_check < remelio_ilgis:

        start_kodono_pozicija = ilgio_check
        while remelis[start_kodono_pozicija] != "ATG" and start_kodono_pozicija < remelio_ilgis-1:
            start_kodono_pozicija += 1

        seka_nuo_start_kodono = remelis[start_kodono_pozicija:]
        sekos_indeksai = len(seka_nuo_start_kodono) - 1

       # print(seka_nuo_start_kodono)

        stop_kodono_pozicija = 0
        # čia jau einu pro kitą sąrašą, tai negaliu naudoti bendro indexo.
        #all ([seka_nuo_start_kodono[stop_kodono_pozicija] != "TGA", seka_nuo_start_kodono[stop_kodono_pozicija] != "TAA", seka_nuo_start_kodono[stop_kodono_pozicija] != "TAG", stop_kodono_pozicija < remelio_ilgis ])

        while all ([seka_nuo_start_kodono[stop_kodono_pozicija] != 'TGA', seka_nuo_start_kodono[stop_kodono_pozicija] != 'TAA', seka_nuo_start_kodono[stop_kodono_pozicija] != 'TAG', stop_kodono_pozicija < sekos_indeksai]):  #(seka_nuo_start_kodono[stop_kodono_pozicija]) != "TGA" and (seka_nuo_start_kodono[stop_kodono_pozicija]) != "TAA" and (seka_nuo_start_kodono[stop_kodono_pozicija]) != "TAG" and stop_kodono_pozicija < remelio_ilgis:
            stop_kodono_pozicija += 1
        
        start_stop_poros = seka_nuo_start_kodono[:stop_kodono_pozicija+1]
        if start_stop_poros[-1] != 'TGA' and start_stop_poros[-1] != 'TAA' and start_stop_poros[-1] != 'TAG':
            break

        gauta_pora = ''.join(start_stop_poros)

        #print("\nstop_kodono_pozicija = " + str(stop_kodono_pozicija) + ", start kodono pozicija = " + str(start_kodono_pozicija))
        #print("start_stop_pora:\n")
        #print(start_stop_poros)
        ilgio_check = start_kodono_pozicija + stop_kodono_pozicija + 1
        remelio_poros.append(gauta_pora)
        
    return remelio_poros


def gauti_ORFus(remelis):
    ORFai = []
    remelio_ilgis = len(remelis)

    ORFo_pradzia = 0
    while ORFo_pradzia < remelio_ilgis:
        #paskutinio stop kodono pozicija arba remelio pradzia:
        #gal nereikia? ORFo_ieskojimo_pradzios_pozicija = ilgio_check
        while all ([remelis[ORFo_pradzia] != 'TGA', remelis[ORFo_pradzia] != 'TAA', remelis[ORFo_pradzia] != 'TAG', ORFo_pradzia < remelio_ilgis - 1]):
            ORFo_pradzia += 1

        seka_nuo_ORFo_pradzios = remelis[ORFo_pradzia:]

        if seka_nuo_ORFo_pradzios[0] != 'TGA' and seka_nuo_ORFo_pradzios[0] != 'TAA' and seka_nuo_ORFo_pradzios[0] != 'TAG':
            break

        sekos_indeksas = len(seka_nuo_ORFo_pradzios) - 1
        start_kodono_pozicija = 0

        ORFo_pabaiga = 1
        while all ([seka_nuo_ORFo_pradzios[ORFo_pabaiga] != 'TGA', seka_nuo_ORFo_pradzios[ORFo_pabaiga] != 'TAA', seka_nuo_ORFo_pradzios[ORFo_pabaiga] != 'TAG', ORFo_pabaiga < sekos_indeksas]):
            if seka_nuo_ORFo_pradzios[ORFo_pabaiga] == 'ATG':
                start_kodono_pozicija = ORFo_pabaiga
            ORFo_pabaiga += 1
        
        if start_kodono_pozicija == 0:
            ORFo_pradzia += ORFo_pabaiga# + 1
            continue

        ORFas = seka_nuo_ORFo_pradzios[:start_kodono_pozicija+1]
        #print(ORFas)
        ORFo_pradzia += ORFo_pabaiga + 1
        ORFai.append(''.join(ORFas))
        
    return ORFai


def pasalinti_trumpesnius_nei_100(fragmentas):

    ilgos_sekos = []
    for i in range(len(fragmentas)):
        if len(fragmentas[i]) > 100:
            ilgos_sekos.append(fragmentas[i])
    return ilgos_sekos



start_message = "----------begin------------\n\n"

def pirmas_punktas(remeliai):

    remeliu_start_stop_poros = []

    for i in range(len(remeliai)):
        remeliu_start_stop_poros.append(gauti_start_stop_poras(remeliai[i]))

    return remeliu_start_stop_poros


#visos_start_stop_poros = pirmas_punktas(gauti_remelius(seka=seka))

#print(len(visos_start_stop_poros))


def antras_punktas(remeliai):
    remeliu_ORFai = []

    for i in range(len(remeliai)):
        remeliu_ORFai.append(gauti_ORFus(remeliai[i]))
    
    return remeliu_ORFai


def trecias_punktas(fragmentai):
    
    ilgi_fragmentai = []

    for i in range(len(fragmentai)):
        ilgi_fragmentai.append(pasalinti_trumpesnius_nei_100(fragmentai[i]))
    
    return ilgi_fragmentai

def visu_ilgu_seku_sujungimas(ilgos_sekos):
    eilute = ""
    for i in range (len(ilgos_sekos)):
        eilute += ''.join(ilgos_sekos[i])
    return eilute

def kodonu_arba_dikodonu_daznis(koduojancios_sekos, kodonai_arba_dikodonai):
    kodonu_arba_dikodonu_skaicius = dict((x,0) for x in kodonai_arba_dikodonai)
    for w in re.findall(r"\w+",koduojancios_sekos):
        if w in kodonu_arba_dikodonu_skaicius:
            kodonu_arba_dikodonu_skaicius[w] += 1
    return kodonu_arba_dikodonu_skaicius

def gauti_genomo_koduojancias_ilgas_sekas(seka):

    visos_start_stop_poros = trecias_punktas(pirmas_punktas(gauti_remelius(seka=seka)))
    start_stop_poru_bendras = visu_ilgu_seku_sujungimas(ilgos_sekos=visos_start_stop_poros)
    visi_ORFai = trecias_punktas(antras_punktas(gauti_remelius(seka=seka)))
    ORFu_bendras = visu_ilgu_seku_sujungimas(ilgos_sekos=visi_ORFai)
    bendras_abieju = start_stop_poru_bendras + ORFu_bendras

    return bendras_abieju

def gauti_kodonu_dazni_sekoje(seka):

    bendras_abieju = gauti_genomo_koduojancias_ilgas_sekas(seka=seka)
    visu_seku_kodonu_konstravimas = ' '.join([bendras_abieju[i:i+3] for i in range(0, len(bendras_abieju), 3)]) # cia vos ne gaunu jau visu kodonu skaiciu, tai galeciau i lista sudeti ir panaudoti dazniui apskaiciuoti
    kodonu_skaicius = len(visu_seku_kodonu_konstravimas.split(" ")) + 1
   # print(kodonu_skaicius)
    kodonu_pasikartojimai = kodonu_arba_dikodonu_daznis(visu_seku_kodonu_konstravimas,kodonai_arba_dikodonai=visi_kodonai)
   # print(kodonu_pasikartojimai)
    kodonu_skaicius = sum(kodonu_pasikartojimai.values(), 0.0)
   # print(kodonu_skaicius)
    kodonu_daznis = {k: v/ kodonu_skaicius * 100 for k, v in kodonu_pasikartojimai.items()}
    #print(kodonu_daznis)

    return kodonu_daznis

def gauti_dikodonu_dazni_sekoje(seka):

    bendras_abieju = gauti_genomo_koduojancias_ilgas_sekas(seka=seka)
    visu_seku_dikodonu_konstravimas = ' '.join([bendras_abieju[i:i+6] for i in range(0, len(bendras_abieju), 3)])
    #print(visu_seku_dikodonu_konstravimas)
    dikodonu_pasikartojimai = kodonu_arba_dikodonu_daznis(visu_seku_dikodonu_konstravimas,kodonai_arba_dikodonai=visi_dikodonai)
   # print(dikodonu_pasikartojimai)
    dikodonu_skaicius = sum(dikodonu_pasikartojimai.values())
    #print(dikodonu_skaicius)
    dikodonu_daznis = {k: v/ dikodonu_skaicius * 100 for k, v in dikodonu_pasikartojimai.items()}
    #print(dikodonu_daznis)
    return dikodonu_daznis

def gauti_dvieju_seku_atstuma(seka1_dazniai,seka2_dazniai):

    atstumai = []
    for (k,v) in seka1_dazniai.items():
        if (k in seka2_dazniai and seka2_dazniai[k] != 0):
            atstumai.append(seka1_dazniai[k] / seka2_dazniai[k])

    # print(atstumai)
    # print(len(atstumai))
    # print(sum(atstumai))
    atstumas = sum(atstumai) / len(atstumai)
    return atstumas


bacterial1 = nuskaityti_genoma("bacterial1.fasta")
bacterial2 = nuskaityti_genoma("bacterial2.fasta")
bacterial3 = nuskaityti_genoma("bacterial3.fasta")
bacterial4 = nuskaityti_genoma("bacterial4.fasta")

mamalian1 = nuskaityti_genoma("mamalian1.fasta")
mamalian2 = nuskaityti_genoma("mamalian2.fasta")
mamalian3 = nuskaityti_genoma("mamalian3.fasta")
mamalian4 = nuskaityti_genoma("mamalian4.fasta")


bacterial1_kodonu_daznis = gauti_kodonu_dazni_sekoje(bacterial1)
bacterial2_kodonu_daznis = gauti_kodonu_dazni_sekoje(bacterial2)
bacterial3_kodonu_daznis = gauti_kodonu_dazni_sekoje(bacterial3)
bacterial4_kodonu_daznis = gauti_kodonu_dazni_sekoje(bacterial4)

bacterial1_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(bacterial1)
bacterial2_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(bacterial2)
bacterial3_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(bacterial3)
bacterial4_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(bacterial4)

mamalian1_kodonu_daznis = gauti_kodonu_dazni_sekoje(mamalian1)
mamalian2_kodonu_daznis = gauti_kodonu_dazni_sekoje(mamalian2)
mamalian3_kodonu_daznis = gauti_kodonu_dazni_sekoje(mamalian3)
mamalian4_kodonu_daznis = gauti_kodonu_dazni_sekoje(mamalian4)

mamalian1_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(mamalian1)
mamalian2_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(mamalian2)
mamalian3_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(mamalian3)
mamalian4_dikodonu_daznis = gauti_dikodonu_dazni_sekoje(mamalian4)


atstumu_matrica = [[0 for y in range(8)] for x in range(8)]

atstumu_matrica[0][0] = "bacterial1"
atstumu_matrica[1][0] = "bacterial2"
atstumu_matrica[2][0] = "bacterial3"
atstumu_matrica[3][0] = "bacterial4"
atstumu_matrica[4][0] = "mamalian1"
atstumu_matrica[5][0] = "mamalian2"
atstumu_matrica[6][0] = "mamalian3"
atstumu_matrica[7][0] = "mamalian4"

#bacterial1
atstumu_matrica[0][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=bacterial2_kodonu_daznis)
atstumu_matrica[0][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=bacterial3_kodonu_daznis)
atstumu_matrica[0][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=bacterial4_kodonu_daznis)
atstumu_matrica[0][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=mamalian1_kodonu_daznis)
atstumu_matrica[0][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=mamalian2_kodonu_daznis)
atstumu_matrica[0][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=mamalian3_kodonu_daznis)
atstumu_matrica[0][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_kodonu_daznis, seka2_dazniai=mamalian4_kodonu_daznis)

#bacterial2
atstumu_matrica[1][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=bacterial1_kodonu_daznis)
atstumu_matrica[1][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=bacterial3_kodonu_daznis)
atstumu_matrica[1][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=bacterial4_kodonu_daznis)
atstumu_matrica[1][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=mamalian1_kodonu_daznis)
atstumu_matrica[1][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=mamalian2_kodonu_daznis)
atstumu_matrica[1][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=mamalian3_kodonu_daznis)
atstumu_matrica[1][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_kodonu_daznis, seka2_dazniai=mamalian4_kodonu_daznis)

#bacterial3
atstumu_matrica[2][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=bacterial1_kodonu_daznis)
atstumu_matrica[2][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=bacterial2_kodonu_daznis)
atstumu_matrica[2][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=bacterial4_kodonu_daznis)
atstumu_matrica[2][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=mamalian1_kodonu_daznis)
atstumu_matrica[2][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=mamalian2_kodonu_daznis)
atstumu_matrica[2][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=mamalian3_kodonu_daznis)
atstumu_matrica[2][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_kodonu_daznis, seka2_dazniai=mamalian4_kodonu_daznis)

#bacterial4
atstumu_matrica[3][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=bacterial1_kodonu_daznis)
atstumu_matrica[3][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=bacterial2_kodonu_daznis)
atstumu_matrica[3][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=bacterial3_kodonu_daznis)
atstumu_matrica[3][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=mamalian1_kodonu_daznis)
atstumu_matrica[3][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=mamalian2_kodonu_daznis)
atstumu_matrica[3][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=mamalian3_kodonu_daznis)
atstumu_matrica[3][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_kodonu_daznis, seka2_dazniai=mamalian4_kodonu_daznis)

#mamalian1
atstumu_matrica[4][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = bacterial1_kodonu_daznis)
atstumu_matrica[4][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = bacterial2_kodonu_daznis)
atstumu_matrica[4][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = bacterial3_kodonu_daznis)
atstumu_matrica[4][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = bacterial4_kodonu_daznis)
atstumu_matrica[4][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = mamalian2_kodonu_daznis)
atstumu_matrica[4][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = mamalian3_kodonu_daznis)
atstumu_matrica[4][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_kodonu_daznis, seka2_dazniai = mamalian4_kodonu_daznis)

#mamalian2
atstumu_matrica[5][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = bacterial1_kodonu_daznis)
atstumu_matrica[5][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = bacterial2_kodonu_daznis)
atstumu_matrica[5][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = bacterial3_kodonu_daznis)
atstumu_matrica[5][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = bacterial4_kodonu_daznis)
atstumu_matrica[5][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = mamalian1_kodonu_daznis)
atstumu_matrica[5][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = mamalian3_kodonu_daznis)
atstumu_matrica[5][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_kodonu_daznis, seka2_dazniai = mamalian4_kodonu_daznis)

#mamalian3
atstumu_matrica[6][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = bacterial1_kodonu_daznis)
atstumu_matrica[6][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = bacterial2_kodonu_daznis)
atstumu_matrica[6][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = bacterial3_kodonu_daznis)
atstumu_matrica[6][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = bacterial4_kodonu_daznis)
atstumu_matrica[6][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = mamalian1_kodonu_daznis)
atstumu_matrica[6][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = mamalian2_kodonu_daznis)
atstumu_matrica[6][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_kodonu_daznis, seka2_dazniai = mamalian4_kodonu_daznis)

#mamalian4
atstumu_matrica[7][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = bacterial1_kodonu_daznis)
atstumu_matrica[7][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = bacterial2_kodonu_daznis)
atstumu_matrica[7][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = bacterial3_kodonu_daznis)
atstumu_matrica[7][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = bacterial4_kodonu_daznis)
atstumu_matrica[7][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = mamalian1_kodonu_daznis)
atstumu_matrica[7][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = mamalian2_kodonu_daznis)
atstumu_matrica[7][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_kodonu_daznis, seka2_dazniai = mamalian3_kodonu_daznis)

#pprint.pprint(atstumu_matrica)


dikodonu_atstumu_matrica = [[0 for y in range(8)] for x in range(8)]

dikodonu_atstumu_matrica[0][0] = "bacterial1"
dikodonu_atstumu_matrica[1][0] = "bacterial2"
dikodonu_atstumu_matrica[2][0] = "bacterial3"
dikodonu_atstumu_matrica[3][0] = "bacterial4"
dikodonu_atstumu_matrica[4][0] = "mamalian1"
dikodonu_atstumu_matrica[5][0] = "mamalian2"
dikodonu_atstumu_matrica[6][0] = "mamalian3"
dikodonu_atstumu_matrica[7][0] = "mamalian4"


#bacterial1
dikodonu_atstumu_matrica[0][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[0][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[0][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[0][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[0][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[0][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=mamalian3_dikodonu_daznis)
dikodonu_atstumu_matrica[0][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial1_dikodonu_daznis, seka2_dazniai=mamalian4_dikodonu_daznis)

#bacterial2
dikodonu_atstumu_matrica[1][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[1][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[1][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[1][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[1][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[1][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=mamalian3_dikodonu_daznis)
dikodonu_atstumu_matrica[1][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial2_dikodonu_daznis, seka2_dazniai=mamalian4_dikodonu_daznis)

#bacterial3
dikodonu_atstumu_matrica[2][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[2][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[2][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[2][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[2][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[2][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=mamalian3_dikodonu_daznis)
dikodonu_atstumu_matrica[2][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial3_dikodonu_daznis, seka2_dazniai=mamalian4_dikodonu_daznis)

#bacterial4
dikodonu_atstumu_matrica[3][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[3][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[3][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[3][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[3][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[3][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=mamalian3_dikodonu_daznis)
dikodonu_atstumu_matrica[3][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=bacterial4_dikodonu_daznis, seka2_dazniai=mamalian4_dikodonu_daznis)

#mamalian1
dikodonu_atstumu_matrica[4][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[4][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[4][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[4][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[4][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[4][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = mamalian3_dikodonu_daznis)
dikodonu_atstumu_matrica[4][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian1_dikodonu_daznis, seka2_dazniai = mamalian4_dikodonu_daznis)

#mamalian2
dikodonu_atstumu_matrica[5][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[5][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[5][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[5][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[5][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[5][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = mamalian3_dikodonu_daznis)
dikodonu_atstumu_matrica[5][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian2_dikodonu_daznis, seka2_dazniai = mamalian4_dikodonu_daznis)

#mamalian3
dikodonu_atstumu_matrica[6][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[6][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[6][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[6][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[6][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[6][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[6][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian3_dikodonu_daznis, seka2_dazniai = mamalian4_dikodonu_daznis)

#mamalian4
dikodonu_atstumu_matrica[7][1] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = bacterial1_dikodonu_daznis)
dikodonu_atstumu_matrica[7][2] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = bacterial2_dikodonu_daznis)
dikodonu_atstumu_matrica[7][3] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = bacterial3_dikodonu_daznis)
dikodonu_atstumu_matrica[7][4] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = bacterial4_dikodonu_daznis)
dikodonu_atstumu_matrica[7][5] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = mamalian1_dikodonu_daznis)
dikodonu_atstumu_matrica[7][6] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = mamalian2_dikodonu_daznis)
dikodonu_atstumu_matrica[7][7] = gauti_dvieju_seku_atstuma(seka1_dazniai=mamalian4_dikodonu_daznis, seka2_dazniai = mamalian3_dikodonu_daznis)


pprint.pprint(dikodonu_atstumu_matrica)

#Paimti visų sekų visų kodonų dažnius, kiekvieną, kur sutampa kodonas, apskaičiuoti dažnių skirtumą ir įdėti į masyvą.
#Dažnių skirtumo suma / masyvo skaičiaus - šių dviejų sekų atstumas?





    


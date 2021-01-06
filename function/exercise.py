#!/usr/.binx/env python

import gzip
import copy
import parse_go_obo


user = input("user:")

dataset = ""

if user == "davide":
    dataset = "Q04206,Q06413,Q9NRY6,Q5VWK5,Q9BQ51,P10620,P01374,Q8N7H5,Q5D1E8,P81172,Q99735,Q15653,Q9Y6Y9,P30203,P31151,A6NGU5,P49675,P48023,P04179,P46663,Q96RT1,P23443,P16581,P05362,Q8N423,Q86V24,Q03181,P22301,Q02556,P35408,P32927,Q00613,P05109,P27338,Q8WUQ7,Q9Y2R2,Q9BQE4,P21731,P36268,Q9UBI6,P19875,Q9Y6Q6,P13716,P01100,P01588,P05089,P20333,P29460,P34810,Q9UBH0,P46531,P20366,Q92847,Q9NZH7,Q9Y616,P08311,P01210,P16615,Q9Y2G9,P45983,Q07325,Q9ULZ2,Q9Y2X7,P33681,Q92887,P25116,Q9UQL6,P17535,P20393,P16671,O15111,P02775,P09211,P58753,Q9UBK2,Q99727,Q9HCY8,Q9NZC2,Q9NZH8,Q9UJ14,P10646,P42830,Q8NFZ5,O43823,Q16822,P15428,P02778,P08571,Q6PD62,Q9NZZ3,Q9UHJ6,Q9ULZ3,Q9Y263,P17302,P35225,Q15842,P17861,P01584,P05412,Q86XR7,P08246,Q9BWU1,Q8TDQ0,Q9P2L0,O43927,P04798,P19876,A0A087X1L8,Q9UHA7,P19320,P05186,Q8WV44,P51617,P29033,O60674,P07101,Q5XLA6,P00519,P30530,Q96RI1,Q9UBC1,P41221,Q9Y6W6,Q96G74,P30793,A6NK06,P16885,P22309,P43119,Q13133,Q53H80,Q96FA3,P19838,O00206,P11684,P00966,P26583,Q96P20,P37840,Q9GZT4,Q9NRQ2,A1A4Y4,P05305,P17676,P59666,P34995,Q8IUC6,Q76LX8,Q96RU8,Q9NZH6,P17181,P09486,P29459,Q9BRX9,P35318,P42081,P21554,Q9H334,P40200,Q6P2Q9,P05177,O95452,P10720,P11802,P08670,Q13443,Q8NHL6,Q02952,P00797,Q15327,O14684,O95832,P06702,P78536,P11766,Q96KP6,Q5EG05,Q99467,P13500,Q9C030,P30536,P49238,P43088,O95751,P09543,Q13651,Q5VWQ8,P30048,O14625,O43541,P59665,P23510,Q9UK39,P42574,P51164,Q99684,Q9Y4K3,P24530,P01583,Q00653,Q92769,Q99836,P32248,Q8WX94,P21730,P12931,P21580,Q53EL6,P35228,O15229,Q07343,Q99966,O00182,P07204,Q9NZQ7,P10145,P04141,Q68CJ6,P02776,P05121,P16109,P32119,P29466,Q08397,Q5EBM0,P49137,Q01524,A6NDB9,Q15744,P07585,P16591,P34972,Q16539,P01148,P35372,P78324,Q16644,P61586,P29474,P51681,Q13151,O15528,P07911,Q14790,P09429,P21802,P09919,P25942,P26651,P35558,Q8WWZ1,Q9NYJ8,Q9BZS1,P49279,Q01523,O43766,Q13007,Q92900,P80162,P19440,P51812,O95999,P22897,P04114,Q9HC29,P18428,P17655,P12838,P05164,O15520,P29475,P31327,Q6P1J9,Q86U42,P03973,Q86VZ5,Q99732,P09341,P36269,P05231,P14902,Q01740,Q8WTV0"
elif user == "darko":
    dataset = "P49810,Q9H0U3,Q7Z3S7,Q92504,Q96RD7,Q13698,B1AL88,Q15413,Q8TDI7,Q96HE7,Q9BSW2,Q8WXH2,Q0D2K0,Q9UPR5,P26678,Q9H8M5,O75439,Q9UQ90,Q96H72,P32418,Q6PML9,O95180,Q8TDI8,Q6ZMH5,Q8TEL6,Q9NVV0,P98194,Q9Y6D0,P23634,Q13061,Q13563,Q16720,Q6P499,Q99726,Q9P0L9,Q12879,P48023,Q8N4V1,P14416,P21817,P06241,Q7RTP0,Q9NY26,Q96P56,Q9Y5I7,O60895,P80075,P29279,Q8NE86,P20020,O15554,Q6P4Q7,Q9UQM7,Q8NFF2,Q00975,Q9HCF6,Q4KMQ2,P28223,Q7Z2W7,Q14643,P62955,O15399,Q9H1D0,Q7RTX7,Q7RTS1,Q9BRQ5,Q9UJZ1,Q8IVJ1,Q01484,O00305,P11836,Q2M3R5,O60840,P13693,P54284,Q8N1S5,Q9HC07,P16615,Q7Z442,Q9H251,Q9UM00,Q08289,Q14571,Q6NXT4,P25116,Q06432,Q8IZS8,Q9BPX6,Q9BR39,Q9P0X4,Q9P246,Q9ULF5,P35638,O75844,P61619,P35232,P43681,Q9H4I9,P31146,P30626,P62942,Q13936,O95202,Q8IWU4,Q9NP94,P10147,O43497,P36544,Q00535,Q08209,O14983,O00555,Q8TC26,Q9HBA0,P78348,Q9H841,Q9H6F2,P32246,Q8TAD4,Q9HB89,Q8WWC4,P21796,Q13454,Q9Y698,Q86XE3,Q71RS6,Q10713,Q16602,Q99572,Q9HD23,Q02641,Q14957,P48995,Q9BV35,Q01668,P05814,P49281,Q9Y5S1,Q01814,P98161,O14863,P16885,Q9UF02,Q8NER1,Q8TCU5,P05771,P49768,P41180,Q9HDC5,Q5W0U4,Q9HC58,Q9Y210,Q9UGM1,Q8NEW0,Q9H300,Q8TD43,P13501,P68106,Q8N8Q9,Q92736,Q96JJ6,Q15043,O60391,Q9BRI3,Q9NY47,Q9UI40,Q9Y4W6,Q9BRY0,Q96QT4,Q9GZZ6,Q96TA2,P06756,P56704,Q93034,Q9P2D0,Q8IZK6,P17787,Q6XR72,P06239,O76027,Q15878,Q9UBN4,P35348,Q5VU97,Q93084,Q8IZF0,P46094,P62166,Q9NUM3,Q9GZQ4,Q9NTG1,Q14393,P63252,P51674,Q9C0K1,Q8N7U6,O60936,P54289,Q03135,Q13433,Q9NZQ8,Q99731,Q9NWR8,O60896,Q9UL62,Q9NZM6,P08575,Q13224,A8K7I4,Q96JW4,O75762,Q8NEC5,Q8NHX9,Q9BX84,Q9HCX4,P07355,Q6P5W5,Q8TDX9,Q96BY9,Q9NQA5,Q13586,P32248,Q8WXS4,O00429,P57103,Q14573,Q9ULQ1,O60721,Q8TDD5,Q8WXS5,P25713,Q96D31,Q96AQ8,Q8IU99,Q13507,Q8NET8,P47992,Q9H244,P35372,O75185,Q504Y0,Q8IUH5,P51681,Q99571,O15528,Q8IUH4,Q7Z4N2,Q05586,Q6J4K2,Q9GZU1,O60894,P08133,Q9BXT2,Q86XQ3,P28335,P49279,Q6NVV3,O00631,O60313,Q6KCM7,Q8IWX8,P08195,Q96SN7,Q99623,P35212,Q9Y6N3,Q86YW0,O94759,Q9Y6M5,Q8IYU8,O00585,O60359,P09038,Q7Z443,O75949,Q14644,Q8N4Y2,Q9UBN1,Q9Y271,A6NMY6,P41595"
dataset = dataset.split(",")


def gen_block(f):
    """
    Genrator function that parses GOA GAF files (https://www.ebi.ac.uk/GOA/downloads)
    The generator yields a block of lines corresponding to the same protein
    UniProtKB       A0A024R1R8      hCG_2014768             GO:0002181      PMID:21873635   IBA     PANTHER:PTN002008372|SGD:S000007246     P       HCG2014768, isoform CRA_a       hCG_2014768     protein taxon:9606      20171102        GO_Central
    UniProtKB       A0A024RBG1      NUDT4B          GO:0003723      GO_REF:0000037  IEA     UniProtKB-KW:KW-0694    F       Diphosphoinositol polyphosphate phosphohydrolase NUDT4B NUDT4B  protein taxon:9606      20191109        UniProt
    UniProtKB       A0A024RBG1      NUDT4B          GO:0005829      GO_REF:0000052  IDA             C       Diphosphoinositol polyphosphate phosphohydrolase NUDT4B NUDT4B  protein taxon:9606      20161204        HPA
    """
    name, old_name = None, None
    chunk = []
    for line in f:
        line = line.decode()
        if line and line[0] != "!":
            _, name, _, _, term, _, ec, _, namespace, protein_name = line.split("\t")[
                :10]
            if name != old_name and old_name:
                # return a set as there can be repetitions, i.e. the same term with different evidence codes
                yield (old_name, set(chunk))
                chunk = []
            old_name = name
            chunk.append(term)
    # Last line
    if old_name:
        yield (old_name, set(chunk))


if __name__ == "__main__":

    # Get ontology data from parse_go_obo lib
    graph = parse_go_obo.parse_obo("../data/function/go.obo")
    ancestors, depth, roots = parse_go_obo.get_ancestors(graph)
    # children = parse_go_obo.get_children(ancestors)

    proteins = {}  # { term : proteins_annotated_with_term } It contains proteins annotated directly with the term or its children
    with gzip.open("../data/function/goa_human.gaf.gz") as f:
        for acc, annotations in gen_block(f):
            # Copy direct annotations
            terms = copy.copy(annotations)
            # Add ancestors
            for term in annotations:
                terms.update(ancestors.get(term, set()))
            # For each term add protein accession to proteins dict
            for term in terms:
                proteins.setdefault(term, set()).add(acc)
    print(graph["GO:0043549"]["def"], len(proteins["GO:0043549"]))
    print(graph["GO:0005739"]["def"], len(proteins["GO:0005739"]))

    terms_set = {}  # { term : count }  dataset proteins
    terms_rest = {}  # { term : count }  other proteins
    proteins_set = 0  # number of dataset proteins
    proteins_rest = 0  # number of remaining proteins
    with gzip.open("../data/function/goa_human.gaf.gz") as f:
        for acc, annotations in gen_block(f):

            # Copy direct annotations
            terms = copy.copy(annotations)

            # Add ancestors
            for term in annotations:
                terms.update(ancestors.get(term, set()))

            # For each term add protein accession to proteins dict
            if acc in dataset:
                proteins_set += 1
                for term in terms:
                    terms_set.setdefault(term, 0)
                    terms_set[term] += 1
            else:
                proteins_rest += 1
                for term in terms:
                    terms_rest.setdefault(term, 0)
                    terms_rest[term] += 1

    data = []
    for term in terms_set:
        ratio_set = (terms_set[term] + 1) / proteins_set  # add pseudo count
        ratio_rest = terms_rest.get(
            term, 1) / proteins_rest  # add pseudo count
        fold_increase = ratio_set / ratio_rest
        data.append((term, terms_set[term], terms_rest.get(
            term, 0), ratio_set, ratio_rest, fold_increase))
    for ele in sorted(data, key=lambda x: x[5], reverse=True)[:20]:
        print("{} {} {} {:.2g} {:.2g} {:.2g} {} {}".format(
            *ele, depth[ele[0]], graph[ele[0]]))

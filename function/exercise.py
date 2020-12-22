#!/usr/binx/env python

import gzip
import copy
import parse_go_obo


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
            _, name, _, _, term, _, ec, _, namespace, protein_name = line.split("\t")[:10]
            if name != old_name and old_name:
                yield (old_name, set(chunk))  # return a set as there can be repetitions, i.e. the same term with different evidence codes
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

    dataset = "P48742,Q06413,Q8NFJ8,Q8IZE3,P21453,P23769,P08476,P57073,Q9HCJ5,Q9NY30,Q9NYK6,O60443,Q9Y6F9,O43293,P35222,Q04741,Q9BTL4,Q8N474,Q9NP95,Q9UBV4,O60469,Q969G9,Q02535,O43915,Q9UG01,Q9H2A3,P23435,Q9UIW0,Q6IQ32,Q03395,O95096,Q05925,Q96T92,Q7RTU3,Q9Y6N9,Q9Y6X8,P04628,Q8N0W4,Q13402,O14905,Q9BZM3,Q9GZT8,O00744,Q2MV58,A8MTQ0,Q9Y6K1,Q12837,Q13516,Q7RTS1,Q9H0D6,P55145,P32243,Q86U70,Q86Z02,Q16619,Q9NZN1,O75084,Q8NDY6,Q16620,Q14332,P78414,Q01851,P56706,Q9UP38,Q9Y231,Q9UL68,O00755,Q9UQL6,O43424,Q04656,Q9Y466,P15692,Q13950,P78415,O94900,P11021,P61968,Q92886,Q969G2,Q01101,Q9H1J7,P20265,P52954,P23468,O43490,O75445,P00492,P18074,P55317,Q76NI1,Q99453,Q9ULW2,Q92858,P04629,Q00535,Q9UBR4,Q9NSC2,P12644,P41134,Q8TE12,Q9C0K0,Q5SQQ9,Q9H1J5,P78412,Q13761,Q03468,Q92753,P56178,Q15303,Q92786,P19320,Q96SQ7,P35716,P17542,Q9BWQ8,Q93097,P04626,Q7RTS3,Q9P2K5,P41221,P61201,Q01814,Q8NCM8,Q99929,O95573,Q6PIY7,Q8N111,Q96NK8,P49768,Q9UMR3,P43354,Q9H2C1,Q08629,Q9Y210,Q96NZ1,Q9UBY8,P09544,P57796,Q6UXX9,Q9Y261,P17676,P54920,Q9GZZ0,P43699,P35398,Q15319,Q9H461,Q9NQX7,Q8TAG6,Q16665,P78413,Q9ULX5,Q9Y4Z2,P56704,O75386,P78504,O95922,P36873,Q07687,P56703,Q96KG9,A6NJT0,P18848,Q9BZI1,Q9H2X6,P40763,Q99697,Q9UPW5,Q96SC8,P14138,O14813,Q9NQ69,P55075,O75581,Q8N158,P61371,Q9UKV0,P08581,Q13467,O60890,Q15751,P78426,P48436,P50458,P34130,Q9NPE2,P42574,Q9UL62,Q8NBI3,Q93098,P04062,P19447,Q7Z553,Q13547,P56705,Q04743,Q9NQ38,Q9NY43,Q9GZT5,P58546,Q9H4S2,O96014,Q13562,Q49AH0,Q9NX09,P21246,P23352,Q9NQZ6,Q96Q05,Q8TAK6,Q7Z4P5,O60663,Q8TDD5,P34925,Q96A47,Q9UPM6,Q8N100,Q68G74,Q58EX2,Q02750,Q15784,Q8NFP4,Q9NZR4,Q8TDC3,Q15078,Q5XKR4,Q8NBP7,P56177,Q14549,Q9ULV1,Q8TDI0,O94788,P19622,P78411,Q13526,O15079,O75364,P07333,Q8IWQ3,Q06945,P61296,O60488,Q15465,Q96ST3,Q01196,Q9BYB0,Q9H4Z2,P50553,Q9NPG1,P28370,Q9HD90,Q9UM54,Q99835,O14904,Q9NSY0,Q13635"
    dataset = dataset.split(",")

    terms_set = {}  # { term : count }  dataset proteins
    terms_rest = {}  #  { term : count }  other proteins
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
        ratio_rest = terms_rest.get(term, 1) / proteins_rest  # add pseudo count
        fold_increase = ratio_set / ratio_rest
        data.append((term, terms_set[term], terms_rest.get(term, 0), ratio_set, ratio_rest, fold_increase))
    for ele in sorted(data, key=lambda x: x[5], reverse=True)[:20]:
        print("{} {} {} {:.2g} {:.2g} {:.2g} {} {}".format(*ele, depth[ele[0]], graph[ele[0]]))


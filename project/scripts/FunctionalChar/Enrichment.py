import gzip
import copy
''' start functions'''
''' parse_go_obo '''
def parse_obo(obo_file):
    # Parse the ontology, exclude obsolete terms
    graph = {}  # { term_id : term_object }
    obj = {}  # { id: term_id, name: definition, is_a: list_of_parents, is_obsolete: True, namespace: namespace }
    with open(obo_file) as f:
        for line in f:
            line = line.strip().split(": ")
            if line and len(line) > 1:
                # print(line)
                k, v = line[:2]
                if k == "id" and v.startswith("GO:"):
                    obj["id"] = v
                elif k == "name":
                    obj["def"] = v
                elif k == "is_a" and v.startswith("GO:"):
                    obj.setdefault("is_a", []).append(v.split()[0])
                elif k == "is_obsolete":
                    obj["is_obsolete"] = True
                elif k == "namespace":
                    obj["namespace"] = v
            else:
                if obj.get("id") and not obj.get("is_obsolete"):
                    if "is_a" not in obj:
                        obj["is_root"] = True
                    graph[obj["id"]] = obj
                    # print(obj)
                obj = {}
    return graph
def get_ancestors(graph):
    """
    Build a dictionary of ancestors
    and calculate term depth (shortest path)
    """

    roots = set()
    for node in graph:
        if graph[node].get("is_root"):
            roots.add(node)

    depth = {}
    ancestors = {}  # { term : list_of_ancestor_terms }
    for node in graph:
        c = 0
        node_ancestors = []
        node_parents = graph[node].get("is_a")

        # Loop parents levels (parents of parents) until no more parents
        while node_parents:
            c += 1

            # Set root
            if node not in depth and roots.intersection(set(node_parents)):
                depth[node] = c

            # Add ancestors
            node_ancestors.extend(node_parents)

            # Update the list of parents (1 level up)
            node_parents = [term for parent in node_parents for term in graph[parent].get("is_a", [])]

        ancestors[node] = set(node_ancestors)
    return ancestors, depth, roots
def get_children(ancestors):
    children = {}  # { node : list_of_children }, leaf terms are not keys
    for node in ancestors:
        for ancestor in ancestors[node]:
            children.setdefault(ancestor, set()).add(node)
    return children
''' function that parses the GOA.gaf files '''
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
''' end functions'''
graph = parse_obo("scripts/FunctionalChar/go.obo")
ancestors, depth, roots = get_ancestors(graph)
children = get_children(ancestors)
with open('data/Family Sequences/ids.txt') as f:
    for line in f:
        rec_id = line.strip().split()
# rec_id -> sequence of proteins id (acc)

terms_set_raw = {}  # { term : count }  no ancestors
terms_rest_raw = {}  # { term : count }  no ancestors
terms_set = {}  # { term : count }
terms_rest = {}  # { term : count }  other proteins
proteins_set = 0  # number of mitochondrial proteins
proteins_rest = 0  # number of remaining proteins

with gzip.open("data/goa_human.gaf.gz") as f:
    for acc, annotations in gen_block(f):
        # example: acc -> W6CW81
        # annotations -> {'GO:0003690', 'GO:0002218', 'GO:0098586', 'GO:0035458', 'GO:0031333', 'GO:0097169', 'GO:0006954', 'GO:0005515'}
        annotations = set(filter(lambda x: graph.get(x), annotations))
        if annotations:
            terms = copy.copy(annotations) 
            for term in annotations:
                terms.update(ancestors.get(term, set()))
                
                if acc in rec_id:
                    proteins_set += 1
                    for term in terms:                     # terms -> all annotations with ancestors 
                        terms_set.setdefault(term, 0)
                        terms_set[term] += 1
                    for term in annotations:               # annotations (w/o ancestors)
                        terms_set_raw.setdefault(term, 0)
                        terms_set_raw[term] += 1
                else:
                    proteins_rest += 1
                    for term in terms:
                        terms_rest.setdefault(term, 0)
                        terms_rest[term] += 1
                    for term in annotations:
                        terms_rest_raw.setdefault(term, 0)
                        terms_rest_raw[term] += 1

print(sum([terms_set_raw[term] for term in terms_set_raw]))
print(sum([terms_rest_raw[term] for term in terms_rest_raw]))
print(len(terms_set_raw))
print(len(terms_rest_raw))
print(len([term for term in terms_set_raw if children.get(term) is None]))
print()
from scipy.stats import fisher_exact
data = []
for term in terms_set:
    ratio_set = (terms_set[term] + 1) / (proteins_set - terms_set[term] + 1)
    ratio_rest = terms_rest.get(term, 1) / (proteins_rest - terms_rest.get(term, 0) + 1)
    fold_increase = ratio_set / ratio_rest
    p_value = fisher_exact([[terms_set[term] + 1, proteins_set - terms_set[term] + 1],
                           [terms_rest.get(term, 1), proteins_rest - terms_rest.get(term, 0) + 1]], 'two-sided')
    data.append((term, terms_set[term] + 1, terms_rest.get(term, 1), ratio_set, ratio_rest, fold_increase,
                 p_value[1], graph[term]["def"], graph[term]["namespace"]))

# Retrieving the go terms that are enriched
enriched_terms = {}
e_t = []
for line in list(filter(lambda k: k[5] > 1.0, data)):
    enriched_terms.setdefault(line[0], 0)
    enriched_terms[line[0]] += line[1]
enriched_terms
from wordcloud import WordCloud
import matplotlib.pyplot as plt
wc = WordCloud(background_color="white",width=1000,height=1000, max_words=1000,relative_scaling=0.5,normalize_plurals=False).generate_from_frequencies(enriched_terms)
plt.imshow(wc)
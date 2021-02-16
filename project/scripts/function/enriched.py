
# ? Calculate the enrichment of each term in the dataset compared to GO annotations available in the SwissProt database
# ? (you can download the entire SwissProt XML here). You can use Fisher’ exact test and verify that both two-tails
# ? and right-tail P-values (or left-tail depending on how you build the confusion matrix) are close to zero.
from go_parser import parse_obo, get_ancestors, get_children, gen_block
from plot import create_word_cloud
from scipy.stats import fisher_exact
import gzip
import copy
import pandas as pd

graph = parse_obo("../../data/function/go.obo")
ancestors, depth, roots = get_ancestors(graph)
children = get_children(ancestors)

df = pd.read_csv("../../data/datasets/family sequences/family_sequences.csv")
prot_list = list(df['GO ID'])

terms_set_raw = {}  # { term : count }  no ancestors
terms_rest_raw = {}  # { term : count }  no ancestors
terms_set = {}  # { term : count }
terms_rest = {}  # { term : count }  other proteins
proteins_set = 0  # number of mitochondrial proteins
proteins_rest = 0  # number of remaining proteins

with gzip.open("../../data/function/goa_human.gaf.gz") as f:
    for acc, annotations in gen_block(f):
        annotations = set(filter(lambda x: graph.get(x), annotations))
        if annotations:
            terms = copy.copy(annotations)
            for term in annotations:
                terms.update(ancestors.get(term, set()))

                if acc in prot_list:
                    proteins_set += 1
                    for term in terms:                     # terms -> all annotations with ancestors
                        terms_set.setdefault(term, 0)
                        terms_set[term] += 1
                    # annotations (w/o ancestors)
                    for term in annotations:
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

data = []
for term in terms_set:
    ratio_set = (terms_set[term] + 1) / (proteins_set - terms_set[term] + 1)
    ratio_rest = terms_rest.get(
        term, 1) / (proteins_rest - terms_rest.get(term, 0) + 1)
    fold_increase = ratio_set / ratio_rest
    p_value = fisher_exact([[terms_set[term] + 1, proteins_set - terms_set[term] + 1],
                            [terms_rest.get(term, 1), proteins_rest - terms_rest.get(term, 0) + 1]], 'two-sided')
    data.append((term, terms_set[term] + 1, terms_rest.get(term, 1), ratio_set, ratio_rest, fold_increase,
                 p_value[1], graph[term]["def"], graph[term]["namespace"]))

# ? enriched terms

f = open('../../data/function/enriched_terms.txt', 'w+')
columns = ['Term', '#Term Set', '#Term Rest', 'Ratio Set', 'Ratio Rest',
           'Fold Increase', 'p_value', 'Definition', 'Sub-Ontology']
for namespace in ['molecular_function', 'biological_process', 'cellular_component']:
    f.write('Sub-Ontology: {}\n'.format(namespace))
    f.write('{:<10} | {:<10} | {:<10} | {:<10} | {:<15} | {:<20} | {:<10} | {:<70} | {:<30}\n'.format(*columns))
    data_sub = list(filter(lambda x: x[8] == namespace, data))
    for ele in sorted(data_sub, key=lambda x: x[5], reverse=True):
        f.write(
            '{:<10} | {:<10} | {:<10} | {:<10.2} | {:<15.2} | {:<20.2} | {:<10.2} | {:<70} | {:<30}\n'.format(*ele))
    f.write('----'*50+'\n')

print(f"Enriched terms: {len(data)}")
print("enriched_terms.txt created")

enriched_terms = {}

for term in list(filter(lambda k: k[5] > 1.0, sorted(data, key=lambda x: x[5], reverse=True))):
    enriched_terms.setdefault(term[0], 0)
    enriched_terms[term[0]] += term[1]

# ? enriched branches

GO_Terms = []
with open('../../data/function/go_annotations_count.txt') as f:
    for line in f:
        GO_Terms.append(line.strip().split()[0])

GO_Terms_Parent = []
for term in GO_Terms:
    if term in children.keys():
        GO_Terms_Parent.append(term)

f = open('../../data/function/enriched_branches.txt', 'w+')
for namespace in ['molecular_function', 'biological_process', 'cellular_component']:
    f.write('Sub-Ontology: {}\n'.format(namespace))
    f.write('{:<10} | {:<10} | {:<10} | {:<10} | {:<15} | {:<20} | {:<10} | {:<70} | {:<30}\n'.format(*columns))
    data_sub = list(filter(lambda x: x[8] == namespace, data))
    data_sub2 = list(filter(lambda x: x[0] in GO_Terms_Parent, data_sub))
    for ele in sorted(data_sub2, key=lambda x: x[5], reverse=True):
        f.write(
            '{:<10} | {:<10} | {:<10} | {:<10.2} | {:<15.2} | {:<20.2} | {:<10.2} | {:<70} | {:<30}\n'.format(*ele))
    f.write('----' * 50 + '\n')
f.close()
print("enriched_branches.txt created")

# ? World Cloud
create_word_cloud(enriched_terms)

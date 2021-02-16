from ete3 import NCBITaxa, TreeStyle, AttrFace

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

taxonomy = {}
''' Import Taxonomic IDs '''
with open('../../data/taxonomy/taxonomic_lineage.txt') as f:
    next(f)
    for line in f:
        taxonomy[line.strip().split(' | ')[0]] = line.strip().split(' | ')[1]
tax_id = []
for rec in taxonomy:
    tax_id.append(taxonomy[rec])

tree = ncbi.get_topology(tax_id)


def my_layout(node):
    if getattr(node, "rank", None):
        rank_face = AttrFace("rank", fsize=7, fgcolor="indianred")
        node.add_face(rank_face, column=0, position="branch-top")
    if node.is_leaf():
        sciname_face = AttrFace("sci_name", fsize=9, fgcolor="steelblue")
        taxid_face = AttrFace('taxid', fsize=9, fgcolor="limegreen")
        node.add_face(sciname_face, column=0, position="branch-right")
        node.add_face(taxid_face, column=1, position="branch-right")
    else:
        sciname_face = AttrFace("sci_name", fsize=7, fgcolor="steelblue")
        taxid_face = AttrFace('taxid', fsize=7, fgcolor="limegreen")
        node.add_face(sciname_face, column=1, position="branch-bottom")
        node.add_face(taxid_face, column=1, position="branch-bottom")


ts = TreeStyle()
ts.layout_fn = my_layout
ts.show_leaf_name = False
tree.render("../../data/taxonomy/tree.jpg", w=1980, tree_style=ts)

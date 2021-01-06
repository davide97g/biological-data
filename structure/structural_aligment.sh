
# TM-align
# https://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.gz
# WARNING the cpp version gives different output files (use the Fortran version whenever possible)

# Pymol
# https://pymolwiki.org/index.php/Windows_Install
# sudo apt install pymol
#



# Align two structures
../.binx/TMalign ../data/structure/davide/pdb1iuj_cut.ent ../data/structure/davide/pdb1lq9_cut.ent -o ../data/structure/davide/1iuj_1lq9.sup > ../data/structure/davide/1iuj_1lq9.out
# Visualize the alignment with Pymol (run in background, redirect stderr and stout to /dev/null)
pymol ../data/structure/1ucd_1ioo.sup_all_atm.pml > /dev/null 2>&1 &

# Align with TMalign two different structures
# ../.binx/TMalign ../data/structure/101m.pdb ../data/structure/1mba.pdb -o ../data/structure/101m_1mba.sup > ../data/structure/101m_1mba.out
# pymol ../data/structure/101m_1mba.sup_all_atm.pml > /dev/null 2>&1 &

# Compare sequence identity with the TMalign output



# Align the 2 structures starting from a sequence alignment. Generate the sequence alignment with
# Needle (https://www.ebi.ac.uk/Tools/psa/emboss_needle/)
# ./.binx/TMalign ../data/structure/davide/pdb1iuj.ent ../data/structure/davide/pdb1lq9.ent -I ../data/structure/1ucd_1ioo_needle.fa -o ../data/structure/1ucd_1ioo_needle.sup > ../data/structure/1ucd_1ioo_needle.out
# pymol ../data/structure/1ucd_1ioo_needle.sup_all_atm.pml  > /dev/null 2>&1 &


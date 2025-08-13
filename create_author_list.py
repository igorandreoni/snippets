# Author: Michael Coughlin, Igor Andreoni
import os
import sys

# EDIT! Create a dictionary where each author name is the key to which is assigned a list of institutions
authors = {'Igor Andreoni': ['Joint Space-Science Institute, University of Maryland, College Park, MD 20742, USA', 'Department of Astronomy, University of Maryland, College Park, MD 20742, USA', 'Astrophysics Science Division, NASA Goddard Space Flight Center, Mail Code 661, Greenbelt, MD 20771, USA'],
           'Michael W. Coughlin': ['School of Physics and Astronomy, University of Minnesota, Minneapolis, Minnesota 55455, USA'],
           'S. Bradley Cenko': ['Joint Space-Science Institute, University of Maryland, College Park, MD 20742, USA', 'Astrophysics Science Division, NASA Goddard Space Flight Center, Mail Code 661, Greenbelt, MD 20771, USA'],
    }

# Create an ordered list of institutions
institution_list_ordered = []
for key in authors.keys():
    instutions = authors[key]
    for instution in instutions:
        if instution in institution_list_ordered: continue
        institution_list_ordered.append(instution)

# Create the author list
author_list = []
author_list_arxiv = []
for key in authors.keys():
    author_institutions = authors[key]
    indices = []
    for author_institution in author_institutions:
        indices.append(institution_list_ordered.index(author_institution))
    author_list.append('%s$^{%s}$' % (key, ",".join([str(x+1) for x in indices])))
    author_list_arxiv.append(key)

# Print the author list for ArXiv
print("List of authors for ArXiv submission:")
print("")
print(", ".join(author_list_arxiv))
print("--------\n")

# Print the author list for the paper
print("List of authors for the paper:")
print("")
print(", ".join(author_list))
print("")

# Print the list of institutions
institution_list = []
for ii, institution in enumerate(institution_list_ordered):
    institution_list.append('$^{%s}$ %s' % (str(ii+1), institution))
print("\n".join(institution_list)) 
print("--------\n")

# Print the total number of authors and institutions
print(f"Total of {len(authors)} authors and {len(institution_list_ordered)} institutions")

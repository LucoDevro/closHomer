import json
import pandas as pd
import sys

path = sys.argv[1]
output = sys.argv[2]

with open(path,"r") as handle:
    cont = json.load(handle)

kegg_annots=[] # will contain the hierarchical annotation of each gene family
for a in cont["children"]:
    an = ' '.join(a["name"].split(' ')[1:]).split(' [')[0]
    try:
        for b in a["children"]:
            bn = ' '.join(b["name"].split(' ')[1:]).split(' [')[0] # parse name
            try:
                for c in b["children"]:
                    cn = ' '.join(c["name"].split(' ')[1:]).split(' [')[0]
                    try:
                        for d in c["children"]:
                            dn, term = d["name"].split('  ')[:2]
                            record = {"A": an, 
                                      "B": bn, 
                                      "C": cn, 
                                      "D": dn,
                                      "term": term}
                            kegg_annots.append(record) # A full 4-level KEGG annotation record for this gene family
                    except KeyError: # Stop drilling down if there is no deeper annotation level (i.e. when there is no 'children' key)
                        continue
            except KeyError:
                continue
    except KeyError:
        continue

# Convert list of dictionary-like records to a dataframe
kegg_annots = pd.DataFrame(kegg_annots)
kegg_annots.to_csv(output, index = False, header = True, sep = "\t")

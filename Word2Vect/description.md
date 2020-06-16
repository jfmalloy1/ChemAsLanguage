# Word2Vec Description
Describes contents of Word2Vec directory

## Purpose
How are the word embeddings of chemistry - both functional groups and chemical fragments - different than those of human language?
These embeddings will also be used to quantify and track changes over time within chemical evolution.

## Files (and descriptions)
### build_cooccurance_matrix.py
Build co-occurance matrix - needed for word2vec embeddings - from a list of SMARTS chemicals. The matrix is outputted to a .csv file.
Each chemical pairing is compared to all molecules within KEGG, and if the molecule pair appears together, this is a co-occurance.

### common_fragments.py
Adopted from Cadeddeu 2014. Obtains a SMARTS representation list of chemical fragments which match Zipf's law.

### frags.txt
List of fragments matching Zipf's Law. In SMARTS representation.

### frags_coocurrance_V1.csv
Un-transposed matrix of fragments

### frags_coocurrance_full.csv
Full matrix of fragments - both triangles (lower and upper) match

### functional_groups_coocurrance.csv
Full matrix of chemical functional groups - taken from Cadaddeu 2014, supplemental

### functional_groups_smarts.txt
List of functional groups in SMARTS representation

### get_kegg_smiles.py
Returns KEGG molecules in SMILES representation

### kegg_data.csv
CSV of KEGG database - contains chemical name, SMILES representation. molecular weight, etc...

### kegg_smiles.txt
List of KEGG molecules in SMILES representation

### transpose_matrix.py
Cleans up frags_coocurrance_V1.csv in order to make it a full matrix. Output is frags_coocurrance_full.csv

# Word2Vec Description
Describes contents of Word2Vec directory

## Purpose
How are the word embeddings of chemistry - both functional groups and chemical fragments - different than those of human language?
These embeddings will also be used to quantify and track changes over time within chemical evolution.

## Files (and descriptions)
### AA_random_V0.png
Image of amino acids (KEGG-defined class) vs random compounds

### br08001.json
Compound classes defined by KEGG (ex: Organic Acids, Lipids, etc...)

### build_KEGG_gensim.py
Builds a Word2Vec model - based on Cadaddeu fragments - using all of KEGG. Saves vectors to vectors_fullKEGG.kv

### build_cooccurance_matrix.py
Build co-occurance matrix - needed for word2vec embeddings - from a list of SMARTS chemicals. The matrix is outputted to a .csv file.
Each chemical pairing is compared to all molecules within KEGG, and if the molecule pair appears together, this is a co-occurance.

### cpd_similarity.py
Embeds KEGG compounds - defined by classes within KEGG - using pre-defined Word2Vec models and TSNE.

### common_fragments.py
Adopted from Cadeddeu 2014. Obtains a SMARTS representation list of chemical fragments which match Zipf's law.

### cpd_classes_V0.png
Initial image of cpd classes within KEGG

### frags.txt
List of fragments matching Zipf's Law. In SMARTS representation.

### frags_coocurrance_V1.csv
Un-transposed matrix of fragments

### frags_coocurrance_full.csv
Full matrix of fragments - both triangles (lower and upper) match

### frags_linesentence.txt
Fragments (shuffled) within each compound of KEGG

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

### load_KEGG_gensim.py
Test to load pre-built vectors (built in build_KEGG_gensim.py)

### transpose_matrix.py
Cleans up frags_coocurrance_V1.csv in order to make it a full matrix. Output is frags_coocurrance_full.csv

### vectors_01.kv
Initial test of saving vectors

### vectors_fullKEGG.kv
Fragment vectors, trained on all KEGG

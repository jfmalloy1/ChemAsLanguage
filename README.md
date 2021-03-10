# Chemistry as Language
========

ChemAsLang is a collection of resources which treats chemistry as language. It breaks sets of chemistry down into component parts (fragments), implements networks of fragments, and trains Word-to-Vector algorithms on fragments to track changes in chemistry over time.

Features
--------

- Utilizes RDKit generate chemical fragments
- Uses Gensim to train Word-to-Vector models of fragments & chemical compounds

Usage
------------
- To generate fragments: upload a set of chemicals (in SMILES or SMARTS format) to the `Fragments/Data` directory, edit `Fragments/common_frags_parallel.py` as necessary.
- To train W2V: edit `Word2Vec/build_KEGG_gensim.py` as necessary.

Progress
-------

This is an ongoing project as part of my graduate research at Arizona State Univeristy. Future plans are to:

- Simplify the workflow for generating fragments from a list of user-specified chemical compounds (estimated data: Summer 2021)
- Train a W2V model as part of the fragment generation process (estimated date: Fall 2021)
- Release this code as a fully documented Python package (estimated data: Spring 2022)

Support
-------

If you have questions, email me at: jfmalloy1@gmail.com

License
-------

The project is licensed under the BSD license.

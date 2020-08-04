import pandas as pd
import numpy as np
import re
import itertools
import random
import os
from gensim.models import Word2Vec
from gensim.test.utils import get_tmpfile
from gensim.models import KeyedVectors

def main():
    word2vec = KeyedVectors.load("vectors_fullKEGG.kv", mmap="r")
    v1 = word2vec.wv["[#6]-[#6]"]
    print(v1)
    
    sim_frags = word2vec.wv.most_similar("[#6]-[#6]")
    print(sim_frags)


if __name__ == "__main__":
    main()

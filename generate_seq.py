#!/usr/bin/env python
# coding: utf-8

import numpy as np
import random
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns


# # Generate a random DNA sequence
n = 1000
def generate_seq(n):
    return ''.join(random.choices('ACTG',k=n))
seq = generate_seq(n)





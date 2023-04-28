# strand_displacement_genome

In this project, we are attempting to detect potential instances of strand displacement in various genomes. In simple terms, strand displacement occurs when a sequence contains the following pattern: A A* A, such that A* is a sequence that is complementary to A. The paired region A A* is displaced by the second A, which is called the invading strand. The pattern can also be A* A A, or A A A*. The scripts in this repository attempt to detect all 3 types.

Currently, the script allows for mismatches between the invading and substrate strands. The next step would be to allow mismatches between the complementary and substrate strands as well.

Below is a summary of each file/directory:

### data/
Contains the E. coli genome (fasta file) and genome annotation (gff file). Used by the bacteria.py script.
### test_data/
Contains fasta files for each RNA family downloaded from the RFAM database. The ID of each family is in the format 'RF#####'. Additionally, each family is tested using a variety of mismatch thresholds. This number is the number of mismatches that are allowed between the invading strand and its complementary strand.
### bacteria.py
Runs the strand displacement script on the E. coli genome. It's called "bacteria.py" but will work on any genome. It uses a sliding window approach, where it scans in windows 150 nucleotides long across the entire genome, determining the longest sequence possible strand displacement sequences in each window. The output is a dataframe with two columns: the length of the longest possible invading strand at each window, and the entire window sequence.
### combine_gf.py
This is a short script used in the analysis notebook that merges the genome annotation and the output of the strand displacement script. This merging is necessary because the GFF file annotates ranges of sequences rather than individual nucleotides.
### find_compl.py
This is a simplified version of the strand displacement problem that I initially worked on. Given a sequence, it finds the longest sub-sequence A such that A* exists in the sequence.
### generate_seq.py
Generates a random sequence with a specified length. Useful for testing out the strand displacement scripts and comparing results on actual genomic sequences to random sequences.
### genome_results.ipynb
A jupyter notebook that analyzes results from the bacteria.py script. It 
### quick_fasta_search3.cpp
Script written by Prof. Petr Sulc to find complementary strands in a genome.
### results.ipynb
Analyzes output of sd_finder.py. Generates plots to compare the RFAM families to random sequences. Additionally, it does the comparison for riboswitch and non-riboswtich families
### sd_finder.py
Initial strand displacement script. Used for computing strand displacement in RNA families and random sequences, but it will work for any DNA sequences.

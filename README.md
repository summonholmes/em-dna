# em-algorithm-python
This is the Expectation-Maximization (EM) algorithm performed on multiple DNA sequences.  As a latent variable model, the EM algorithm is an example of unsupervised learning.  The primary objectives of this program are to discover optimal alignment positions and motifs for multiple DNA sequences; and to provide a clear, efficient EM implementation.  

Optimal alignment positions and motifs are the latent variables requiring detection; meaning that they are known to exist, are known to possess high importance, but are unable to be identified directly.

Since the starting positions of the DNA motifs are always randomized, no two runs are ever alike.  However, when provided the same input more than once, the program results are (almost) always the same.

The Expectation (E) step composes the majority of source code.  Total positional counts, frequencies, odds, then log-odds calculations compose the E step and derive the log-likelihood function.  The Maximization Step scores on the log-likelihood function derived from the previous Expectation Step.  Taking the score for every possible contiguous motif, for each sequence, the motif positions that generate the maximum scores on each sequence are identified, and become the input of the next E step.  The Expectation and Maximization Steps are then repeated X number of times.  After enough iterations, convergence occurs to identify the optimal alignment positions and motifs.

A brief description of the five classes:
1. EM_Input: Loads all user data including the number of random alignments, the number of iterations, and the FASTA file.  This class had an interactive mode which is commented out.  The defaults include a motif width of 6, 50 total random alignments, and 50 iterations.

2. EM_Core: The first iteration or 'for' loop over the number of random alignments.  Previously, the entire program could be thought of as one gigantic, quintuple, nested 'for' loop.  However, this approach has been reduced to a triple 'for' loop via numpy vectorization and indexing.  This is the 'nucleas' of the program, which encompasses the entire EM process and iterates for the specified number of total random alignments.

3. EM_Count: Performs all counting operations regarding bases, motifs, and normalization.  The only area where looping occurs is over the motif width when storing positional counts.  The rest of this class aims to be as Pythonic and vectorized as possible.

4. EM_Matrix: Consists of the E step, and generates the log odds matrix for the M step.  Numpy matrices are preferred here for their vector capabilities.  Child of EM_Count.

5. EM_Run: Consists of remaining inner two 'for' loops from mentioned above in EM_Core, and the M step.  Scores are calculated using a numpy character array for each motif and the log odds matrix.  The entire process repeats to convergence, and the final results are printed in 'main.py'.  Note that the function 'update_log_odds' is the entire Expectation Step and recycles all functions from EM_Count & EM_Matrix.  Child of EM_Matrix.

## Getting Started
The program requires few dependences and should be trivial to set up.  However, an in-depth understanding of the EM algorithm requires some knowledge of bioinformatics, data science, and machine learning.

### Dependenciese:
* python3  
* python3-biopython  
* python3-numpy
* python3-pprint  

### Usage:
#### Mac/Linux/Unix
```
$ python3 /path/to/main.py
```
#### Windows
```
C:\path\to\main.py C:\path\to\python3.exe main.py
```
### Example Run:
```
shanekimble@Shanes-MacBook-Pro  ~/em-algorithm-python   master  python main.py
Please specify the width of the motif: 6
Use 50 initial rand starting aligns? ('1' for yes or '0' for no): 1
Use 500 iters to perform the E-M steps? ('1' for yes or '0' for no): 0
Please specify a value for the number of iters: 50
Please specify the path of the fasta file: example.fasta
{'final_scores_seqs_posits_motifs': {'final_motif': 'TTATCA',
                                     'final_pos': 18,
                                     'final_score': 411.85700000000003,
                                     'final_seq': 0,
                                     'final_sum_scores': 4301.683},
 'max_scores_posits_motifs': {'max_motifs_matrix': ['TTATCA',
                                                    'CGGTCA',
                                                    'CTATCA',
                                                    'CAACCA',
                                                    'TTACAA',
                                                    'CTATCT',
                                                    'TTATCA',
                                                    'TTATCA',
                                                    'CTATAA',
                                                    'CTATCT',
                                                    'TGGTCA',
                                                    'TTGTAA',
                                                    'TTATCT',
                                                    'TTATCT',
                                                    'TTATCA',
                                                    'CTATAA',
                                                    'TTATCC',
                                                    'GGATAA',
                                                    'TGATAA',
                                                    'CTAGTA',
                                                    'GTATCC',
                                                    'GTATCC',
                                                    'CTATCT',
                                                    'CTATCC',
                                                    'TTATCT',
                                                    'CTATCG',
                                                    'TTATCA',
                                                    'CTATCT',
                                                    'TTGTCA'],
                              'max_posits_matrix': [18,
                                                    22,
                                                    15,
                                                    29,
                                                    19,
                                                    18,
                                                    20,
                                                    2,
                                                    17,
                                                    14,
                                                    21,
                                                    33,
                                                    20,
                                                    2,
                                                    10,
                                                    3,
                                                    13,
                                                    12,
                                                    16,
                                                    33,
                                                    14,
                                                    2,
                                                    23,
                                                    26,
                                                    17,
                                                    15,
                                                    19,
                                                    15,
                                                    2],
                              'max_scores_matrix': [411.85700000000003,
                                                    11.096,
                                                    298.58499999999998,
                                                    1.863,
                                                    11.616,
                                                    158.79300000000001,
                                                    411.85700000000003,
                                                    411.85700000000003,
                                                    90.510000000000005,
                                                    158.79300000000001,
                                                    15.305999999999999,
                                                    25.829999999999998,
                                                    219.03200000000001,
                                                    219.03200000000001,
                                                    411.85700000000003,
                                                    90.510000000000005,
                                                    115.36,
                                                    5.3369999999999997,
                                                    22.423999999999999,
                                                    2.085,
                                                    27.454999999999998,
                                                    27.454999999999998,
                                                    158.79300000000001,
                                                    83.632999999999996,
                                                    219.03200000000001,
                                                    35.851999999999997,
                                                    411.85700000000003,
                                                    158.79300000000001,
                                                    85.212999999999994]}}
```

## License
This project is licensed under the [GNU General Public License GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Questions?
**Shane Kimble** - shanekimble12@hotmail.com - [summonholmes](https://github.com/summonholmes)

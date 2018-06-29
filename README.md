# em-algorithm-python
This is the Expectation-Maximization (EM) algorithm performed on multiple DNA sequences.  As a latent variable model, the EM algorithm is an example of unsupervised learning.  The primary objective of this program is to discover optimal alignment positions and motifs for multiple DNA sequences.  

Optimal alignment positions and motifs are the latent variables requiring detection; meaning that they are known to exist, are known to possess high importance, but are unable to be identified directly.  

Since the starting positions of the DNA motifs are always randomized, no two runs are ever alike.  However, when provided the same input more than once, the program results are (almost) always the same.

The Expectation Step composes the vast majority of program code.  Construction of the log-likelihood function is performed matrix-by-matrix, requiring more preparation.  The Maximization Step scores on the log-likelihood function derived from the previous Expectation Step.  Using the maximum score plus the current information about optimal positions and motifs, the Expectation and Maximization Steps are repeated X number of times.  After enough iterations, convergence occurs to identify the optimal alignment positions and motifs.

A brief description of the five classes within 'em_oop/main.py':
1. EM_Input: Loads all user data including the number of random alignments, the number of iterations, and the FASTA file.  Parent of all classes.

2. EM_Core: The first iteration or 'for' loop over the number of random alignments.  The entire program can be thought of as one gigantic, quintuple, nested 'for' loop.  Therefore, this is the 'nucleas' of the program.  Child of EM_Input.

3. EM_Count: Performs all counting operations regarding bases, motifs, and normalization.  Child of EM_Core.

4. EM_Matrix: Consists of the first Expectation Step with the end goal of generating the log odds matrix.  Child of EM_Count.

5. EM_Run: Consists of remaining four inner layers of the 'for' loop from EM_Core.  Maximization is performed, scores are calculated against the log odds matrix, and the entire process repeats to convergence.  Note that the function 'exp_max_update_log_odds' is the entire Expectation Step and recycles all functions from EM_Count & EM_Matrix.  Child of EM_Matrix.

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

# em-algorithm-python
This is the Expectation-Maximization (EM) algorithm performed on multiple DNA sequences.  As a latent variable model, the EM algorithm is an example of unsupervised learning.  The primary objective of this program is to discover optimal alignment positions and motifs for multiple DNA sequences.  

Optimal alignment positions and motifs are the latent variables requiring detection; meaning that they are known to exist, are known to possess high importance, but are unable to be identified directly.  

Since the starting positions of the DNA motifs are always randomized, no two runs are ever alike.  However, when provided the same input more than once, the program results are (almost) always the same.

The Expectation Step composes the vast majority of program code.  Construction of log-likelihood is performed matrix-by-matrix, requiring more preparation.  The Maximization Step scores on the log-likelihood derived from the previous Expectation Step.  Using the maximum score and its current information about optimal positions and motifs, the Expectation and Maximization Steps are repeated X number of times.  After enough iterations, convergence occurs to identify the optimal alignment positions and motifs.

A brief description of the five classes within 'em_oop/main.py':
1. EM_Input: Loads all user data including the number of random alignments, the number of iterations, and the FASTA file.  Parent of all classes.

2. EM_Core: The first iteration or 'for' loop over the number of random alignments.  The entire program can be thought of as one gigantic, quintuple, nested 'for' loop.  Therefore, this is the 'nucleas' of the program.  Child of EM_Input.

3. EM_Count: Performs all counting operations regarding bases, motifs, and normalization.  Child of EM_Core.

4. EM_Matrix: Consists of the first Expectation Step with the end goal of generating the log odds matrix.  Child of EM_Count.

5. EM_Run: Consists of remaining four inner layers of the 'for' loop from EM_Core.  Maximization is performed, scores are calculated against the log odds matrix, and the entire process repeats to convergence.  Note that the function 'exp_max_update_log_odds' is the entire Expectation Step and recycles all functions from EM_Count & EM_Matrix.  Child of EM_Matrix.

## Getting Started
The program requires few dependences and should be trivial to set up.  However, an in-depth understanding of the EM algorithm requires some knowledge of bioinformatics, data science, and machine learning.  I'd recommend using the 'main.py' file in the 'em_oop' folder, because the process is more easily understood.  However, you may also run 'main.py' in 'em_functional' file if you prefer functional programming.

### Dependenciese:
* python3  
* python3-biopython  
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
 summonholmes@10x-Orange-G  ~/Documents/em-algorithm-python   master  conda-python main.py
Welcome to my Python implementation of Expectation-Maximization
Please specify the width of the motif: 6
Use 50 initial random starting alignments? ('1' for yes or '0' for no): 1
Use 500 iterations to perform the E-M steps? ('1' for yes or '0' for no): 1
Please specify the path of the fasta file: example.fasta

Please be patient.  Now preparing for E-M...
Progress: 100.0%
DONE!  First Horizontal Dictionary: Max scored motif (Alsofound below)
Vertical Dictionary: Max set of positions, scores, and motifs

{'max_final_sco_seq_pos_mot': {'max_final_motif': 'TTATCT',
                               'max_final_position': 20,
                               'max_final_score': 349.706,
                               'max_final_sequence': 12,
                               'sum_score_max_motif': 4265.537},
 'scores_pos_motifs': {'max_motifs': ['TTATCA',
                                      'CTGACT',
                                      'CTATCA',
                                      'ATAACT',
                                      'CTGATT',
                                      'CTATCT',
                                      'TTATCA',
                                      'TTATCA',
                                      'CTATAA',
                                      'CTATCT',
                                      'ATTTCA',
                                      'TTGTAA',
                                      'TTATCT',
                                      'TTATCT',
                                      'TTATCA',
                                      'CTATAA',
                                      'TTATCC',
                                      'AAGTAA',
                                      'TTGATA',
                                      'ATAACA',
                                      'CTGTAT',
                                      'CTGTAT',
                                      'CTATCT',
                                      'TTGTCT',
                                      'TTATCT',
                                      'TAGTCT',
                                      'TTATCA',
                                      'CTATCT',
                                      'TTGTCA'],
                       'max_pos': [18,
                                   28,
                                   15,
                                   19,
                                   15,
                                   18,
                                   20,
                                   2,
                                   17,
                                   14,
                                   31,
                                   33,
                                   20,
                                   2,
                                   10,
                                   3,
                                   13,
                                   28,
                                   15,
                                   26,
                                   12,
                                   0,
                                   23,
                                   4,
                                   17,
                                   29,
                                   19,
                                   15,
                                   2],
                       'max_scores': [278.783,
                                      25.67,
                                      173.766,
                                      17.802,
                                      4.491,
                                      217.972,
                                      278.783,
                                      278.783,
                                      56.493,
                                      217.972,
                                      9.87,
                                      55.793,
                                      349.706,
                                      349.706,
                                      278.783,
                                      56.493,
                                      36.504,
                                      1.269,
                                      5.744,
                                      14.192,
                                      43.622,
                                      43.622,
                                      217.972,
                                      215.269,
                                      349.706,
                                      18.405,
                                      278.783,
                                      217.972,
                                      171.611]}}
```

## License
This project is licensed under the [GNU General Public License GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Questions?
**Shane Kimble** - shanekimble12@hotmail.com - [summonholmes](https://github.com/summonholmes)

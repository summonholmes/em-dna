# em-dna
This is the Expectation-Maximization (EM) algorithm performed on multiple DNA sequences.  As a latent variable model, the EM algorithm is an example of unsupervised learning.  The primary objectives of this program are to discover optimal alignment positions and motifs for multiple DNA sequences; and to provide a clear, efficient EM implementation.  

Optimal alignment positions and motifs are the latent variables requiring detection; meaning that they are known to exist, are known to possess high importance, but are unable to be identified directly.

Since the starting positions of the DNA motifs are always randomized, no two runs are ever alike.  However, when provided the same input more than once, the program results are (almost) always the same.

The Expectation (E) step composes the majority of source code.  Total positional counts, frequencies, odds, then log-odds calculations compose the E step and derive the log-likelihood function.  The Maximization Step scores on the log-likelihood function derived from the previous Expectation Step.  Taking the score for every possible contiguous motif, for each sequence, the motif positions that generate the maximum scores on each sequence are identified, and become the input of the next E step.  The E and M steps are then repeated X number of times.  After enough iterations, convergence occurs to identify the optimal alignment positions and motifs.

An implementation in Julia is also available.  The Julia implementation was constructed to test the claims regarding Julia's superior performance.  Despite being as identical to the Python code blueprint as possible (minus the object-oriented structure), the Julia code runs far slower.  Perhaps I can improve it when I learn more about Julia's best practices, or Julia isn't what it claims to be.

A brief description of the six classes:
1. EM_Input: Loads all user data including the number of random alignments, the number of iterations, and the FASTA file.  This class has an interactive mode, but is commented out in favor of default parameters.  The defaults include a motif width of 6, 50 total random alignments, and 50 iterations.

2. EM_Prep: Consolidates all consistent and reusable information.  This includes all bases as a string, all base counts, all contiguous motifs, and cumulative sums.

2. EM_Core: The first iteration or 'for' loop over the number of random alignments.  Previously, the entire program could be thought of as one gigantic, quintuple, nested 'for' loop.  However, this approach has been reduced to a double 'for' loop via numpy vectorization and indexing.  This is the 'nucleas' of the program, which encompasses the entire EM process and iterates for the specified number of total random alignments.

3. EM_Count: Performs all counting operations regarding bases, motifs, and normalization.  The only area where looping occurs is over the motif width when storing positional counts.  The rest of this class aims to be as Pythonic and vectorized as possible.

4. EM_Matrix: Consists of the E step, and generates the log odds matrix for the M step.  Numpy matrices are preferred here for their vector capabilities.  Child of EM_Count.

5. EM_Run: Consists of remaining inner two 'for' loops from mentioned above in EM_Core, and the M step.  Scores are calculated using a numpy character array for each motif and the log odds matrix.  The entire process repeats to convergence, and the final results are printed in 'main.py'.  Note that the function 'update_log_odds' is the entire Expectation Step and recycles all functions from EM_Count & EM_Matrix.  Child of EM_Matrix.

## Getting Started
This program requires few dependences and should be trivial to set up.  However, an in-depth understanding of the EM algorithm requires some knowledge of bioinformatics, data science, and machine learning.

### Dependencies:
* python3  
* python3-numpy
* python3-pprint  
* python3-termcolor

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
![alt text](https://raw.githubusercontent.com/summonholmes/em-dna/master/example.png)
```
Results:
{'final_motif': 'TTATCT',
 'final_pos': 20,
 'final_score': 401.64821219980882,
 'final_seq': 12,
 'final_sum_scores': 4852.1552401637646,
 'motifs_set': ['TTATCA',
                'CTGACT',
                'CTATCA',
                'ATAACT',
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
                'AGATAT',
                'TGATAA',
                'ATAACA',
                'CTGTAT',
                'CTGTAT',
                'CTATCT',
                'TTGTCT',
                'TTATCT',
                'CTATCG',
                'TTATCA',
                'CTATCT',
                'TTGTCA'],
 'posits_set': [18,
                28,
                15,
                19,
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
                20,
                16,
                26,
                12,
                0,
                23,
                4,
                17,
                15,
                19,
                15,
                2],
 'scores_set': [349.55201551685258,
                11.048262627379067,
                208.43428487091478,
                10.19243588705069,
                8.934927563323658,
                239.49871310500956,
                349.55201551685258,
                349.55201551685258,
                87.115543742405634,
                239.49871310500956,
                16.479781956304635,
                53.931154636101589,
                401.64821219980882,
                401.64821219980882,
                349.55201551685258,
                87.115543742405634,
                45.618802496192338,
                4.3536481630627533,
                18.658492127411794,
                8.8704154509530913,
                36.951416551025055,
                36.951416551025055,
                239.49871310500956,
                148.26793906959458,
                401.64821219980882,
                29.494915407020887,
                349.55201551685258,
                239.49871310500956,
                129.03669271786535]}
```

## License
This project is licensed under the [GNU General Public License GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Questions?
**Shane Kimble** - shanekimble12@hotmail.com - [summonholmes](https://github.com/summonholmes)

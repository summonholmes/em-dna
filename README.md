# em-dna
This is the Expectation-Maximization (EM) algorithm performed on multiple DNA sequences.  As a latent variable model, the EM algorithm is an example of unsupervised learning.  The primary objectives of this program are to discover optimal alignment positions and motifs for multiple DNA sequences; and to provide a clear, efficient EM implementation.  

Optimal alignment positions and motifs are the latent variables requiring detection; meaning that they are known to exist, are known to possess high importance, but are unable to be identified directly.

Since the starting positions of the DNA motifs are always randomized, no two runs are ever alike.  However, when provided the same input more than once, the program results are (almost) always the same.

The Expectation (E) step composes the majority of source code.  Total positional counts, frequencies, odds, then log-odds calculations compose the E step and derive the log-likelihood function.  The Maximization Step scores on the log-likelihood function derived from the previous Expectation Step.  Taking the score for every possible contiguous motif, for each sequence, the motif positions that generate the maximum scores on each sequence are identified, then become the input of the next E step.  The E and M steps are then repeated X number of times.  After enough iterations, convergence occurs to identify the optimal alignment positions and motifs. 

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
* python3-pandas
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
![alt text](https://raw.githubusercontent.com/summonholmes/em-dna/master/Images/example.png)
```
    Final Scores  Final Positions Final Motifs
0     303.827137               18       TTATCA
1       4.534022               28       CTGACT
2     230.734903               15       CTATCA
3      15.733361               19       ATAACT
4       4.545193               19       TTACAA
5     263.652353               18       CTATCT
6     303.827137               20       TTATCA
7     303.827137                2       TTATCA
8      33.331416               17       CTATAA
9     263.652353               14       CTATCT
10     20.255142               31       ATTTCA
11     20.255142               13       ATTTCA
12    347.172181               20       TTATCT
13    347.172181                2       TTATCT
14    303.827137               10       TTATCA
15     51.522972               25       GTATCC
16    143.810374               13       TTATCC
17      1.287469               18       GTAGAT
18      5.407760               29       CTTCCT
19     13.769024               26       ATAACA
20     51.522972               14       GTATCC
21     51.522972                2       GTATCC
22    263.652353               23       CTATCT
23    109.213657               26       CTATCC
24    347.172181               17       TTATCT
25     37.783513               15       CTATCG
26    303.827137               19       TTATCA
27    263.652353               15       CTATCT
28     38.807012                2       TTGTCA
('Final Max Score', 347.1721810115789)
('Final Sum Scores', 4449.298542515307)
('Final Sequence', 12)
('Final Position', 20)
('Final Motif', 'TTATCT')
```

## License
This project is licensed under the [GNU General Public License GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Questions?
**Shane Kimble** - shanekimble12@hotmail.com - [summonholmes](https://github.com/summonholmes)

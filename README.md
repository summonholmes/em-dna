# em-algorithm-python
This is the Expectation-Maximization (EM) machine learning algorithm performed on multiple DNA sequences, with the intention of finding the optimal alignment position.  The program will ("almost") always provide the same final output, although no two runs are ever alike.  This is because the starting positions of the DNA motif length are always randomized.  After enough iterations, convergence occurs to identify the optimal set of scores, position, and motifs.

## Getting Started
The program requires few dependences and should be trivial to set up.  However, an in-depth understanding of the EM algorithm requires some knowledge of bioinformatics, data science, and machine learning.

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

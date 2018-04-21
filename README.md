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

Now preparing for E-M...
M = Max Motif
S = Max Score
P = Max Position
# = Max Sequence #
SS = Max Sum Scores

Please be patient...

DONE!  First horizontal array: max information.
Vertical arrays: positions, scores, and motifs.

MAX: Motif      Score    Pos    Seq     SumScore
[['TTATCT', 390.722, 20, 12, 4867.459],
 [[18,
   28,
   15,
   19,
   19,
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
   22,
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
  [362.541,
   9.741,
   227.859,
   12.96,
   8.179,
   245.572,
   362.541,
   362.541,
   84.566,
   245.572,
   10.593,
   43.231,
   390.722,
   390.722,
   362.541,
   84.566,
   47.505,
   2.393,
   7.923,
   12.025,
   29.283,
   29.283,
   245.572,
   125.54,
   390.722,
   48.168,
   362.541,
   245.572,
   116.485],
  ['TTATCA',
   'CTGACT',
   'CTATCA',
   'ATAACT',
   'TTACAA',
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
   'ATATTG',
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
   'TTGTCA']]]
```

## License
This project is licensed under the [GNU General Public License GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Questions?
**Shane Kimble** - shanekimble12@hotmail.com - [summonholmes](https://github.com/summonholmes)

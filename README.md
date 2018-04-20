# em-algorithm-python
Demonstration of the Expectaeion Maximization (EM) algorithm to compare multiple DNA sequences read from a FASTA file.
This is a machine learning algorithm performed on multiple DNA sequences to find the optimal alignment position.  The program will ("almost") always provide the same final output, although each run is unique.  This is because the starting positions of the DNA motif length are always randomized.  After enough iterations, convergence occurs to identify the optimal set of scores, position, and motifs.

## Getting Started
The program requires few dependences, and should be trivial to set up.  However, an in-depth understanding of the EM algorithm requires some knowledge of bioinformatics, data science, and machine learning.

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
Welcome to a Python implementation of Expectation-Maximization
Please specify the width of the motif: 6
```
```
Use 50 initial random starting alignments? ('1' for yes or '0' for no): 1
```
```
Use 500 iterations to perform the E-M steps? ('1' for yes or '0' for no): 0
Please specify a value for the number of iterations: 50
```
```
Please specify the path of the fasta file: example.fasta
```
```
Now preparing for E-M...
M = Max Motif
S = Max Score
P = Max Position
# = Max Sequence #
SS = Max Sum Scores

Please be patient...
```
```
DONE!  First horizontal array: max information.
Vertical arrays: positions, scores, and motifs.

MAX: Motif      Score    Pos    Seq     SumScore
[['TTATCA', 229.603, 18, 0, 2505.927],
 [[18,
   33,
   15,
   ...
```

## License
This project is licensed under the [GNU General Public License GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Questions?
**Shane Kimble** - shanekimble12@hotmail.com - [summonholmes](https://github.com/summonholmes)

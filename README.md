# FYS3150_Project_2

This program consists of 10 functions in addition to the main function. We first have
three functions that deal with arrays.

double ** sparse( int n)
double* zeros(int n)
void clear_memory(double ** & A, int n)

The function "sparse" returns an array of dimension n times n. 
The function "zeros" returns an array of dimension n.
The function "clear_memory" deletes an n times n array from memory.

Then there are 3 functions that take a matrix and do something with it.
We have:

void get_max(double** &A, int n, int &k, int &l, double &maxval) 
void print(double** & A, int n)
void jacobi_rot_algorithm(double ** & A, double ** & U, double * & lambda, int N)

The function "get_max" sets the integer maxval equal to the pivot element (entry with greatest absolute value) of A.
Moreover, it sets k and l equal to the row number and column number corresponding to the pivot element.
  The function "print" simply prints the array A of dimension n times n to the terminal.
The function "jacobi_rot_algorithm" finds eigenvectors and associated eigenvalues for a given matrix A.
It then stores the eigenvectors in the matrix U as column vectors, and the eigenvalues are stored in the array lambda.

There are two unit tests:

void unit_test1(void)
void unit_test2(void)

"unit_test1" checks the validity of the function "get_max" and prints the result to the terminal.,
"unit_test2" checks the validity of "jacobi_rot_algorithm".

also there is a function "count_iterations" which sets up a problem for "jacobi_rot_algorithm" to solve,
for an n times n matrix A determined by a given n.

Lastly, there is a function 

void comparison(int N)

This guy compares our implementation of Jacobi's algorithm to the eigenvalue/eigenvector solver of armadillo
in terms of speed. The results are then stored as files with filenames of the type "jacobi_vs_armadillo_N_equals_100".


The main program consists of three parts, separated by lines of equalsigns (============).
The first part compares the tridiagonal toeplitz matrix eigenvalues found using "jacobi_rot_algorithm" to the analytical ones.
The second part finds the 10 lowest eigenvalues of a dimensionless schrödingers equation with harmonic oscillator potential,
and the last part finds the 10 lowest eigenvalues of a dimensionless schrödingers with two particles with variable oscillation.
The results of each part is printed to the terminal.


For example, the output of the unaltered program becomes:

"Relative error of eigenvalues of buckling beam (N = 10)

analytical:     numerical:      relative error:
38.196601125015363      38.196601125010609      1.2444905468622788e-13
82.442949541515418      82.442949541509336      7.3775208819326013e-14
138.19660112502623      138.19660112501043      1.1434775034017053e-13
200.00000000002066      199.99999999999682      1.192290710605325e-13
261.80339887501304      261.80339887500139      4.4510120634563753e-14
317.55705045851795      317.55705045848185      1.1366641340349709e-13
361.8033988750089       361.80339887498963      5.326074617798446e-14
390.2113032590421       390.21130325903056      2.9571706232927114e-14
1.6911946456991202e-306 4.2279662430882073e-307 0.75000120454257857

Solving Schr÷dingers for harmonic oscillator potential with (N=150 and infinity = 10)
The 10 smalles eigenvalues (energies):

2.8516915643062219
6.771190294406777
10.705018895874323
14.643823583845606
18.584072728105806
22.523993232766585
26.462546401104049
30.399061157079519
34.333073358429353
38.264245383984552

Solving Schr÷dingers for 2-particles with (N=150 and infinity = 10 and omega_r=5 )
The 10 smalles eigenvalues (energies):

20.507780622666342
41.688408495639848
62.488732479579838
83.021971591830038
103.33637881101346
123.4577494244657
143.40160866355285
163.17798372678314
182.79363608299028
202.25323608802762

success!"

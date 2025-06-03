# Grover adaptive serach based multivariable boolean function optimization

We implement the grover adaptiver search based optimization of interger-valued boolean functions in $n$ varialbes.

The code and its complexity analysis is orgainzes as follows:

  1. The file "QC_MP1_gates.py" contains implementation of all the basic gates used. The gates are implemented only using mcx and $1$-qubit gates: Hadamard and the phase gates, and possible ancillas in some cases.
  2. The file "QC_MP1_marker.py" contains the implementation of the marker-oracle and the function-encoding routines. The function-encoding routine encodes a multivriable boolean function in a quantum circuit. The marker-oracle routine returns the unitary $U_{f,t,C}$ for a given objective function $f:\mathbb{F}_2^n\to \mathbb{Z}$, a threshold $t\in\mathbb{Z}$ and a constraint function $C:\mathbb{F}_2^n\to \mathbb{Z}$.
  3. The file "QC Mini Project 1.ipynb" demonstrates the marker-oracle function $U_{f,t,C}$ for various choices of $n$, $f$, $t$, $C$, and the input $x$ to $f$ and $C$.
  4. The file "QC MP1 Analysis of Complexity.pdf" containes a detailed description of the design and the analysis of space and time complexity of the ciruit. 

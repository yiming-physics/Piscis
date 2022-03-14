# Piscis 

## Introduction

<font size=4>**Pisicis**</font> [’bi:ʃıs] is a  simple Python library for quantum computing based on state vector method. If you have requirements of speed, I recommend you to use other mature toolkits such as QuEST (C/C++), CuQuantum (C/C++) and Qiskit (Python). However, for the purpose of a quick and convenient verification of quantum algorithms with few qubits, Piscis is useful, and this is the motivation for me to write this toolkit.

## Features

In Pisics, after you create a quantum circuit and add some 
quantum gates, 
you can print the state directly by 

```Python
my_circuit.print_state()
```



For example, the following codes 

```Python
from piscis import QCirc

my_circuit1 = QCirc(qnum=3)

my_circuit1.H(0)
my_circuit1.X(2)
my_circuit1.CNOT(0, 1)

my_circuit1.print_state()
```
output 

```angular2html
(0.7071)|100〉+ (0.7071)|111〉
```

You can also print the state distribution by
```Python
my_circuit1.print_state_distribution()
```
to get
```angular2html
prob		|	state	
0.50000		|	|100〉
0.50000		|	|111〉
```

In Piscis, you can use the non-unitary projectors $|0\rangle\langle0|$ and $|1\rangle\langle1|$.

For example, the following codes 
```Python
my_circuit2 = QCirc(qnum=2)

my_circuit2.H(0)
my_circuit2.projector_to_one(0)

my_circuit2.print_state()
```
output
```angular2html
(1.0)|01〉
```

In Picis, you can get the probability or amplitude (up to a global phase factor) of a specified qubit 
to be zero or one.

For example, the following codes
```Python
from math import pi

my_circuit3 = QCirc(qnum=4)
my_circuit3.Rx(2, pi/7)

prob = my_circuit3.get_prob_of_zero(2)
print("The probability of the 2nd qubit to be 0 is {}.".format(prob))
```
output
```angular2html
The probability of the 2nd qubit to be 0 is 0.95048443.
```


Standard quantum computing process (measure and read) and basic 
quantum gates are also available in Piscis. 

For more functions, see APIs.

## APIs
See *APIs.md*.

## Contact

If you find any bug, or have any question/advice, you may 
email me by **dym0107@163.com**.


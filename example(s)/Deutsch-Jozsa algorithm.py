# Deutsch-Jozsa Algorithm with Piscis
import piscis as ps
import numpy as np

def main():
    # Create a circuit with 4 qubits.
    qubit_num = 4
    my_circuit = ps.QCirc(qubit_num)

    # Assign the first register.
    reg1 = range(0, 3)

    # DIY a quantum gate by unitary X.
    X = np.array([[0, 1], [1, 0]])
    my_circuit.U_2by2(X, 3)

    # Apply Hadamard gates to all qubits.
    my_circuit.H_all()

    # f(x) in this demo
    my_circuit.CNOT(0, 3)       # '0' is the controlled while '3' is the target.
    my_circuit.CX(1, 3)       # CX and CNOT are interchangeable.
    my_circuit.CU_2by2(X, 2, 3)  # DIY a controlled unitary.

    # Apply Hadamard gates to register 1.
    my_circuit.H_multi(reg1)

    # Depict the state distribution.
    my_circuit.depict_state_distribution()

    # Measure all qubits in register 1.
    my_circuit.measure_multi(reg1)

    # Judge whether f(x) is constant.
    prob = 1
    for i in reg1:
        prob *= my_circuit.get_prob_of_zero(i)
    print("f is constant.") if prob == 1 else print("f is balanced.")

if __name__ == "__main__":
    main()

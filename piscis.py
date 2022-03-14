''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
* Source Code of Piscis V-0.0.0.1
* Author: @Yiming_Ding
* Email: dym0107@163.com
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
import cmath
from math import pi, cos, sin, ceil, floor
import numpy as np
import random
import matplotlib.pyplot as plt

'''''''''''''''''''''''''''''''''''''''''
Methods in 'QCirc':
    * Basic Modules (hidden)
    * State Initialization
    * Statistical Tools
    * Non-unitary Operations (Projectors & Measurements)
    * 1-Qubit Unitary Gates
    * 2-Qubit Unitary Gates
    * 3-Qubit Unitary Gates
    * Multi-Controlled Unitary Gates
    * Error Reporting (hidden)
'''''''''''''''''''''''''''''''''''''''''

class QCirc():
    def __init__(self, qnum):
        self.__positive_integer_check(qnum)
        self.__qnum = qnum
        self.__dim = 2**qnum
        self.__state_vec = None
        if self.__state_vec == None:
            self.__state_vec = (np.zeros((self.__dim, 1))).astype(complex)
            self.__state_vec[0, 0] = 1
        self.__measured_record = np.zeros(self.__dim)

    # Accuracy control
    @property
    def __round_to(self):
        return [4, 8]

    '''***********************************************************************
    * Basic Modules
    ***********************************************************************'''

    # Reuturn a Hermitian form of an unitary matrix U
    def __get_U_dagger(self, U):
        U_dagger = np.conj(U).T
        return U_dagger

    # Check whether matrix U is unitary.
    def __check_unitary(self, U):
        m, n = len(U), len(U[0])
        if (m != n):
            self.__report_error_input_U(0)
        else:
            I = np.eye(m)
            U_dagger = self.__get_U_dagger(U)
            if (np.allclose(I, np.matmul(U, U_dagger)) == False):
                self.__report_error_input_U(0)

    # For a specified label, return the value of this bit (little-endian).
    def __get_specified_bin(self, N, bin_len, label):
        _N = N
        search_bin_state = []
        for j in range(bin_len):
            search_bin_state.append(_N % 2)
            _N = _N // 2
        sp_bin = search_bin_state[label]
        return sp_bin

    # For an integer N, print its binary form in Dirac notation.
    def __print_binary_state(self, N, bin_len):
        _N = N
        bin_state = []
        print("|", end="")
        for j in range(bin_len):
            bin_state.append(_N % 2)
            _N = _N // 2
        for j in range(bin_len):
            print(bin_state[bin_len - 1 - j], end="")
        print("〉", end="")

    '''***********************************************************************
    * State Initialization
    ***********************************************************************'''

    # Normalize a state.
    def __state_normalization(self):
        sum = 0
        for j in range(len(self.__state_vec)):
            sum += self.__state_vec[j] * self.__state_vec[j].conjugate()
        if sum.imag != 0:
            raise Exception("invalid state vector.")
        factor = 1 / (sum.real ** 0.5)
        self.__state_vec = factor * self.__state_vec
        return self.__state_vec

    # Initialize a state by a set of specified amplitudes.
    def ini_set_amp(self, amp_list):
        self.__set_amp_error(len(amp_list))
        amp_list = amp_list.astype(complex)
        for i in range(self.__dim):
            self.__state_vec[i, 0] = amp_list[i]
        self.__state_vec = self.__state_normalization()

    '''***********************************************************************
    * Statistical Tools
    ***********************************************************************'''

    # Print the state, except for a global phase factor.
    def print_state(self):
        plus_sig, N = 0, 0
        for i in range(self.__dim):
            if (np.allclose(self.__state_vec[i, 0].real, 0) == True) and (np.allclose(self.__state_vec[i, 0].imag, 0) == True):
                pass
            else:
                if plus_sig != 0:
                    print("+ ", end="")
                if self.__state_vec[i, 0].real != 0:
                    print("({}".format(round(self.__state_vec[i, 0].real, self.__round_to[0])), end="")
                else:
                    print("(", end="")
                if self.__state_vec[i, 0].imag > 0:
                    print("+{}i".format(round(self.__state_vec[i, 0].imag, self.__round_to[0])), end="")
                elif self.__state_vec[i, 0].imag < 0:
                    print("{}i".format(round(self.__state_vec[i, 0].imag, self.__round_to[0])), end="")
                print(")", end="")
                self.__print_binary_state(N, self.__qnum)
                plus_sig = 1
            N += 1
        print("")

    # Depict the state distribution.
    def depict_state_distribution(self):
        self.__depict_error()
        figure_scale = 10
        plt.figure(figsize=(figure_scale, figure_scale), dpi=100)
        prob_set = []
        pick_from_dim = []
        plt.title("State Distribution", fontsize=30)
        for j in range(self.__dim):
            if (np.allclose(self.__state_vec[j, 0].real, 0) == True) and (np.allclose(self.__state_vec[j, 0].imag, 0) == True):
                pass
            else:
                pick_from_dim.append(j)
                prob_set.append((self.__state_vec[j, 0] * self.__state_vec[j, 0].conjugate()).real)

        plt.bar(pick_from_dim, prob_set, align='center', color='steelblue', alpha=0.8)
        # plt.xticks(pick_from_dim, range(self.__dim))
        plt.xticks(pick_from_dim, pick_from_dim)
        plt.tick_params(labelsize=25)
        plt.ylim([0, 1])
        plt.show()

    # Print the state distribution.
    def print_state_distribution(self):
        print("\nprob\t\t|\tstate\t")
        for j in range(self.__dim):
            if (np.allclose(self.__state_vec[j, 0].real, 0) == True) and (np.allclose(self.__state_vec[j, 0].imag, 0) == True):
                pass
            else:
                print("%.5f" % (self.__state_vec[j, 0] * self.__state_vec[j, 0].conjugate()).real, end="\t")
                print("\t|\t", end="")
                self.__print_binary_state(j, self.__qnum)
                print("")

    # Output the circuit after measurements.
    def output(self):
        s = "|"
        for j in range(self.__dim):
            if self.__measured_record[j] == 1:
                s += "1" if self.get_prob_of_zero(j) == 0 else "0"
        s += "〉"
        self.__output_error(len(s))
        print(s)

    '''***********************************************************************
    * Non-unitary Operations (Projectors & Measurements)
    ***********************************************************************'''

    # Project a specified qubit to |1〉.
    def projector_to_one(self, label):
        self.__test_error_label(label)
        N = 0
        for j in range(self.__dim):
            sp_bin_for_N = self.__get_specified_bin(N, self.__qnum, label)
            if sp_bin_for_N != 1:
                self.__state_vec[j, 0] = 0
            N += 1
        self.__state_vec = self.__state_normalization()

    # Project a specified qubit to |0〉.
    def projector_to_zero(self, label):
        self.__test_error_label(label)
        N = 0
        for j in range(self.__dim):
            sp_bin_for_N = self.__get_specified_bin(N, self.__qnum, label)
            if sp_bin_for_N != 0:
                self.__state_vec[j, 0] = 0
            N += 1
        self.__state_vec = self.__state_normalization()

    # For a specified label, return the probability of this bit to be |1〉.
    def get_prob_of_one(self, label):
        self.__test_error_label(label)
        N = 0
        prob = 0
        for j in range(self.__dim):
            sp_bin_for_N = self.__get_specified_bin(N, self.__qnum, label)
            if sp_bin_for_N == 1:
                prob += self.__state_vec[j, 0] * self.__state_vec[j, 0].conjugate()
            N += 1
        return round(prob.real, self.__round_to[1])

    # For a specified label, return the amplitude of this bit to be |1〉.
    def get_amp_of_one(self, label):
        self.__test_error_label(label)
        N = 0
        amplitude = 0
        for j in range(self.__dim):
            sp_bin_for_N = self.__get_specified_bin(N, self.__qnum, label)
            if sp_bin_for_N == 1:
                amplitude += self.__state_vec[j, 0]
            N += 1
        return round(amplitude.real, self.__round_to[1]) + 1j * round(amplitude.imag, self.__round_to[1])

    # For a specified label, return the probability of this bit to be |1〉.
    def get_prob_of_zero(self, label):
        self.__test_error_label(label)
        N = 0
        prob = 0
        for j in range(self.__dim):
            sp_bin_for_N = self.__get_specified_bin(N, self.__qnum, label)
            if sp_bin_for_N == 0:
                prob += self.__state_vec[j, 0] * self.__state_vec[j, 0].conjugate()
            N += 1
        return round(prob.real, self.__round_to[1])

    # For a specified label, return the amplitude of this bit to be |0〉.
    def get_amp_of_zero(self, label):
        self.__test_error_label(label)
        N = 0
        amplitude = 0
        for j in range(self.__dim):
            sp_bin_for_N = self.__get_specified_bin(N, self.__qnum, label)
            if sp_bin_for_N == 0:
                amplitude += self.__state_vec[j, 0]
            N += 1
        return round(amplitude.real, self.__round_to[1]) + 1j * round(amplitude.imag, self.__round_to[1])

    # Add a measurement to a qubit on the circuit.
    def measure(self, label):
        self.__measured_record[label] = 1
        sample_num = 1000
        zero_num = int(self.get_prob_of_zero(label) * sample_num)
        rand = random.randint(1, sample_num)
        if rand <= zero_num:
            self.projector_to_zero(label)
        else:
            self.projector_to_one(label)

    # Add several measurements to a set of qubits on the circuit.
    def measure_multi(self, label_list):
        for j in label_list:
            self.measure(j)

    # Add measurements to all qubits on the circuit
    def measure_all(self):
        for i in range(self.__qnum):
            self.measure(i)

    '''***********************************************************************
    * 1-Qubit Unitary Gates
    ***********************************************************************'''

    # General definition of a single qubit operation.
    def __single_qubit_operation(self, u, label):
        self.__test_error_label(label)
        label = self.__qnum - 1 - label
        U = 1
        for j in range(self.__qnum):
            if j == label:
                U = np.kron(U, u)
            else:
                U = np.kron(U, np.eye(2))
        self.__state_vec = np.dot(U, self.__state_vec)
        return self.__state_vec

    # Pauli X operator.
    def X(self, label):
        x = np.array([[0, 1], [1, 0]])
        self.__single_qubit_operation(x, label)

    # Pauli X operators on several qubits.
    def X_multi(self, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.X(label)

    # Pauli X operators on all qubits.
    def X_all(self):
        for label in range(self.__qnum):
            self.X(label)

    # Pauli Y operator.
    def Y(self, label):
        y = np.array([[0, -1j], [1j, 0]])
        self.__single_qubit_operation(y, label)

    # Pauli Y operators on several qubits.
    def Y_multi(self, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.Y(label)

    # Pauli Y operators on all qubits.
    def Y_all(self):
        for label in range(self.__qnum):
            self.Y(label)

    # Pauli Z operator.
    def Z(self, label):
        z = np.array([[1, 0], [0, -1]])
        self.__single_qubit_operation(z, label)

    # Pauli Z operators on several qubits.
    def Z_multi(self, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.Z(label)

    # Pauli Z operators on all qubits.
    def Z_all(self):
        for label in range(self.__qnum):
            self.Z(label)

    # Hadamard gate.
    def H(self, label):
        h = np.array([[1 / 2 ** 0.5, 1 / 2 ** 0.5], [1 / 2 ** 0.5, -1 / 2 ** 0.5]])
        self.__single_qubit_operation(h, label)

    # Hadamard gates on several qubits.
    def H_multi(self, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.H(label)

    # Hadamard gates on all qubits.
    def H_all(self):
        for label in range(self.__qnum):
            self.H(label)

    # S gate.
    def S(self, label):
        s = np.array([[1, 0], [0, 1j]])
        self.__single_qubit_operation(s, label)

    # S gates on several qubits.
    def S_multi(self, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.S(label)

    # S gates on all qubits.
    def S_all(self):
        for label in range(self.__qnum):
            self.S(label)

    # T gate.
    def T(self, label):
        t = np.array([[1, 0], [0, cmath.exp(1j * pi / 4)]])
        self.__single_qubit_operation(t, label)

    # T gates on several qubits.
    def T_multi(self, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.T(label)

    # T gates on all qubits.
    def T_all(self):
        for label in range(self.__qnum):
            self.T(label)

    # Rotation X gate.
    def Rx(self, label, theta):
        rx = np.array([[cos(theta / 2), -1j * sin(theta / 2)], [1j * sin(theta / 2), cos(theta) / 2]])
        self.__single_qubit_operation(rx, label)

    # Rotation X gates on several qubits.
    def Rx_multi(self, label_list, theta):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.Rx(label, theta)

    # Rotation X gates on all qubits.
    def Rx_all(self, theta):
        for label in range(self.__qnum):
            self.Rx(label, theta)

    # Rotation Y gate.
    def Ry(self, label, theta):
        ry = np.array([[cos(theta / 2), -sin(theta / 2)], [sin(theta / 2), cos(theta / 2)]])
        self.__single_qubit_operation(ry, label)

    # Rotation Y gates on several qubits.
    def Ry_multi(self, label_list, theta):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.Ry(label, theta)

    # Rotation Y gates on all qubit.
    def Ry_all(self, theta):
        for label in range(self.__qnum):
            self.Ry(label, theta)

    # Rotation Z gate.
    def Rz(self, label, theta):
        rz = np.array([[1, 0], [0, cmath.exp(1j * theta)]])
        self.__single_qubit_operation(rz, label)

    # Rotation Z gates on several qubits.
    def Rz_multi(self, label_list, theta):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.Rz(label, theta)

    # Rotation Z gates on all qubits.
    def Rz_all(self, theta):
        for label in range(self.__qnum):
            self.Rz(label, theta)

    # U gate (DIY).
    def U_2by2(self, U, label):
        m, n = len(U), len(U[0])
        if ((m != 2) or (n != 2)):
            self.__report_error_input_U(2)
        self.__check_unitary(U)
        self.__single_qubit_operation(U, label)

    # U gates (DIY) on several qubits.
    def U_2by2_multi(self, U, label_list):
        self.__test_error_multi_label_list(label_list)
        for label in label_list:
            self.U_2by2(U, label)

    # U gates (DIY) on all qubits.
    def U_2by2_all(self, U):
        for label in range(self.__qnum):
            self.U_2by2(U, label)

    '''***********************************************************************
    * 2-Qubit Unitary Gates
    ***********************************************************************'''

    # General definition of a controlled operation.
    def __controlled_operation(self, u, ctr_label, tgt_label):
        self.__test_error_label_CU(ctr_label, tgt_label)

        ctr_label = self.__qnum - 1 - ctr_label
        tgt_label = self.__qnum - 1 - tgt_label

        matrix1 = np.array([[1, 0], [0, 0]])
        matrix2 = np.array([[0, 0], [0, 1]])

        U = np.zeros((self.__dim, self.__dim)).astype(complex)

        u1, u2 = 1, 1
        for i in range(self.__qnum):
            if i == ctr_label:
                u1 = np.kron(u1, matrix1)
            else:
                u1 = np.kron(u1, np.eye(2))
        U += u1

        for j in range(self.__qnum):
            if j == ctr_label:
                u2 = np.kron(u2, matrix2)
            elif j == tgt_label:
                u2 = np.kron(u2, u)
            else:
                u2 = np.kron(u2, np.eye(2))
        U += u2

        self.__state_vec = np.dot(U, self.__state_vec)
        return self.__state_vec

    # CNOT/Controlled Pauli X gate.
    def CNOT(self, ctr_label, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__controlled_operation(x, ctr_label, tgt_label)

    def CX(self, ctr_label, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__controlled_operation(x, ctr_label, tgt_label)

    # CNOT/Controlled Pauli X gates on several qubits.
    def CNOT_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CNOT(ctr_label, tgt_label)
    def CX_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CX(ctr_label, tgt_label)

    # Controlled Pauli Y gate.
    def CY(self, ctr_label, tgt_label):
        y = np.array([[0, -1j], [1j, 0]])
        self.__controlled_operation(y, ctr_label, tgt_label)

    # Controlled Pauli Y gates on several qubits.
    def CY_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CY(ctr_label, tgt_label)

    # Controlled Pauli Z gate.
    def CZ(self, ctr_label, tgt_label):
        z = np.array([[1, 0], [0, -1]])
        self.__controlled_operation(z, ctr_label, tgt_label)

    # Controlled Pauli Z gates on several qubits.
    def CZ_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CZ(ctr_label, tgt_label)

    # Controlled Hadamard gate.
    def CH(self, ctr_label, tgt_label):
        h = np.array([[1 / 2 ** 0.5, 1 / 2 ** 0.5], [1 / 2 ** 0.5, -1 / 2 ** 0.5]])
        self.__controlled_operation(h, ctr_label, tgt_label)

    # Controlled Hadamard gates on several qubits.
    def CH_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CH(ctr_label, tgt_label)

    # Controlled S gate.
    def CS(self, ctr_label, tgt_label):
        s = np.array([[1, 0], [0, 1j]])
        self.__controlled_operation(s, ctr_label, tgt_label)

    # Controlled S gates on several qubits.
    def CS_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CS(ctr_label, tgt_label)

    # Controlled T gate.
    def CT(self, ctr_label, tgt_label):
        t = np.array([[1, 0], [0, cmath.exp(1j * pi / 4)]])
        self.__controlled_operation(t, ctr_label, tgt_label)

    # Controlled T gates on several qubits.
    def CT_multi(self, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CT(ctr_label, tgt_label)

    # Controlled rotation X gate.
    def CRx(self, ctr_label, tgt_label, theta):
        rx = np.array([[cos(theta / 2), -1j * sin(theta / 2)], [1j * sin(theta / 2), cos(theta) / 2]])
        self.__controlled_operation(rx, ctr_label, tgt_label)

    # Controlled rotation X gates on several qubits.
    def CRx_multi(self, ctr_label, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CRx(ctr_label, tgt_label, theta)

    # Controlled rotation Y gate.
    def CRy(self, ctr_label, tgt_label, theta):
        ry = np.array([[cos(theta / 2), -sin(theta / 2)], [sin(theta / 2), cos(theta / 2)]])
        self.__controlled_operation(ry, ctr_label, tgt_label)

    # Controlled rotation X gates on several qubits.
    def CRy_multi(self, ctr_label, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CRy(ctr_label, tgt_label, theta)

    # Controlled rotation Z gate.
    def CRz(self, ctr_label, tgt_label, theta):
        rz = np.array([[1, 0], [0, cmath.exp(1j * theta)]])
        self.__controlled_operation(rz, ctr_label, tgt_label)

    # Controlled rotation Z gates on several qubits.
    def CRz_multi(self, ctr_label, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CRz(ctr_label, tgt_label, theta)

    # Controlled U gate (DIY).
    def CU_2by2(self, U, ctr_label, tgt_label):
        m, n = len(U), len(U[0])
        if ((m != 2) or (n != 2)):
            self.__report_error_input_U(2)
        self.__check_unitary(U)
        self.__controlled_operation(U, ctr_label, tgt_label)

    # Controlled U gates (DIY) on several qubits.
    def CU_2by2_multi(self, U, ctr_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CU_2by2(U, ctr_label, tgt_label)

    # SWAP gate.
    def SWAP(self, label1, label2):
        self.__test_error_SWAP(label1, label2)
        label1 = self.__qnum - 1 - label1
        label2 = self.__qnum - 1 - label2

        matrix1 = np.array([[1, 0], [0, 0]])
        matrix2 = np.array([[0, 0], [0, 1]])

        sig_p = np.array([[0, 1], [0, 0]])
        sig_m = np.array([[0, 0], [1, 0]])

        U = np.zeros((self.__dim, self.__dim)).astype(complex)

        u1, u2, u3, u4 = 1, 1, 1, 1

        for j in range(self.__qnum):
            if (j == label1) or (j == label2):
                u1 = np.kron(u1, matrix1)
            else:
                u1 = np.kron(u1, np.eye(2))
        U += u1

        for j in range(self.__qnum):
            if (j == label1) or (j == label2):
                u2 = np.kron(u2, matrix2)
            else:
                u2 = np.kron(u2, np.eye(2))
        U += u2

        for j in range(self.__qnum):
            if j == label1:
                u3 = np.kron(u3, sig_p)
            elif j == label2:
                u3 = np.kron(u3, sig_m)
            else:
                u3 = np.kron(u3, np.eye(2))
        U += u3

        for j in range(self.__qnum):
            if j == label1:
                u4 = np.kron(u4, sig_m)
            elif j == label2:
                u4 = np.kron(u4, sig_p)
            else:
                u4 = np.kron(u4, np.eye(2))
        U += u4

        self.__state_vec = np.dot(U, self.__state_vec)

    # SWAP gates swapping several qubits.
    def SWAP_multi(self, label_list1, label_list2):
        self.__SWAP_multi_test(label_list1, label_list2)
        for j, label in enumerate(label_list1):
            self.SWAP(label, label_list2[j])

    '''***********************************************************************
    * 3-Qubit Unitary Gates
    ***********************************************************************'''

    # General definition of a controlled-controlled operation.
    def __controlled_controlled_operation(self, u, ctr1_label, ctr2_label, tgt_label):
        self.__test_error_label_CCU(ctr1_label, ctr2_label, tgt_label)

        ctr1_label = self.__qnum - 1 - ctr1_label
        ctr2_label = self.__qnum - 1 - ctr2_label
        tgt_label = self.__qnum - 1 - tgt_label

        matrix1 = np.array([[1, 0], [0, 0]])
        matrix2 = np.array([[0, 0], [0, 1]])

        U = np.zeros((self.__dim, self.__dim)).astype(complex)

        u1, u2, u3 = 1, 1, 1

        for i in range(self.__qnum):
            if i == ctr1_label:
                u1 = np.kron(u1, matrix1)
            else:
                u1 = np.kron(u1, np.eye(2))
        U += u1

        for j in range(self.__qnum):
            if j == ctr1_label:
                u2 = np.kron(u2, matrix2)
            elif j == ctr2_label:
                u2 = np.kron(u2, matrix1)
            else:
                u2 = np.kron(u2, np.eye(2))
        U += u2

        for j in range(self.__qnum):
            if (j == ctr1_label) or (j == ctr2_label):
                u3 = np.kron(u3, matrix2)
            elif j == tgt_label:
                u3 = np.kron(u3, u)
            else:
                u3 = np.kron(u3, np.eye(2))
        U += u3

        self.__state_vec = np.dot(U, self.__state_vec)
        return self.__state_vec

    # CCNOT/Toffoli/Controlled-controlled Pauli X Gate.
    def CCNOT(self, ctr1_label, ctr2_label, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__controlled_controlled_operation(x, ctr1_label, ctr2_label, tgt_label)

    def Toffoli(self, ctr1_label, ctr2_label, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__controlled_controlled_operation(x, ctr1_label, ctr2_label, tgt_label)

    def CCX(self, ctr1_label, ctr2_label, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__controlled_controlled_operation(x, ctr1_label, ctr2_label, tgt_label)

    # CCNOT/Toffoli/Controlled-controlled Pauli X Gates on several qubits.
    def CCNOT_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCNOT(ctr1_label, ctr2_label, tgt_label)

    def Toffoli_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCNOT(ctr1_label, ctr2_label, tgt_label)

    def CCX_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCX(ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled Pauli Y Gate.
    def CCY(self, ctr1_label, ctr2_label, tgt_label):
        y = np.array([[0, -1j], [1j, 0]])
        self.__controlled_controlled_operation(y, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled Pauli Y Gates on several qubits.
    def CCY_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCY(ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled Pauli Z Gate.
    def CCZ(self, ctr1_label, ctr2_label, tgt_label):
        z = np.array([[1, 0], [0, -1]])
        self.__controlled_controlled_operation(z, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled Pauli Z Gates on several qubits.
    def CCZ_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCZ(ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled Hadamard Gate.
    def CCH(self, ctr1_label, ctr2_label, tgt_label):
        h = np.array([[1 / 2 ** 0.5, 1 / 2 ** 0.5], [1 / 2 ** 0.5, -1 / 2 ** 0.5]])
        self.__controlled_controlled_operation(h, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled Hadamard Gates on several qubits.
    def CCH_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCH(ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled S Gate.
    def CCS(self, ctr1_label, ctr2_label, tgt_label):
        s = np.array([[1, 0], [0, 1j]])
        self.__controlled_controlled_operation(s, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled S Gates on several qubits.
    def CCS_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCS(ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled T Gate.
    def CCT(self, ctr1_label, ctr2_label, tgt_label):
        t = np.array([[1, 0], [0, cmath.exp(1j * pi / 4)]])
        self.__controlled_controlled_operation(t, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled T Gates on several qubits.
    def CCT_multi(self, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCT(ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled rotation X gate.
    def CCRx(self, ctr1_label, ctr2_label, tgt_label, theta):
        rx = np.array([[cos(theta / 2), -1j * sin(theta / 2)], [1j * sin(theta / 2), cos(theta) / 2]])
        self.__controlled_controlled_operation(rx, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled rotation X gates on several qubits.
    def CCRx_multi(self, ctr1_label, ctr2_label, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCRx(ctr1_label, ctr2_label, tgt_label, theta)

    # Controlled-controlled rotation Y gate.
    def CCRy(self, ctr1_label, ctr2_label, tgt_label, theta):
        ry = np.array([[cos(theta / 2), -sin(theta / 2)], [sin(theta / 2), cos(theta / 2)]])
        self.__controlled_controlled_operation(ry, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled rotation Y gates on several qubits.
    def CCRy_multi(self, ctr1_label, ctr2_label, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCRy(ctr1_label, ctr2_label, tgt_label, theta)

    # Controlled-controlled rotation Z gate.
    def CCRz(self, ctr1_label, ctr2_label, tgt_label, theta):
        rz = np.array([[1, 0], [0, cmath.exp(1j * theta)]])
        self.__controlled_controlled_operation(rz, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled rotation Z gates on several qubits.
    def CCRz_multi(self, ctr1_label, ctr2_label, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCRz(ctr1_label, ctr2_label, tgt_label, theta)

    # Controlled SWAP/Fradkin gate.
    def CSWAP(self, ctr_label, tgt1_label, tgt2_label):
        self.__test_error_label_CSWAP(ctr_label, tgt1_label, tgt2_label)

        ctr_label = self.__qnum - 1 - ctr_label
        tgt1_label = self.__qnum - 1 -tgt1_label
        tgt2_label = self.__qnum - 1 - tgt2_label

        matrix1 = np.array([[1, 0], [0, 0]])
        matrix2 = np.array([[0, 0], [0, 1]])
        sig_p = np.array([[0, 1], [0, 0]])
        sig_m = np.array([[0, 0], [1, 0]])

        U = np.zeros((self.__dim, self.__dim)).astype(complex)

        u1, u2, u3, u4, u5 = 1, 1, 1, 1, 1

        for j in range(self.__qnum):
            if j == ctr_label:
                u1 = np.kron(u1, matrix1)
            else:
                u1 = np.kron(u1, np.eye(2))
        U += u1

        for j in range(self.__qnum):
            if j == ctr_label:
                u2 = np.kron(u2, matrix2)
            elif (j == tgt1_label) or (j == tgt2_label):
                u2 = np.kron(u2, matrix1)
            else:
                u2 = np.kron(u2, np.eye(2))
        U += u2

        for j in range(self.__qnum):
            if (j == ctr_label) or (j == tgt1_label) or (j == tgt2_label):
                u3 = np.kron(u3, matrix2)
            else:
                u3 = np.kron(u3, np.eye(2))
        U += u3

        for j in range(self.__qnum):
            if j == ctr_label:
                u4 = np.kron(u4, matrix2)
            elif j == tgt1_label:
                u4 = np.kron(u4, sig_p)
            elif j == tgt2_label:
                u4 = np.kron(u4, sig_m)
            else:
                u4 = np.kron(u4, np.eye(2))
        U += u4

        for j in range(self.__qnum):
            if j == ctr_label:
                u5 = np.kron(u5, matrix2)
            elif j == tgt1_label:
                u5 = np.kron(u5, sig_m)
            elif j == tgt2_label:
                u5 = np.kron(u5, sig_p)
            else:
                u5 = np.kron(u5, np.eye(2))
        U += u5

        self.__state_vec = np.dot(U, self.__state_vec)

    def Fradkin(self, ctr_label, tgt1_label, tgt2_label):
        self.__test_error_label_CSWAP(ctr_label, tgt1_label, tgt2_label)

        ctr_label = self.__qnum - 1 - ctr_label
        tgt1_label = self.__qnum - 1 -tgt1_label
        tgt2_label = self.__qnum - 1 - tgt2_label

        matrix1 = np.array([[1, 0], [0, 0]])
        matrix2 = np.array([[0, 0], [0, 1]])
        sig_p = np.array([[0, 1], [0, 0]])
        sig_m = np.array([[0, 0], [1, 0]])

        U = np.zeros((self.__dim, self.__dim)).astype(complex)

        u1, u2, u3, u4, u5 = 1, 1, 1, 1, 1

        for j in range(self.__qnum):
            if j == ctr_label:
                u1 = np.kron(u1, matrix1)
            else:
                u1 = np.kron(u1, np.eye(2))
        U += u1

        for j in range(self.__qnum):
            if j == ctr_label:
                u2 = np.kron(u2, matrix2)
            elif (j == tgt1_label) or (j == tgt2_label):
                u2 = np.kron(u2, matrix1)
            else:
                u2 = np.kron(u2, np.eye(2))
        U += u2

        for j in range(self.__qnum):
            if (j == ctr_label) or (j == tgt1_label) or (j == tgt2_label):
                u3 = np.kron(u3, matrix2)
            else:
                u3 = np.kron(u3, np.eye(2))
        U += u3

        for j in range(self.__qnum):
            if j == ctr_label:
                u4 = np.kron(u4, matrix2)
            elif j == tgt1_label:
                u4 = np.kron(u4, sig_p)
            elif j == tgt2_label:
                u4 = np.kron(u4, sig_m)
            else:
                u4 = np.kron(u4, np.eye(2))
        U += u4

        for j in range(self.__qnum):
            if j == ctr_label:
                u5 = np.kron(u5, matrix2)
            elif j == tgt1_label:
                u5 = np.kron(u5, sig_m)
            elif j == tgt2_label:
                u5 = np.kron(u5, sig_p)
            else:
                u5 = np.kron(u5, np.eye(2))
        U += u5

        self.__state_vec = np.dot(U, self.__state_vec)

    # Controlled SWAP/Fradkin gates on several qubits.
    def CSWAP_multi(self, ctr_label, tgt_label_list1, tgt_label_list2):
        self.__SWAP_multi_test(tgt_label_list1, tgt_label_list2)
        for j, label in enumerate(tgt_label_list1):
            self.CSWAP(ctr_label, label, tgt_label_list2[j])

    def Fradkin_multi(self, ctr_label, tgt_label_list1, tgt_label_list2):
        self.__SWAP_multi_test(tgt_label_list1, tgt_label_list2)
        for j, label in enumerate(tgt_label_list1):
            self.CSWAP(ctr_label, label, tgt_label_list2[j])

    # Controlled-controlled U gate (DIY).
    def CCU_2by2(self, U, ctr1_label, ctr2_label, tgt_label):
        m, n = len(U), len(U[0])
        if ((m != 2) or (n != 2)):
            self.__report_error_input_U(2)
        self.__check_unitary(U)
        self.__controlled_controlled_operation(U, ctr1_label, ctr2_label, tgt_label)

    # Controlled-controlled U gates (DIY) on several qubits.
    def CCU_2by2_multi(self, U, ctr1_label, ctr2_label, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CCU_2by2(U, ctr1_label, ctr2_label, tgt_label)

    '''***********************************************************************
    * Multi-Controlled Unitary Gates
    ***********************************************************************'''

    # General definition of a multi-controlled operation.
    def __multi_controlled_operation(self, u, ctr_label_list, tgt_label):
        self.__test_error_label_CnU(ctr_label_list, tgt_label)
        n = len(ctr_label_list)

        if n == 1:
            self.__controlled_operation(u, ctr_label_list[0], tgt_label)

        else:
            tgt_label = self.__qnum - 1 - tgt_label

            for k in range(n):
                ctr_label_list[k] = self.__qnum - 1 - ctr_label_list[k]

            A = np.array([[1, 0], [0, 0]])
            B = np.array([[0, 0], [0, 1]])

            U = np.zeros((self.__dim, self.__dim)).astype(complex)

            term1, term_last = 1, 1

            for i in range(self.__qnum):
                if i == ctr_label_list[0]:
                    term1 = np.kron(term1, A)
                else:
                    term1 = np.kron(term1, np.eye(2))
            U += term1

            for j in range(self.__qnum):
                if j in ctr_label_list:
                    term_last = np.kron(term_last, B)
                elif j == tgt_label:
                    term_last = np.kron(term_last, u)
                else:
                    term_last = np.kron(term_last, np.eye(2))
            U += term_last

            # remaining (n - 1) terms
            for j in range(n - 1):
                term = 1
                for k in range(self.__qnum):
                    if k in ctr_label_list[: j + 1]:
                        term = np.kron(term, B)
                    elif k == ctr_label_list[j + 1]:
                        term = np.kron(term, A)
                    else:
                        term = np.kron(term, np.eye(2))
                U += term

            self.__state_vec = np.dot(U, self.__state_vec)

        return self.__state_vec

    # Multi-controlled NOT/Pauli X gate.
    def CnNOT(self, ctr_label_list, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__multi_controlled_operation(x, ctr_label_list, tgt_label)

    def CnX(self, ctr_label_list, tgt_label):
        x = np.array([[0, 1], [1, 0]])
        self.__multi_controlled_operation(x, ctr_label_list, tgt_label)

    # Multi-controlled NOT/Pauli X gates on several qubit.
    def CnNOT_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnNOT(ctr_label_list, tgt_label)

    def CnX_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnX(ctr_label_list, tgt_label)

    # Multi-controlled Pauli Y gate.
    def CnY(self, ctr_label_list, tgt_label):
        y = np.array([[0, -1j], [1j, 0]])
        self.__multi_controlled_operation(y, ctr_label_list, tgt_label)

    # Multi-controlled Pauli Y gates on several qubits.
    def CnY_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnY(ctr_label_list, tgt_label)

    # Multi-controlled Pauli Z gate.
    def CnZ(self, ctr_label_list, tgt_label):
        z = np.array([[1, 0], [0, -1]])
        self.__multi_controlled_operation(z, ctr_label_list, tgt_label)

    # Multi-controlled Pauli Y gates on several qubits.
    def CnZ_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnZ(ctr_label_list, tgt_label)

    # Multi-controlled Hadamard gate.
    def CnH(self, ctr_label_list, tgt_label):
        h = np.array([[1 / 2 ** 0.5, 1 / 2 ** 0.5], [1 / 2 ** 0.5, -1 / 2 ** 0.5]])
        self.__multi_controlled_operation(h, ctr_label_list, tgt_label)

    # Multi-controlled Hadamard gates on several qubits.
    def CnH_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnH(ctr_label_list, tgt_label)

    # Multi-controlled S gate.
    def CnS(self, ctr_label_list, tgt_label):
        s = np.array([[1, 0], [0, 1j]])
        self.__multi_controlled_operation(s, ctr_label_list, tgt_label)

    # Multi-controlled S gates on several qubits.
    def CnS_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnS(ctr_label_list, tgt_label)

    # Multi-controlled T gate.
    def CnT(self, ctr_label_list, tgt_label):
        t = np.array([[1, 0], [0, cmath.exp(1j * pi / 4)]])
        self.__multi_controlled_operation(t, ctr_label_list, tgt_label)

    # Multi-controlled T gates on several qubits.
    def CnT_multi(self, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnT(ctr_label_list, tgt_label)

    # Multi-controlled rotation X gate.
    def CnRx(self, ctr_label_list, tgt_label, theta):
        rx = np.array([[cos(theta / 2), -1j * sin(theta / 2)], [1j * sin(theta / 2), cos(theta) / 2]])
        self.__multi_controlled_operation(rx, ctr_label_list, tgt_label)

    # Multi-controlled rotation X gates on several qubits.
    def CnRx_multi(self, ctr_label_list, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnRx(ctr_label_list, tgt_label, theta)

    # Multi-controlled rotation Y gate.
    def CnRy(self, ctr_label_list, tgt_label, theta):
        ry = np.array([[cos(theta / 2), -sin(theta / 2)], [sin(theta / 2), cos(theta / 2)]])
        self.__multi_controlled_operation(ry, ctr_label_list, tgt_label)

    # Multi-controlled rotation Y gates on several qubits.
    def CnRy_multi(self, ctr_label_list, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnRy(ctr_label_list, tgt_label, theta)

    # Multi-controlled rotation Z gate.
    def CnRz(self, ctr_label_list, tgt_label, theta):
        rz = np.array([[1, 0], [0, cmath.exp(1j * theta)]])
        self.__multi_controlled_operation(rz, ctr_label_list, tgt_label)

    # Multi-controlled rotation Z gates on several qubits.
    def CnRz_multi(self, ctr_label_list, tgt_label_list, theta):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnRz(ctr_label_list, tgt_label, theta)

    # Multi-controlled U gate (DIY).
    def CnU_2by2(self, U, ctr_label_list, tgt_label):
        m, n = len(U), len(U[0])
        if ((m != 2) or (n != 2)):
            self.__report_error_input_U(2)
        self.__check_unitary(U)
        self.__multi_controlled_operation(U, ctr_label_list, tgt_label)

    # Multi-controlled U gates (DIY) on several qubits.
    def CnU_2by2_multi(self, U, ctr_label_list, tgt_label_list):
        self.__test_error_multi_label_list(tgt_label_list)
        for tgt_label in tgt_label_list:
            self.CnU_2by2(U, ctr_label_list, tgt_label)

    # Multi-controlled SWAP gate.
    def CnSWAP(self, ctr_label_list, tgt1_label, tgt2_label):
        self.__test_error_label_CnSWAP(ctr_label_list, tgt1_label, tgt2_label)
        n = len(ctr_label_list)

        tgt1_label = self.__qnum - 1 - tgt1_label
        tgt2_label = self.__qnum - 1 -tgt2_label
        for k in range(n):
            ctr_label_list[k] = self.__qnum - 1 - ctr_label_list[k]

        A = np.array([[1, 0], [0, 0]])
        B = np.array([[0, 0], [0, 1]])

        sig_1 = np.array([[0, 1], [0, 0]])
        sig_2 = np.array([[0, 0], [1, 0]])

        U = np.zeros((self.__dim, self.__dim)).astype(complex)
        pink1, pink2, pink3, pink4, pink5 = 1, 1, 1, 1, 1

        # pink 1
        for j in range(self.__qnum):
            if j == ctr_label_list[0]:
                pink1 = np.kron(pink1, A)
            else:
                pink1 = np.kron(pink1, np.eye(2))
        U += pink1

        # pink 2
        for j in range(self.__qnum):
            if j in ctr_label_list:
                pink2 = np.kron(pink2, B)
            elif (j == tgt1_label) or (j == tgt2_label):
                pink2 = np.kron(pink2, A)
            else:
                pink2 = np.kron(pink2, np.eye(2))
        U += pink2

        # pink 3
        for j in range(self.__qnum):
            if j in ctr_label_list:
                pink3 = np.kron(pink3, B)
            elif j == tgt1_label:
                pink3 = np.kron(pink3, sig_1)
            elif j == tgt2_label:
                pink3 = np.kron(pink3, sig_2)
            else:
                pink3 = np.kron(pink3, np.eye(2))
        U += pink3

        # pink 4
        for j in range(self.__qnum):
            if j in ctr_label_list:
                pink4 = np.kron(pink4, B)
            elif j == tgt1_label:
                pink4 = np.kron(pink4, sig_2)
            elif j == tgt2_label:
                pink4 = np.kron(pink4, sig_1)
            else:
                pink4 = np.kron(pink4, np.eye(2))
        U += pink4

        # pink 5
        for j in range(self.__qnum):
            if (j in ctr_label_list) or (j == tgt1_label) or (j == tgt2_label):
                pink5 = np.kron(pink5, B)
            else:
                pink5 = np.kron(pink5, np.eye(2))

        # greens
        for j in range(n - 1):
            term = 1
            for k in range(self.__qnum):
                if k in ctr_label_list[: j + 1]:
                    term = np.kron(term, B)
                elif k == ctr_label_list[j + 1]:
                    term = np.kron(term, A)
                else:
                    term = np.kron(term, np.eye(2))
            U += term

        self.__state_vec = np.dot(U, self.__state_vec)

        return self.__state_vec

    # Multi-controlled SWAP gates on several qubits.
    def CnSWAP_multi(self, ctr_label_list, tgt_label_list1, tgt_label_list2):
        self.__SWAP_multi_test(tgt_label_list1, tgt_label_list2)
        for j, label in enumerate(tgt_label_list1):
            self.CnSWAP(ctr_label_list, label, tgt_label_list2[j])

    '''***********************************************************************
    * Error Reporting
    ***********************************************************************'''

    # Check whether the state vector is normalized.
    def __state_normalization_check(self, state_vec):
        sum = 0
        for j in range(len(state_vec)):
            sum += state_vec[j, 0] * state_vec[j, 0].conjugate()
        print("sum = {}".format(sum))
        if np.allclose(1.0, sum.real) == False:
            raise Exception("invalid state vector.")

    # Check whether n is an integer.
    def __positive_integer_check(self, n):
        if (ceil(n) != floor(n)) or (n <= 0):
            raise Exception("the number of qubit(s) must be a positive integer.")

    # Check whether n is a natural number.
    def __natural_number_check(self, n):
        if (ceil(n) != floor(n)) or (n < 0):
            raise Exception("invalid label(s) of qubit(s), which must be natural number(s).")

    # When the input of U is non-unitary, throw the error.
    def __report_error_input_U(self, type):
        if type == 0:
            raise Exception("invalid input of U, which must be unitary.")
        else:
            raise Exception("invalid input of U, which must be a {} × {} unitary matrix.".format(type, type))

    # Error testing in the definition of SWAP gate.
    def __test_error_SWAP(self, label1, label2):
        self.__test_error_label(label1)
        self.__test_error_label(label2)
        if label1 == label2:
            raise Exception("two qubits to be swapped must be different from each other.")

    # Check whether the input label is valid.
    def __test_error_label(self, label):
        self.__natural_number_check(label)
        if label >= self.__qnum:
            raise Exception("invalid label(s) of qubit(s), which cannot exceed the maximum label {}.".format(self.__qnum - 1))

    # Error testing in the definitions of multi-SWAP operation.
    def __SWAP_multi_test(self, label_list1, label_list2):
        n1 = len(label_list1)
        n2 = len(label_list2)
        if n1 != n2:
            raise Exception("two lists of qubits to be swapped must have the same length.")
        for j, label in enumerate(label_list1):
            if label == label_list2[j]:
                raise Exception("qubits to be swapped cannot be the same one.")
            for k in range(j + 1, n1):
                if label == label_list1[k]:
                    raise Exception("qubits in the 1st list must be different from each other.")
        for j, label in enumerate(label_list2):
            for k in range(j + 1, n1):
                if label == label_list2[k]:
                    raise Exception("qubits in the 2nd list must be different from each other.")

    # Error testing in the definitions of controlled operations.
    def __test_error_controlled(self, ctr_label):
        self.__natural_number_check(ctr_label)
        if ctr_label >= self.__qnum:
            raise Exception("invalid label(s) of controlled qubit(s), which cannot exceed the maximum label {}.".format(self.__qnum - 1))

    def __test_error_target(self, tgt_label):
        self.__natural_number_check(tgt_label)
        if tgt_label >= self.__qnum:
            raise Exception("invalid label of target qubit, which cannot exceed the maximum label {}.".format(self.__qnum - 1))

    def __test_error_label_CU(self, ctr_label, tgt_label):
        self.__test_error_controlled(ctr_label)
        self.__test_error_target(tgt_label)
        if ctr_label == tgt_label:
            raise Exception("the controlled qubit(s) and target qubit cannot be the same one.")

    def __test_error_label_CCU(self, ctr1_label, ctr2_label, tgt_label):
        self.__test_error_controlled(ctr1_label)
        self.__test_error_controlled(ctr2_label)
        self.__test_error_target(tgt_label)
        if (ctr1_label == tgt_label) or (ctr2_label == tgt_label):
            raise Exception("the controlled qubit(s) and target qubit cannot be the same one.")
        if ctr1_label == ctr2_label:
            raise Exception("the controlled qubits must be the different from each other")

    def __test_error_label_CnU(self, ctr_label_list, tgt_label):
        for ctr_label in ctr_label_list:
            self.__test_error_controlled(ctr_label)
        self.__test_error_target(tgt_label)
        for j, ctr_label in enumerate(ctr_label_list):
            if ctr_label_list[j] == tgt_label:
                raise Exception("the controlled qubit(s) and target qubit cannot be the same one.")
            if ctr_label in ctr_label_list[j + 1:]:
                raise Exception("the controlled qubits must be the different from each other.")

    def __test_error_label_CSWAP(self, ctr_label, tgt1_label, tgt2_label):
        self.__test_error_controlled(ctr_label)
        self.__test_error_target(tgt1_label)
        self.__test_error_target(tgt2_label)
        self.__test_error_SWAP(tgt1_label, tgt2_label)
        if (ctr_label == tgt1_label) or (ctr_label == tgt2_label):
            raise Exception("the controlled qubit(s) and target qubit(s) cannot be the same one.")

    def __test_error_label_CnSWAP(self, ctr_label_list, tgt1_label, tgt2_label):
        for j, ctr_label in enumerate(ctr_label_list):
            self.__test_error_controlled(ctr_label)
            if (ctr_label == tgt1_label) or (ctr_label == tgt2_label):
                raise Exception("the controlled qubit(s) and target qubit(s) cannot be the same one.")
            if ctr_label in ctr_label_list[j + 1:]:
                raise Exception("the controlled qubits must be the different from each other.")
        self.__test_error_target(tgt1_label)
        self.__test_error_target(tgt2_label)
        self.__test_error_SWAP(tgt1_label, tgt2_label)

    # Check whether a set of labels are valid.
    def __test_error_multi_label_list(self, tgt_label_list):
        for j, tgt_label in enumerate(tgt_label_list):
            for k in range(j + 1, len(tgt_label_list)):
                if tgt_label == tgt_label_list[k]:
                    raise Exception("labels in list of target qubits must be different from each other.")

    # The restriction for depicting.
    def __depict_error(self):
        if self.__qnum >= 5:
            raise Exception("number of qubits cannot greater than 4 in Piscis when depicting.\n You can use 'print_state_distribution' instead.")

    # Error when there is no measurement before output.
    def __output_error(self, length):
        if length <= 2:
            raise Exception("you should set at least 1 measurement on the circuit before outputing the result.")

    # When initialize the state by a set of amplitudes, make sure that the length of input is equal to that of the state vector.
    def __set_amp_error(self, length):
        if length != self.__dim:
            raise Exception("list of amplitudes must have the length {}.".format(self.__dim))
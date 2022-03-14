## APIs, Piscis V0.0.0.1 

### One-qubit unitary gates (methods in QCirc):

- Pauli X gate: `X(label)`
- Pauli Y gate: `Y(label)`
- Pauli Z gate: `Z(label)`
- Hadmard gate: `H(label)`
- S gate: `S(label)`
- T gate: `T(label)`
- Rotation X gate: `Rx(label, theta)`
- Rotation Y gate: `Ry(label, theta)`
- Rotation Z gate: `Rz(label, theta)`
- Define a 2 × 2 unitary gate: `U_2by2(np_array_2by2_unitary, label)`

Each of the single qubit gates has two variants, `multi` and `all` such as `X_multi(label_list)` and `X_all()`.

### Two-qubit unitary gates (methods in QCirc):
(`ctr` means `controlled` and `tgt` means `target`)
- CNOT / Controlled Pauli X gate : `CNOT(ctr_label, tgt_label)` / `CX(ctr_label, tgt_label)`
- Controlled Pauli Y gate: `CY(ctr_label, tgt_label)`
- Controlled Pauli Z gate: `CZ(ctr_label, tgt_label)`
- Controlled Hadamrd gate: `CH(ctr_label, tgt_label)`
- Controlled S gate: `CS(ctr_label, tgt_label)`
- Controlled T gate: `CT(ctr_label, tgt_label)`
- Controlled rotation X gate: `CRx(ctr_label, tgt_label, theta)`
- Controlled rotation Y gate: `CRy(ctr_label, tgt_label, theta)`
- Controlled rotation Z gate: `CRz(ctr_label, tgt_label, theta)`
- Define a controlled-unitary (2 × 2) gate: `CU_2by2(np_array_2by2_unitary, ctr_label, tgt_label)`
- SWAP gate: `SWAP(qubit1_label, qubit2_label)`

The controlled operations and SWAP operation have a variant `multi` such as `CNOT_multi(ctr_label, tgt_label_list)`.

### Three-qubit unitary gates (methods in QCirc):
- Toffoli / CCNOT / controlled-controlled Pauli X gate: `Toffoli(ctr1_label, ctr2_label, tgt_label)` / `CCNOT(ctr1_label, ctr2_label, tgt_label)` /`CCX(ctr1_label, ctr2_label, tgt_label)`
- Controlled-controlled Pauli Y gate: `CCY(ctr1_label, ctr2_label, tgt_label)`
- Controlled-controlled Pauli Z gate: `CCZ(ctr1_label, ctr2_label, tgt_label)`
- Controlled-controlled Hadamard gate: `CCH(ctr1_label, ctr2_label, tgt_label)`
- Controlled-controlled S gate: `CCS(ctr1_label, ctr2_label, tgt_label)`
- Controlled-controlled T gate: `CCT(ctr1_label, ctr2_label, tgt_label)`
- Controlled-controlled rotation X gate: `CCRx(ctr1_label, ctr2_label, tgt_label, theta)`
- Controlled-controlled rotation Y gate: `CCRy(ctr1_label, ctr2_label, tgt_label, theta)`
- Controlled-controlled rotation Z gate: `CCRz(ctr1_label, ctr2_label, tgt_label, theta)`
- Define a controlled-controlled unitary (2 × 2) gate: `CCU_2by2(np_array_2by2_unitary, ctr1_label, ctr2_label, tgt_label)`
- Fradkin / CSWAP gate: `Fradkin(ctr_label, tgt1_label, tgt2_label)` / `CSWAP(ctr_label, tgt1_label, tgt2_label)` 

The controlled-controlled operations and Fradkin / CSWAP operation have a variant `multi` such as `CSWAP_multi(ctr_label, tgt_label_list1, tgt2_label_list2)`.

### Multi-controlled unitary gates (methods in QCirc):
- Multi-controlled NOT / Pauli X gate: `CnNOT(ctr_label_list, tgt_label)` / `CnX(ctr_list, tgt_label)`
- Multi-controlled Pauli Y gate: `CnY(ctr_label_list, tgt_label)`
- Multi-controlled Pauli Z gate: `CnZ(ctr_label_list, tgt_label)`
- Multi-controlled Hadamard gate: `CnH(ctr_label_list, tgt_label)`
- Multi-controlled S gate: `CnS(ctr_label_list, tgt_label)`
- Multi-controlled T gate: `CnT(ctr_label_list, tgt_label)`
- Multi-controlled rotation X gate: `CnRx(ctr_label_list, tgt_label, theta)`
- Multi-controlled rotation Y gate: `CnRy(ctr_label_list, tgt_label, theta)`
- Multi-controlled rotation Z gate: `CnRz(ctr_label_list, tgt_label, theta)`
- Define a multi-controlled unitary (2 × 2) gate: `CnU_2by2(np_array_2by2_unitary, ctr_list, tgt_label)`
- Multi controlled SWAP gate: `CnSWAP(U,ctr_qubit_list,tgt1_qubit_label, tgt2_qubit_label)`

The multi-controlled operations have a variant `multi` such as `CnY_multi(ctr_label_list, tgt_label_list)`.

### Non-unitary operations (methods in QCirc):
- Measure a qubit: `measure(label)`
- Measure some qubits: `measure_multi(label_list)`
- Measure all qubits: `measure_all()`
- Project a specified qubit to zero: `projector_to_zero(label)`
- Project a specified qubit to one: `projector_to_one(label)`


### Statistical tools
- Output (read) the circuit after meausrements: `output()`
- Print the state vector: `print_state()`
- Print the state distribution: `print_state_distribution()`
- Depict the state distribuition (`qnum`< 5): `depict_state_distribution()`
- Return the probability (float-type) of a specified qubit to be zero: `get_prob_of_zero(label)`
- Return the probability (float-type) of a specified qubit to be one: `get_prob_of_one(label)`
- Return the amplitude (complex-type) of a specified qubit to be zero: `get_amp_of_zero(label)`
- Return the amplitude (complex-type) of a specified qubit to be one: `get_amp_of_one(label)`

### State initialization
- Initialize a state by a set of amplitudes: `ini_set_amp(amp_list)`
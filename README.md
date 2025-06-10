# Tranzor

> WARNING: This project is a preliminary mock-up and has not been extensively tested,
so it may encounter bugs.

Tranzor is a generator for transversal Clifford logical circuits using the color code.

Given a *logical-level* circuit represented in `stim`, `tranzor` compiles each qubit in the circuit
into a logical qubit encoded with the superdense color code, and each gate into a transversal gate
supported by the color code. A syndrome extraction round is explicitly represented using the `I` instruction
in the circuit.

## Installation

```sh
pip install tranzor
```

## Examples

### S Gate

```python
import stim
from tranzor import SuperDenseColorCode, compile_circuit

logical_circuit = stim.Circuit("""
QUBIT_COORDS(0, 0) 0
RX 0
TICK
I 0
TICK
S_DAG 0
TICK
I 0
TICK
S 0
TICK
I 0
TICK
MX 0
OBSERVABLE_INCLUDE(0) rec[-1]
""")

color_code = SuperDenseColorCode(d=3)
out = compile_circuit(logical_circuit, color_code)

```

compiles to the following [Crumble circuit](https://algassert.com/crumble#circuit=Q(0,0)0;Q(1,0)1;Q(1,1)2;Q(1,2)3;Q(2,0)4;Q(2,1)5;Q(2,2)6;Q(2,3)7;Q(3,0)8;Q(3,1)9;Q(3,2)10;Q(4,0)11;Q(4,1)12;R_4_6_12;RX_1_3_9_0_2_5_7_8_10_11;TICK;CX_1_4_9_12_3_6;TICK;CX_2_1_10_9_7_6_5_4;TICK;CX_0_1_5_9_10_6_8_4;TICK;CX_8_9_2_3_5_6_11_12;TICK;CX_1_2_9_10_6_7_4_5;TICK;CX_1_0_9_5_6_10_4_8;TICK;CX_9_8_3_2_6_5_12_11;TICK;CX_1_4_9_12_3_6;TICK;M_4_6_12;MX_1_3_9;DT(1,0,0)rec[-3];DT(1,2,0)rec[-2];DT(3,1,0)rec[-1];TICK;R_4_6_12;RX_1_3_9;TICK;CX_1_4_9_12_3_6;TICK;CX_2_1_10_9_7_6_5_4;TICK;CX_0_1_5_9_10_6_8_4;TICK;CX_8_9_2_3_5_6_11_12;TICK;CX_1_2_9_10_6_7_4_5;TICK;CX_1_0_9_5_6_10_4_8;TICK;CX_9_8_3_2_6_5_12_11;TICK;CX_1_4_9_12_3_6;TICK;M_4_6_12;MX_1_3_9;DT(2,0,1)rec[-6]_rec[-11]_rec[-12];DT(2,2,1)rec[-5]_rec[-12];DT(4,1,1)rec[-4];DT(1,0,1)rec[-3]_rec[-9];DT(1,2,1)rec[-2]_rec[-8];DT(3,1,1)rec[-1]_rec[-7];TICK;S_2_8_10;S_DAG_0_5_7_11;TICK;R_4_6_12;RX_1_3_9;TICK;CX_1_4_9_12_3_6;TICK;CX_2_1_10_9_7_6_5_4;TICK;CX_0_1_5_9_10_6_8_4;TICK;CX_8_9_2_3_5_6_11_12;TICK;CX_1_2_9_10_6_7_4_5;TICK;CX_1_0_9_5_6_10_4_8;TICK;CX_9_8_3_2_6_5_12_11;TICK;CX_1_4_9_12_3_6;TICK;M_4_6_12;MX_1_3_9;DT(2,0,2)rec[-6]_rec[-11]_rec[-12];DT(2,2,2)rec[-5]_rec[-12];DT(4,1,2)rec[-4];DT(1,0,2)rec[-3]_rec[-9]_rec[-11]_rec[-12];DT(1,2,2)rec[-2]_rec[-8]_rec[-12];DT(3,1,2)rec[-1]_rec[-7];TICK;S_0_5_7_11;S_DAG_2_8_10;TICK;R_4_6_12;RX_1_3_9;TICK;CX_1_4_9_12_3_6;TICK;CX_2_1_10_9_7_6_5_4;TICK;CX_0_1_5_9_10_6_8_4;TICK;CX_8_9_2_3_5_6_11_12;TICK;CX_1_2_9_10_6_7_4_5;TICK;CX_1_0_9_5_6_10_4_8;TICK;CX_9_8_3_2_6_5_12_11;TICK;CX_1_4_9_12_3_6;TICK;M_4_6_12;MX_1_3_9;DT(2,0,3)rec[-6]_rec[-11]_rec[-12];DT(2,2,3)rec[-5]_rec[-12];DT(4,1,3)rec[-4];DT(1,0,3)rec[-3]_rec[-9]_rec[-11]_rec[-12];DT(1,2,3)rec[-2]_rec[-8]_rec[-12];DT(3,1,3)rec[-1]_rec[-7];TICK;MX_0_2_5_7_8_10_11;DT(1,0,4)rec[-3]_rec[-5]_rec[-6]_rec[-7]_rec[-10];DT(1,2,4)rec[-2]_rec[-4]_rec[-5]_rec[-6]_rec[-9];DT(3,1,4)rec[-1]_rec[-2]_rec[-3]_rec[-5]_rec[-8];OI(0)rec[-1]_rec[-2]_rec[-3]_rec[-4]_rec[-5]_rec[-6]_rec[-7]_rec[-17]_rec[-18]).

### H and CNOT Gate

```python
import stim
from tranzor import SuperDenseColorCode, compile_circuit

logical_circuit = stim.Circuit("""
QUBIT_COORDS(0, 0) 0
QUBIT_COORDS(1, 0) 1
RZ 0 1
TICK
H 0
TICK
I 0 1
TICK
CX 0 1
TICK
I 0 1
TICK
MZ 0 1
OBSERVABLE_INCLUDE(0) rec[-1] rec[-2]
""")

color_code = SuperDenseColorCode(d=3)
out = compile_circuit(logical_circuit, color_code)
```

compiles to the following [Crumble circuit](https://algassert.com/crumble#circuit=Q(0,0)0;Q(1,0)1;Q(1,1)2;Q(1,2)3;Q(2,0)4;Q(2,1)5;Q(2,2)6;Q(2,3)7;Q(3,0)8;Q(3,1)9;Q(3,2)10;Q(4,0)11;Q(4,1)12;Q(5,0)13;Q(6,0)14;Q(6,1)15;Q(6,2)16;Q(7,0)17;Q(7,1)18;Q(7,2)19;Q(7,3)20;Q(8,0)21;Q(8,1)22;Q(8,2)23;Q(9,0)24;Q(9,1)25;R_4_6_12_0_2_5_7_8_10_11_17_19_25_13_15_18_20_21_23_24;RX_1_3_9_14_16_22;TICK;CX_1_4_9_12_3_6_14_17_22_25_16_19;TICK;CX_2_1_10_9_7_6_5_4_15_14_23_22_20_19_18_17;TICK;CX_0_1_5_9_10_6_8_4_13_14_18_22_23_19_21_17;TICK;CX_8_9_2_3_5_6_11_12_21_22_15_16_18_19_24_25;TICK;CX_1_2_9_10_6_7_4_5_14_15_22_23_19_20_17_18;TICK;CX_1_0_9_5_6_10_4_8_14_13_22_18_19_23_17_21;TICK;CX_9_8_3_2_6_5_12_11_22_21_16_15_19_18_25_24;TICK;CX_1_4_9_12_3_6_14_17_22_25_16_19;TICK;M_4_6_12_17_19_25;MX_1_3_9_14_16_22;DT(2,0,0)rec[-12];DT(2,2,0)rec[-11];DT(4,1,0)rec[-10];DT(7,0,0)rec[-9];DT(7,2,0)rec[-8];DT(9,1,0)rec[-7];TICK;H_0_8_2_5_11_10_7;TICK;R_4_6_12_17_19_25;RX_1_3_9_14_16_22;TICK;CX_1_4_9_12_3_6_14_17_22_25_16_19;TICK;CX_2_1_10_9_7_6_5_4_15_14_23_22_20_19_18_17;TICK;CX_0_1_5_9_10_6_8_4_13_14_18_22_23_19_21_17;TICK;CX_8_9_2_3_5_6_11_12_21_22_15_16_18_19_24_25;TICK;CX_1_2_9_10_6_7_4_5_14_15_22_23_19_20_17_18;TICK;CX_1_0_9_5_6_10_4_8_14_13_22_18_19_23_17_21;TICK;CX_9_8_3_2_6_5_12_11_22_21_16_15_19_18_25_24;TICK;CX_1_4_9_12_3_6_14_17_22_25_16_19;TICK;M_4_6_12_17_19_25;MX_1_3_9_14_16_22;DT(1,0,1)rec[-12]_rec[-18];DT(1,2,1)rec[-11]_rec[-17];DT(3,1,1)rec[-10]_rec[-16];DT(7,0,1)rec[-9]_rec[-20]_rec[-21];DT(7,2,1)rec[-8]_rec[-21];DT(9,1,1)rec[-7];DT(2,0,1)rec[-6]_rec[-23]_rec[-24];DT(1.5,1,1)rec[-5]_rec[-24];DT(3,1,2)rec[-4];DT(6,0,1)rec[-3]_rec[-15];DT(6,2,1)rec[-2]_rec[-14];DT(8,1,1)rec[-1]_rec[-13];TICK;CX[transversal]_0_13_8_21_2_15_5_18_11_24_10_23_7_20;TICK;R_4_6_12_17_19_25;RX_1_3_9_14_16_22;TICK;CX_1_4_9_12_3_6_14_17_22_25_16_19;TICK;CX_2_1_10_9_7_6_5_4_15_14_23_22_20_19_18_17;TICK;CX_0_1_5_9_10_6_8_4_13_14_18_22_23_19_21_17;TICK;CX_8_9_2_3_5_6_11_12_21_22_15_16_18_19_24_25;TICK;CX_1_2_9_10_6_7_4_5_14_15_22_23_19_20_17_18;TICK;CX_1_0_9_5_6_10_4_8_14_13_22_18_19_23_17_21;TICK;CX_9_8_3_2_6_5_12_11_22_21_16_15_19_18_25_24;TICK;CX_1_4_9_12_3_6_14_17_22_25_16_19;TICK;M_4_6_12_17_19_25;MX_1_3_9_14_16_22;DT(2,0,3)rec[-12]_rec[-23]_rec[-24];DT(2,2,3)rec[-11]_rec[-24];DT(4,1,3)rec[-10];DT(7,0,3)rec[-9]_rec[-20]_rec[-21]_rec[-23]_rec[-24];DT(7,2,3)rec[-8]_rec[-21]_rec[-24];DT(9,1,3)rec[-7];DT(1,0,3)rec[-6]_rec[-15]_rec[-18];DT(1,2,3)rec[-5]_rec[-14]_rec[-17];DT(3,1,3)rec[-4]_rec[-13]_rec[-16];DT(6,0,3)rec[-3]_rec[-15];DT(6,2,3)rec[-2]_rec[-14];DT(8,1,3)rec[-1]_rec[-13];TICK;M_0_2_5_7_8_10_11_13_15_18_20_21_23_24;DT(2,0,4)rec[-10]_rec[-12]_rec[-13]_rec[-14]_rec[-25]_rec[-26];DT(2,3,4)rec[-9]_rec[-11]_rec[-12]_rec[-13]_rec[-26];DT(2,1,4)rec[-8]_rec[-9]_rec[-10]_rec[-12];DT(7,0,4)rec[-3]_rec[-5]_rec[-6]_rec[-7]_rec[-22]_rec[-23];DT(7,3,4)rec[-2]_rec[-4]_rec[-5]_rec[-6]_rec[-23];DT(7,1,4)rec[-1]_rec[-2]_rec[-3]_rec[-5];OI(0)rec[-1]_rec[-2]_rec[-3]_rec[-4]_rec[-5]_rec[-6]_rec[-7]_rec[-8]_rec[-9]_rec[-10]_rec[-11]_rec[-12]_rec[-13]_rec[-14]_rec[-21]_rec[-22]_rec[-24]_rec[-25]_rec[-33]_rec[-34]_rec[-45]_rec[-46]).

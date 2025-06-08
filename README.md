# Tranzor

> WARNING: This project is a preliminary mock-up and has not been extensively tested,
so it may encounter bugs.

Tranzor is a generator for transversal Clifford logical circuits using the color code.

Given a *logical-level* circuit represented in `stim`, `tranzor` compiles each qubit in the circuit
into a logical qubit encoded with the superdense color code, and each gate into a transversal gate
supported by the color code. A syndrome extraction round is explicitly represented using the `I` instruction
in the circuit.

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


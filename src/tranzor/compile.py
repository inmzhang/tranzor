import functools
from itertools import zip_longest
from typing import Iterable

import stim

from tranzor.utils import (
    Moment,
    iter_stim_circuit_without_repeat_by_moments,
    QubitMap,
    MeasurementRecordsMap,
    has_measurement,
    has_reset,
    has_I,
    is_measurement,
    is_reset,
    has_only_measurement_or_is_virtual,
    has_only_reset_or_is_virtual,
    LogicalCoordsMap,
)
from tranzor.base_code import BaseCode


COMMON_SUPPORTED_GATES = [
    "RZ",
    "RX",
    "R",
    "MZ",
    "MX",
    "M",
    "I",
    "OBSERVABLE_INCLUDE",
    "TICK",
    "QUBIT_COORDS",
    "CNOT",
    "CX",
]


def compile_circuit(circuit: stim.Circuit, code: BaseCode) -> stim.Circuit:
    """Compiles a logical circuit into a physical circuit encoded by the given code.

    The logical gates in the logical circuit should be transversal and supported
    by the underlying code.

    Args:
        circuit: The logical circuit to compile.
        code: The quantum error correction code used for encoding the logical qubit
            and implementing the transversal logical gates.

    Returns:
        A stim.Circuit object representing the physical circuit that implements
        the logical circuit on the code.
    """
    circuit = circuit.flattened().without_noise()

    coords_map = _resolve_logical_qubit_coords_map(circuit, code)
    circuit = _filter_qubit_coords_instructions(circuit)
    _validate_logical_circuit(circuit, code, coords_map.logical_indices())
    obs_map = _resolve_observable_include(circuit)

    moments = _split_circuit_and_add_observable_tag(circuit)
    compiled_circuit = coords_map.build_coords_circuit()
    qubit_map = QubitMap.from_circuit(compiled_circuit)
    measurement_record = MeasurementRecordsMap()
    for moment_idx, moment in enumerate(moments):
        moment_inner: list[list[Moment]] = []
        # record detection region sources that back-propagate in the circuit
        detectors_source: list[stim.PauliString] = []
        for logical_inst in moment.circuit:
            assert isinstance(logical_inst, stim.CircuitInstruction)
            if logical_inst.name == "OBSERVABLE_INCLUDE":
                continue
            for target_group in logical_inst.target_groups():
                if len(target_group) == 2:
                    t1, t2 = target_group[0].qubit_value, target_group[1].qubit_value
                    assert t1 is not None and t2 is not None
                    moment_inner.append(
                        list(
                            code.two_qubit_gate(
                                logical_inst.name,
                                coords_map.c2i(t1),
                                coords_map.c2i(t2),
                            )
                        )
                    )
                    continue
                logical_qubit = target_group[0].qubit_value
                assert logical_qubit is not None
                if is_reset(logical_inst):
                    basis = "X" if logical_inst.name == "RX" else "Z"
                    src = stim.PauliString(circuit.num_qubits)
                    src[logical_qubit] = basis
                    detectors_source.append(src)
                    obs_basis = _split_obs_tag(logical_inst.tag)
                    assert all(basis == b for _, b in obs_basis.items())
                    moment_inner.append(
                        list(
                            code.reset(basis, coords_map.c2i(logical_qubit), obs_basis)
                        )
                    )
                elif is_measurement(logical_inst):
                    basis = "X" if logical_inst.name == "MX" else "Z"
                    src = stim.PauliString(circuit.num_qubits)
                    src[logical_qubit] = basis
                    detectors_source.append(src)
                    moment_inner.append(
                        list(
                            code.measure(
                                basis,
                                coords_map.c2i(logical_qubit),
                                include_observable=obs_map.get(logical_qubit, ()),
                            )
                        )
                    )
                elif logical_inst.name == "I":
                    src_x = stim.PauliString(circuit.num_qubits)
                    src_x[logical_qubit] = "X"
                    detectors_source.append(src_x)

                    src_z = stim.PauliString(circuit.num_qubits)
                    src_z[logical_qubit] = "Z"
                    detectors_source.append(src_z)

                    obs_basis = _split_obs_tag(logical_inst.tag)
                    moment_inner.append(
                        list(
                            code.syndrome_extraction(
                                coords_map.c2i(logical_qubit), obs_basis
                            )
                        )
                    )
                else:
                    moment_inner.append(
                        list(
                            code.single_qubit_gate(
                                logical_inst.name, coords_map.c2i(logical_qubit)
                            )
                        )
                    )
        circuit_for_current_moment = stim.Circuit()
        for slice in zip_longest(*moment_inner, fillvalue=Moment(stim.Circuit())):
            merged_moment = functools.reduce(
                lambda a, b: Moment(a.circuit + b.circuit),
                slice,
                Moment(stim.Circuit()),
            )
            if merged_moment.circuit:
                circuit_for_current_moment.append(merged_moment.circuit)
                circuit_for_current_moment.append("TICK", [], [])
        measurement_record_for_current_moment = MeasurementRecordsMap.from_circuit(
            circuit_for_current_moment, qubit_map
        )
        if detectors_source:
            circuit_for_current_moment = circuit_for_current_moment[
                :-1
            ]  # Remove the last TICK
            for source_pauli in detectors_source:
                source_qubit = next(
                    q for q in range(len(source_pauli)) if source_pauli[q] != 0
                )
                source_basis = source_pauli[source_qubit]
                correlated_detectors = (
                    None
                    if has_reset(moment.circuit)
                    else _resolve_detector_correlation(
                        source_pauli, moments[:moment_idx]
                    )
                )
                detector_circuit = code.detectors(
                    current_logical_qubit=source_qubit,
                    basis="X" if source_basis == 1 else "Z",
                    coords_map=coords_map,
                    measurement_record_before_current_moment=measurement_record,
                    measurement_record_for_current_moment=measurement_record_for_current_moment,
                    is_data_qubit_readout=has_measurement(moment.circuit),
                    correlated_detectors=correlated_detectors,
                )
                circuit_for_current_moment += detector_circuit
            circuit_for_current_moment.append("SHIFT_COORDS", [], [0, 0, 1])
            circuit_for_current_moment.append("TICK", [], [])
        compiled_circuit += circuit_for_current_moment
        measurement_record = measurement_record.with_added_measurements(
            measurement_record_for_current_moment
        )

    return compiled_circuit


def _validate_logical_circuit(
    circuit: stim.Circuit, code: BaseCode, all_qubits: set[int]
) -> None:
    """Validates that the logical circuit is compatible with the given code.

    Validation checks include:
        - Instructions between two ticks do not operate on overlapping qubits.
        - Every qubits starts with a reset instruction in the first moment.
        - Instructions should be Reset/Measure/I or transversal logical gates
          supported by the code.
        - A moment that includes R/M instructions does not include any other
          type of instruction.
        - A moment includes I instruction does not include any other type of
          instruction. And the I instruction should operate on all the qubits.
    """
    # Moment data structure has ensured that the instructions will not clash
    moments = iter_stim_circuit_without_repeat_by_moments(circuit)
    # The first moment should be a reset moment.
    first_moment = next(moments)
    if not has_only_reset_or_is_virtual(first_moment.circuit):
        raise ValueError(
            f"The first moment should be a reset moment, but got:\n {first_moment}"
        )
    # Check instructions are supported
    for instruction in circuit:
        if instruction.name in COMMON_SUPPORTED_GATES:
            continue
        supported_transversal_gates = code.supported_transversal_gates()
        if instruction.name not in supported_transversal_gates:
            raise ValueError(
                f"Unsupported logical gate {instruction.name} in the circuit.",
                " Supported gates:",
                f"{COMMON_SUPPORTED_GATES + sorted(supported_transversal_gates)}",
            )
    for moment in moments:
        # Check exclusive reset
        if has_reset(moment.circuit) and not has_only_reset_or_is_virtual(
            moment.circuit
        ):
            raise ValueError(
                f"Moment {moment} includes reset instructions, but also "
                "includes other type of instructions, which is not allowed."
            )
        # Check exclusive measure
        if has_measurement(moment.circuit) and not has_only_measurement_or_is_virtual(
            moment.circuit
        ):
            raise ValueError(
                f"Moment {moment} includes measurement instructions, but also "
                "includes other type of instructions, which is not allowed."
            )
        # Check exclusive I (syndrome extraction)
        if has_I(moment.circuit):
            if any(inst.name != "I" for inst in moment.circuit):
                raise ValueError(
                    f"Moment {moment} includes I instructions, but also "
                    "includes other type of instructions, which is not allowed."
                )
            targets: set[int] = set(
                target.qubit_value
                for inst in moment.circuit
                for target in inst.targets_copy()  # type: ignore
            )
            if targets != all_qubits:
                raise ValueError(
                    f"Moment '{moment.circuit}' includes I instructions, but does not "
                    "operate on all the qubits. I instruction represents syndrome extraction, "
                    "which should operate on all the logical qubits at the same time."
                )


def _resolve_observable_include(circuit: stim.Circuit) -> dict[int, list[int]]:
    """Returns a mapping from the measurement index in the circuit to the list of
    logical observable indices that include the measurement."""
    num_measurements = 0
    obs_map: dict[int, list[int]] = {}
    for inst in circuit:
        assert isinstance(inst, stim.CircuitInstruction)
        if is_measurement(inst):
            num_measurements += len(inst.targets_copy())
        elif inst.name == "OBSERVABLE_INCLUDE":
            obs_idx = int(inst.gate_args_copy()[0])
            targets = [t.value + num_measurements for t in inst.targets_copy()]
            assert all(t >= 0 for t in targets)
            for t in targets:
                obs_map.setdefault(t, []).append(obs_idx)
    return obs_map


def _resolve_logical_qubit_coords_map(
    circuit: stim.Circuit, code: BaseCode
) -> LogicalCoordsMap:
    shift_unit = code.bounding_box

    local_index_map: dict[int, dict[complex, int]] = {}
    global_coords_map: dict[int, complex] = {}

    logical_qubit_coords = circuit.get_final_qubit_coordinates()
    if (
        not logical_qubit_coords
        or len(logical_qubit_coords) != circuit.num_qubits
        or any(len(coord) != 2 for coord in logical_qubit_coords.values())
    ):
        raise ValueError(
            "The circuit should have specified coordinates for all qubits, "
            "and the coordinates should be 2D coordinates (x, y)."
        )

    global_idx = 0
    for logical_idx, logical_coords in logical_qubit_coords.items():
        for local_coords in _sort_complex(code.used_qubits):
            local_index_map.setdefault(logical_idx, {})[local_coords] = global_idx
            global_coords_map[global_idx] = complex(
                logical_coords[0] * shift_unit.real + local_coords.real,
                logical_coords[1] * shift_unit.imag + local_coords.imag,
            )
            global_idx += 1

    return LogicalCoordsMap(local_index_map, global_coords_map)


def _filter_qubit_coords_instructions(circuit: stim.Circuit) -> stim.Circuit:
    """Filters out QUBIT_COORDS instructions from the circuit."""
    filtered_circuit = stim.Circuit()
    for inst in circuit:
        assert isinstance(inst, stim.CircuitInstruction)
        if inst.name != "QUBIT_COORDS":
            filtered_circuit.append(inst)
    return filtered_circuit


def _sort_complex(values: Iterable[complex], reverse=False) -> list[complex]:
    return sorted(values, key=lambda z: (z.real, z.imag), reverse=reverse)


def _resolve_detector_correlation(
    source_pauli: stim.PauliString,
    previous_moments: list[Moment],
) -> stim.PauliString:
    propagated_pauli = source_pauli.copy()
    correlated_detectors: stim.PauliString = stim.PauliString(len(source_pauli))
    for moment in previous_moments[::-1]:
        if has_reset(moment.circuit):
            for q in moment.qubits_indices:
                if propagated_pauli[q] != 0:
                    correlated_detectors[q] = propagated_pauli[q]
                    propagated_pauli[q] = 0
            if propagated_pauli.weight == 0:
                break
            continue
        if has_I(moment.circuit):
            for q in range(len(propagated_pauli)):
                pauli = propagated_pauli[q]
                if pauli != 0:
                    correlated_detectors[q] = pauli
                    propagated_pauli[q] = 0
            break
        propagated_pauli = propagated_pauli.after(moment.circuit.inverse())
    assert propagated_pauli.weight == 0
    return correlated_detectors


def _split_circuit_and_add_observable_tag(
    logical_circuit: stim.Circuit,
) -> list[Moment]:
    """Superdense color code cycle will include additional measurements in the observable.

    This function add observable basis tags to the I instructions.
    """
    detecting_regions = logical_circuit.detecting_regions(targets="L")
    tagged_circuit: list[Moment] = []

    for tick, moment in enumerate(
        iter_stim_circuit_without_repeat_by_moments(logical_circuit)
    ):
        if has_reset(moment.circuit) or has_I(moment.circuit):
            moment_circuit = stim.Circuit()
            for inst in moment.circuit:
                assert isinstance(inst, stim.CircuitInstruction)
                for target in inst.targets_copy():
                    qubit = target.qubit_value
                    assert qubit is not None
                    tags: list[str] = []
                    for obs, regions in detecting_regions.items():
                        if tick not in regions:
                            continue
                        pauli = regions[tick][qubit]
                        if pauli != 0:
                            tags.append(f"{obs} {'IXYZ'[pauli]}")
                    splitted_inst = stim.CircuitInstruction(
                        inst.name, [target], inst.gate_args_copy(), tag=",".join(tags)
                    )
                    moment_circuit.append(splitted_inst)
            tagged_circuit.append(Moment(moment_circuit))
        else:
            tagged_circuit.append(moment)
    return tagged_circuit


def _split_obs_tag(tag: str) -> dict[int, str]:
    res: dict[int, str] = {}
    if not tag:
        return res
    for part in tag.split(","):
        obs, pauli = part.split(" ")
        res[int(obs[1:])] = pauli
    return res

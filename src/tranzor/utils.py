from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
import functools
from typing import Callable, Iterable, Iterator, Sequence, cast

import numpy as np
import stim


@dataclass(frozen=True)
class LogicalCoordsMap:
    local_index_map: dict[int, dict[complex, int]]
    """A mapping from the index of the logical qubit to a map that maps the local
    qubit coordinates in that logical qubit to a global qubit index."""
    global_coords_map: dict[int, complex]
    """A mapping from the global qubit index to the global coordinates of the qubit."""

    def logical_indices(self) -> set[int]:
        """Returns the set of logical qubit indices."""
        return set(self.local_index_map.keys())

    def physical_indices(self) -> set[int]:
        """Returns the set of global qubit indices."""
        return set(self.global_coords_map.keys())

    def c2i(self, logical_qubit: int) -> dict[complex, int]:
        return self.local_index_map[logical_qubit]

    def build_coords_circuit(self) -> stim.Circuit:
        """Builds a stim.Circuit that contains QUBIT_COORDS instructions for all
        global qubit indices."""
        coords_circuit = stim.Circuit()
        for global_idx, global_coords in self.global_coords_map.items():
            coords_circuit.append(
                "QUBIT_COORDS",
                [global_idx],
                [global_coords.real, global_coords.imag],
            )
        return coords_circuit


@dataclass(frozen=True)
class QubitMap:
    """Represent a bijection between qubits and their associated indices."""

    i2q: dict[int, complex] = field(default_factory=dict)

    def __post_init__(self) -> None:
        qubit_counter = Counter(self.i2q.values())
        if len(qubit_counter) != len(self.i2q):
            duplicated_qubits = frozenset(
                q for q in qubit_counter if qubit_counter[q] > 1
            )
            raise ValueError(
                f"Found qubit(s) with more than one index: {duplicated_qubits}."
            )

    @staticmethod
    def from_circuit(circuit: stim.Circuit) -> QubitMap:
        qubit_coordinates = circuit.get_final_qubit_coordinates()
        qubits: dict[int, complex] = {}
        for qi, coords in qubit_coordinates.items():
            if len(coords) != 2:
                raise ValueError(
                    "Qubits should be defined on exactly 2 spatial dimensions. "
                    f"Found {qi} -> {coords} defined on {len(coords)} spatial dimensions."
                )
            qubits[qi] = complex(coords[0], coords[1])
        return QubitMap(qubits)

    @functools.cached_property
    def q2i(self) -> dict[complex, int]:
        return {q: i for i, q in self.i2q.items()}

    def to_circuit(self) -> stim.Circuit:
        """Get a circuit with only ``QUBIT_COORDS`` instructions representing
        ``self``."""
        ret = stim.Circuit()
        for qi, qubit in sorted(self.i2q.items(), key=lambda t: t[0]):
            ret.append("QUBIT_COORDS", qi, (qubit.real, qubit.imag))
        return ret

    def __getitem__(self, qubit: complex) -> int:
        return self.q2i[qubit]


@dataclass(frozen=True)
class MeasurementRecordsMap:
    """A mapping from measurements appearing in a circuit to their record offsets.

    This class stores record offsets which are, by essence, relative to a certain
    position in a circuit. This means that this class and the measurement offsets
    it stores are meaningless without knowledge about the circuit containing the
    represented measurements and the position(s) in the circuit at which the
    instance at hand is valid.
    """

    mapping: dict[complex, list[int]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        all_measurement_records_indices: list[int] = []
        for qubit, measurement_record_offsets in self.mapping.items():
            # Check that the provided measurement record offsets are negative.
            nonnegative_offsets = [
                offset for offset in measurement_record_offsets if offset >= 0
            ]
            if nonnegative_offsets:
                raise ValueError(
                    "Invalid mapping from qubit offsets to measurement record "
                    f"offsets. Found positive offsets ({nonnegative_offsets}) for "
                    f"qubit {qubit}."
                )
            # Check that measurement record offsets are sorted
            if measurement_record_offsets != sorted(measurement_record_offsets):
                raise ValueError(
                    "Got measurement record offsets that are not in sorted "
                    f"order: {measurement_record_offsets}. This is not supported."
                )
            all_measurement_records_indices.extend(measurement_record_offsets)
        # Check that a given measurement record offset only appears once.
        deduplicated_indices = np.unique(all_measurement_records_indices)
        if len(deduplicated_indices) != len(all_measurement_records_indices):
            raise ValueError(
                "At least one measurement record offset has been found twice "
                "in the provided offsets."
            )

    @staticmethod
    def from_circuit(
        circuit: stim.Circuit, qubit_map: QubitMap | None = None
    ) -> MeasurementRecordsMap:
        current_measurement_record = -circuit.num_measurements
        if qubit_map is None:
            qubit_map = QubitMap.from_circuit(circuit)
        measurement_records: dict[complex, list[int]] = {}
        for instruction in circuit:
            assert isinstance(instruction, stim.CircuitInstruction)
            if stim.gate_data(instruction.name).produces_measurements:
                for (qi,) in instruction.target_groups():
                    qubit = qubit_map.i2q[qi.value]
                    measurement_records.setdefault(qubit, []).append(
                        current_measurement_record
                    )
                    current_measurement_record += 1
        assert current_measurement_record == 0
        return MeasurementRecordsMap(measurement_records)

    def __getitem__(self, qubit: complex) -> Sequence[int]:
        return self.mapping[qubit]

    def __contains__(self, qubit: complex) -> bool:
        return qubit in self.mapping

    def with_added_measurements(
        self, mrecords_map: MeasurementRecordsMap, repetitions: int = 1
    ) -> MeasurementRecordsMap:
        num_measurements_without_repetition = sum(
            len(offsets) for offsets in mrecords_map.mapping.values()
        )
        num_added_measurements = repetitions * num_measurements_without_repetition
        records = {
            q: [o - num_added_measurements for o in offsets]
            for q, offsets in self.mapping.items()
        }
        for q, offsets in mrecords_map.mapping.items():
            record = records.setdefault(q, [])
            for i in range(repetitions - 1, -1, -1):
                record.extend(
                    [o - i * num_measurements_without_repetition for o in offsets]
                )
        return MeasurementRecordsMap(records)


class Moment:
    """A collection of instructions that can be executed in parallel."""

    def __init__(self, circuit: stim.Circuit) -> None:
        """Initialize a :class:`Moment` instance.

        Args:
            circuit: collection of instructions representing the :class:`Moment`.
                It should represent a valid moment, see the class documentation
                for a detailed explanation of the pre-conditions.
        """
        Moment.check_is_valid_moment(circuit)
        self._circuit: stim.Circuit = circuit
        self._used_qubits: set[int]
        self._used_qubits = set(_count_qubit_accesses(circuit).keys())

    @property
    def circuit(self) -> stim.Circuit:
        return self._circuit

    @property
    def qubits_indices(self) -> set[int]:
        return self._used_qubits

    @staticmethod
    def check_is_valid_moment(circuit: stim.Circuit) -> None:
        """Check if the provided circuit can be considered a valid moment."""
        if circuit.num_ticks > 0:
            raise ValueError(
                "Cannot initialize a Moment with a stim.Circuit instance "
                "containing at least one TICK instruction."
            )
        qubit_usage = _count_qubit_accesses(circuit)
        multi_used_qubits = [
            qi for qi, usage_count in qubit_usage.items() if usage_count > 1
        ]
        if multi_used_qubits:
            raise ValueError(
                "Moment instances cannot be initialized with a stim.Circuit "
                "instance containing gates applied on the same qubit. Found "
                "multiple gates applied on the following qubits: "
                f"{multi_used_qubits}."
            )
        if any(isinstance(inst, stim.CircuitRepeatBlock) for inst in circuit):
            raise ValueError(
                "Moment instances should no contain any instance "
                "of stim.CircuitRepeatBlock."
            )

    def append(
        self,
        name_or_instr: str | stim.CircuitInstruction,
        targets: Iterable[int | stim.GateTarget] | None = None,
        args: Iterable[float] | None = None,
    ) -> None:
        """Append an instruction to the :class:`Moment`."""
        if targets is None:
            targets = tuple()
        if args is None:
            args = tuple()

        instruction: stim.CircuitInstruction
        if isinstance(name_or_instr, str):
            instruction = stim.CircuitInstruction(name_or_instr, targets, args)
        else:
            instruction = name_or_instr

        if _is_annotation_instruction(instruction):
            self._circuit.append(instruction)
            return

        # Checking Moment invariant
        instruction_qubits = Moment._get_used_qubit_indices(
            targets if isinstance(name_or_instr, str) else name_or_instr.targets_copy()
        )
        overlapping_qubits = self._used_qubits.intersection(instruction_qubits)
        if overlapping_qubits:
            raise ValueError(
                f"Cannot add {instruction} to the Moment due to qubit(s) "
                f"{overlapping_qubits} being already in use."
            )
        self._used_qubits.update(instruction_qubits)
        self._circuit.append(instruction)

    @staticmethod
    def _get_used_qubit_indices(
        targets: Iterable[int | stim.GateTarget],
    ) -> list[int]:
        qubits: list[int] = []
        for target in targets:
            if isinstance(target, int):
                qubits.append(target)
                continue
            # isinstance(target, stim.GateTarget)
            if target.is_qubit_target:
                assert isinstance(target.qubit_value, int)  # type checker is happy
                qubits.append(target.qubit_value)
        return qubits


NOISE_INSTRUCTION_NAMES: frozenset[str] = frozenset(
    [
        "CORRELATED_ERROR",
        "DEPOLARIZE1",
        "DEPOLARIZE2",
        "E",
        "ELSE_CORRELATED_ERROR",
        "HERALDED_ERASE",
        "HERALDED_PAULI_CHANNEL_1",
        "PAULI_CHANNEL_1",
        "PAULI_CHANNEL_2",
        "X_ERROR",
        "Y_ERROR",
        "Z_ERROR",
    ]
)
ANNOTATION_INSTRUCTION_NAMES: frozenset[str] = frozenset(
    ["DETECTOR", "MPAD", "OBSERVABLE_INCLUDE", "QUBIT_COORDS", "SHIFT_COORDS", "TICK"]
)


def _count_qubit_accesses(circuit: stim.Circuit) -> dict[int, int]:
    """Count the number of times a given qubit is used by an instruction that
    is not an annotation."""
    counter: defaultdict[int, int] = defaultdict(int)
    for instruction in circuit:
        if isinstance(instruction, stim.CircuitRepeatBlock):
            for qi, count in _count_qubit_accesses(instruction.body_copy()).items():
                counter[qi] += count * instruction.repeat_count
        else:
            if (
                instruction.name
                in NOISE_INSTRUCTION_NAMES | ANNOTATION_INSTRUCTION_NAMES
            ):
                continue
            for target in instruction.targets_copy():
                # Ignore targets that are not qubit targets.
                if not target.is_qubit_target:
                    continue
                qi = cast(int, target.qubit_value)
                counter[qi] += 1
    return counter


def _is_annotation_instruction(instruction: stim.CircuitInstruction) -> bool:
    return instruction.name in ANNOTATION_INSTRUCTION_NAMES


def iter_stim_circuit_without_repeat_by_moments(
    circuit: stim.Circuit, collected_before_use: bool = True
) -> Iterator[Moment]:
    """Iterate over the ``stim.Circuit`` by moments."""
    copy_func: Callable[[stim.Circuit], stim.Circuit] = (
        (lambda c: c.copy()) if collected_before_use else (lambda c: c)
    )
    cur_moment = stim.Circuit()
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            raise ValueError(
                "Found an instance of stim.CircuitRepeatBlock which is "
                "explicitly not supported by this method."
            )
        elif inst.name == "TICK":
            yield Moment(copy_func(cur_moment))
            cur_moment.clear()
        else:
            cur_moment.append(inst)
    # No need to copy the last moment
    yield Moment(cur_moment)


def is_measurement(instruction: stim.CircuitInstruction) -> bool:
    return instruction.name in ["M", "MX", "MZ"]


def is_reset(instruction: stim.CircuitInstruction) -> bool:
    return stim.gate_data(instruction.name).is_reset  # type: ignore


def is_annotation(instruction: stim.CircuitInstruction) -> bool:
    return instruction.name in ANNOTATION_INSTRUCTION_NAMES


def is_noisy_gate(instruction: stim.CircuitInstruction) -> bool:
    return instruction.name in ANNOTATION_INSTRUCTION_NAMES


def has_I(circuit: stim.Circuit) -> bool:
    """Returns True if the circuit contains I instructions."""
    return any(inst.name == "I" for inst in circuit)


def has_measurement(moment: stim.Circuit) -> bool:
    return any(is_measurement(inst) for inst in moment)  # type:ignore


def is_virtual_instruction(inst: stim.CircuitInstruction) -> bool:
    return is_annotation(inst) or is_noisy_gate(inst)


def has_only_measurement_or_is_virtual(moment: stim.Circuit) -> bool:
    for inst in moment:
        if is_virtual_instruction(inst):  # type: ignore
            continue
        if not is_measurement(inst):  # type: ignore
            return False
    return True


def has_reset(moment: stim.Circuit) -> bool:
    return any(is_reset(inst) for inst in moment)  # type:ignore


def has_only_reset_or_is_virtual(moment: stim.Circuit) -> bool:
    for inst in moment:
        if is_virtual_instruction(inst):  # type: ignore
            continue
        if not is_reset(inst):  # type: ignore
            return False
    return True

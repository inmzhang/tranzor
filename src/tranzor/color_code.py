from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Literal, override

import stim

from tranzor.base_code import BaseCode
from tranzor.utils import MeasurementRecordsMap, Moment

if TYPE_CHECKING:
    from tranzor.compile import LogicalCoordsMap


class SuperDenseColorCode(BaseCode):
    def __init__(self, d: int):
        assert d % 2 == 1 and d >= 3, (
            "d must be an odd integer greater than or equal to 3"
        )
        super().__init__(d)

        base_width = self.base_width
        potential_x_mqs = [
            x + 1j * y
            for x in range(-1, base_width, 2)
            for y in range((x // 2) % 2, base_width, 2)
        ]
        self._data_qubits: set[complex] = set()
        self._x_mqs: set[complex] = set()
        self._z_mqs: set[complex] = set()

        dq_shifts: list[complex] = [-1, +1j, +1j + 1, +2, -1j + 1, -1j]
        for mq in potential_x_mqs:
            dq_in_bounds = [
                mq + shift for shift in dq_shifts if self._is_in_bounds(mq + shift)
            ]
            if len(dq_in_bounds) < 4:
                continue
            self._data_qubits.update(dq_in_bounds)
            self._x_mqs.add(mq)
            self._z_mqs.add(mq + 1)

    @classmethod
    @override
    def supported_transversal_gates(cls) -> frozenset[str]:
        return frozenset(["S", "S_DAG", "H", "CX", "CNOT"])

    @property
    @override
    def data_qubits(self) -> frozenset[complex]:
        return frozenset(self._data_qubits)

    @property
    def measure_qubits(self) -> frozenset[complex]:
        return frozenset(self._x_mqs | self._z_mqs)

    @property
    @override
    def used_qubits(self) -> frozenset[complex]:
        return frozenset(self.data_qubits | self.measure_qubits)

    @override
    def single_qubit_gate(
        self, gate: str, qubit_map: dict[complex, int]
    ) -> Iterable[Moment]:
        moment = stim.Circuit()
        if gate == "H":
            moment.append("H", [qubit_map[dq] for dq in self.data_qubits], [])
        elif gate in ["S", "S_DAG"]:
            s_targets: set[complex] = set()
            s_dag_targets: set[complex] = set()
            for mq in self._x_mqs:
                s_targets.update(
                    mq + shift
                    for shift in [-1, 1 + 1j, 1 - 1j]
                    if (mq + shift) in self.data_qubits
                )
                s_dag_targets.update(
                    mq + shift
                    for shift in [-1j, 1j, 2]
                    if (mq + shift) in self.data_qubits
                )
            if gate == "S_DAG":
                s_targets, s_dag_targets = s_dag_targets, s_targets
            moment.append("S", sorted(qubit_map[q] for q in s_targets), [])
            moment.append("S_DAG", sorted(qubit_map[q] for q in s_dag_targets), [])
        else:
            raise ValueError(f"Unsupported gate: {gate}")
        return [Moment(moment)]

    @property
    @override
    def bounding_box(self) -> complex:
        height = max(self._data_qubits, key=lambda q: q.imag).imag + 1
        return complex(self.base_width, height)

    def mq_to_dqs(self, mq: complex) -> list[complex]:
        assert mq in self.measure_qubits

        if mq in self._x_mqs:
            shift = [-1, 1j, -1j, 1 + 1j, 1 - 1j, 2]
        else:
            shift = [-2, 1j, -1j, -1 + 1j, -1 - 1j, 1]
        return [mq + s for s in shift if (mq + s) in self.data_qubits]

    @property
    def base_width(self) -> int:
        return 2 * self.d - 1

    def _is_in_bounds(self, q: complex) -> bool:
        base_width = self.base_width
        if q.real < 0 or q.imag < 0 or q.real >= base_width or q.imag >= base_width:
            return False
        if q.imag * 2 > q.real * 3:
            return False
        if q.imag * 2 > (base_width - q.real) * 3:
            return False
        return True

    @override
    def reset(
        self,
        basis: str,
        qubit_map: dict[complex, int],
        observable_basis: dict[int, str] | None = None,
    ) -> list[Moment]:
        moments = self.syndrome_extraction(qubit_map, observable_basis=observable_basis)
        moments[0].append(f"R{basis}", sorted(qubit_map[q] for q in self._data_qubits))
        return moments

    @override
    def measure(
        self,
        basis: str,
        qubit_map: dict[complex, int],
        include_observable: Iterable[int] = (),
    ) -> list[Moment]:
        circuit = stim.Circuit()
        circuit.append(
            f"M{basis}",
            sorted(qubit_map[q] for q in self._data_qubits),
            [],
        )
        for obs_idx in include_observable:
            circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(-i - 1) for i in range(len(self._data_qubits))],
                [obs_idx],
            )
        return [Moment(circuit)]

    @override
    def syndrome_extraction(
        self,
        qubit_map: dict[complex, int],
        observable_basis: dict[int, str] | None = None,
    ) -> list[Moment]:
        used_qubits = self.used_qubits
        moment_circuits: list[stim.Circuit] = [stim.Circuit() for _ in range(10)]

        def do_cxs(
            out: stim.Circuit,
            centers: Iterable[complex],
            d_control: complex,
            d_target: complex,
        ) -> None:
            targets: list[int] = []
            for c in centers:
                if c + d_control in used_qubits and c + d_target in used_qubits:
                    targets.append(qubit_map[c + d_control])
                    targets.append(qubit_map[c + d_target])

            out.append("CX", targets, [])

        moment_circuits[0].append("RX", sorted(qubit_map[q] for q in self._x_mqs), [])
        moment_circuits[0].append("RZ", sorted(qubit_map[q] for q in self._z_mqs), [])

        do_cxs(moment_circuits[1], self._x_mqs, 0, 1)
        do_cxs(moment_circuits[2], self._x_mqs, 1j, 0)
        do_cxs(moment_circuits[2], self._z_mqs, 1j, 0)
        do_cxs(moment_circuits[3], self._x_mqs, -1, 0)
        do_cxs(moment_circuits[3], self._z_mqs, 1, 0)
        do_cxs(moment_circuits[4], self._x_mqs, -1j, 0)
        do_cxs(moment_circuits[4], self._z_mqs, -1j, 0)
        do_cxs(moment_circuits[5], self._x_mqs, 0, 1j)
        do_cxs(moment_circuits[5], self._z_mqs, 0, 1j)
        do_cxs(moment_circuits[6], self._x_mqs, 0, -1)
        do_cxs(moment_circuits[6], self._z_mqs, 0, 1)
        do_cxs(moment_circuits[7], self._x_mqs, 0, -1j)
        do_cxs(moment_circuits[7], self._z_mqs, 0, -1j)
        do_cxs(moment_circuits[8], self._x_mqs, 0, 1)

        x_mq_indices = sorted(qubit_map[q] for q in self._x_mqs)
        z_mq_indices = sorted(qubit_map[q] for q in self._z_mqs)

        moment_circuits[9].append("MX", x_mq_indices, [])
        moment_circuits[9].append("MZ", z_mq_indices, [])
        if observable_basis:
            for obs_idx, basis in observable_basis.items():
                if basis == "X":
                    continue
                moment_circuits[9].append(
                    "OBSERVABLE_INCLUDE",
                    [
                        stim.target_rec(
                            z_mq_indices.index(qubit_map[q]) - len(z_mq_indices)
                        )
                        for q in self._z_mqs
                        if q.imag > 0
                    ],
                    [obs_idx],
                )

        return [Moment(circuit) for circuit in moment_circuits]

    @override
    def detectors(
        self,
        current_logical_qubit: int,
        basis: Literal["X", "Z"],
        coords_map: LogicalCoordsMap,
        measurement_record_before_current_moment: MeasurementRecordsMap,
        measurement_record_for_current_moment: MeasurementRecordsMap,
        is_data_qubit_readout: bool = False,
        correlated_detectors: stim.PauliString | None = None,
    ) -> stim.Circuit:
        local_map = coords_map.c2i(current_logical_qubit)
        mqs = self._x_mqs if basis == "X" else self._z_mqs

        num_measurements_in_current_moment = sum(
            len(offsets)
            for offsets in measurement_record_for_current_moment.mapping.values()
        )
        detector_circuit = stim.Circuit()
        for mq in mqs:
            targets: list[int] = []
            if is_data_qubit_readout:
                dqs = [
                    coords_map.global_coords_map[local_map[q]]
                    for q in self.mq_to_dqs(mq)
                ]
                targets.extend(
                    measurement_record_for_current_moment[dq][-1] for dq in dqs
                )
            else:
                targets.append(
                    measurement_record_for_current_moment[
                        coords_map.global_coords_map[local_map[mq]]
                    ][-1]
                )
            if correlated_detectors is not None:
                for logical_qubit in range(len(correlated_detectors)):
                    pauli = correlated_detectors[logical_qubit]
                    if pauli == 0:
                        continue
                    correlated_mqs: list[complex] = []
                    mqx, mqz = (mq - 1, mq) if basis == "Z" else (mq, mq + 1)
                    if pauli == 1:
                        correlated_mqs = [mqx]
                    else:
                        for q in [mqz - 2j if mqz.imag > 0 else mqz, mqz + 2j]:
                            if q in self.measure_qubits:
                                correlated_mqs.append(q)
                        if pauli == 2:
                            correlated_mqs.append(mqx)
                    for q in correlated_mqs:
                        global_index = coords_map.c2i(logical_qubit)[q]
                        global_coords = coords_map.global_coords_map[global_index]
                        targets.append(
                            measurement_record_before_current_moment[global_coords][-1]
                            - num_measurements_in_current_moment
                        )
            detector_circuit.append(
                "DETECTOR", [stim.target_rec(t) for t in targets], [mq.real, mq.imag, 0]
            )
        return detector_circuit

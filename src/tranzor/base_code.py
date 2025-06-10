from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import TYPE_CHECKING, Iterable, Literal

import stim

from tranzor.utils import Moment, MeasurementRecordsMap

if TYPE_CHECKING:
    from tranzor.utils import LogicalCoordsMap


TRANSVERSAL_TAG = "transversal"


class BaseCode(metaclass=ABCMeta):
    """Base class for concrete implementations of CSS codes that admits certain
    transversal logical gates."""

    def __init__(self, d: int):
        self.d = d

    @property
    @abstractmethod
    def bounding_box(self) -> complex:
        """Returns the bounding box of the code."""
        pass

    @classmethod
    @abstractmethod
    def supported_transversal_gates(cls) -> frozenset[str]:
        """Returns a set of supported transversal logical gates. The names should
        be a valid `stim` instruction name."""
        pass

    @property
    @abstractmethod
    def data_qubits(self) -> frozenset[complex]:
        """Returns the data qubits of the code."""
        pass

    @property
    @abstractmethod
    def used_qubits(self) -> frozenset[complex]:
        """Returns all the used qubits of the code."""
        pass

    @abstractmethod
    def single_qubit_gate(
        self, gate: str, qubit_map: dict[complex, int]
    ) -> Iterable[Moment]:
        """Returns stim circuit moments for a single qubit transversal gate."""
        pass

    def two_qubit_gate(
        self,
        gate: str,
        first_qubit_map: dict[complex, int],
        second_qubit_map: dict[complex, int],
    ) -> Iterable[Moment]:
        """Returns stim circuit moments for a two qubit transversal gate."""
        assert gate in self.supported_transversal_gates()
        moment = stim.Circuit()
        for dq in self.data_qubits:
            first_dq, second_dq = first_qubit_map[dq], second_qubit_map[dq]
            moment.append(gate, [first_dq, second_dq], [], tag=TRANSVERSAL_TAG)
        return [Moment(moment)]

    @abstractmethod
    def reset(
        self,
        basis: str,
        qubit_map: dict[complex, int],
        observable_basis: dict[int, str] | None = None,
    ) -> Iterable[Moment]:
        """Returns stim circuit moments for resetting/initializing the code in
        a certain basis."""
        pass

    @abstractmethod
    def measure(
        self,
        basis: str,
        qubit_map: dict[complex, int],
        include_observable: Iterable[int] = (),
    ) -> Iterable[Moment]:
        """Returns stim circuit moments for measuring the code in a certain basis,
        optionally including a observable annotation."""
        pass

    @abstractmethod
    def syndrome_extraction(
        self,
        qubit_map: dict[complex, int],
        observable_basis: dict[int, str] | None = None,
    ) -> Iterable[Moment]:
        """Returns stim circuit moments for syndrome extraction."""
        pass

    @abstractmethod
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
        pass

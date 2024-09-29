from surf.src.pqubit import pqubit

class surf_step:
    """
    Data structure for storing information related to each time step in a circuit.
    """

    def __init__(self, pqubits: list[pqubit], rounds: int, observables: list[(str,list[pqubit],bool)], basis: str):
        self.pqubits: list[pqubit] = pqubits
        self.rounds: int = rounds
        self.observables: list[(str,list[pqubit],bool)] = observables
        self.basis: str = basis

    def get_pqubits(self) -> list[pqubit]:
        return self.pqubits

    def get_rounds(self) -> int:
        return self.rounds

    def get_observables(self) -> list[list[pqubit]]:
        return self.observables

    def set_pqubits(self, pqubits: list[pqubit]):
        self.pqubits = pqubits

    def set_rounds(self, rounds: int):
        self.rounds = rounds

    def set_observables(self, observables: list[list[pqubit]]):
        self.observables = observables

    def add_observable(self, observable: list[pqubit]):
        self.observables.append(observable)

    def get_basis(self) -> str:
        return self.basis

    def set_basis(self, basis: str):
        self.basis = basis





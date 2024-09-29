import warnings

from surf.src.pqubit import pqubit
from surf.src.qubit_topology import qubit_topology
from math import ceil


class surf_patch(qubit_topology):
    """
    Surface-code patch to simulate in stim.
    """

    def __init__(self, dx: int, dy: int, r: int, basis: str, init_state: str, origin: tuple[float, float],
                 edge_types=None, uniform_error_rate=1e-3, id='', pqubit_id_offset=0):
        """
        :param dx,dy: Code distance in x or y directions of the tile
        :param r: Number of QEC rounds
        :param basis: Basis to initialize the patch in
        :param init_state: State within basis to initialize patch in
        :param origin: Location of bottom left corner of the patch
        :param edge_types: Dict of stabilizer edge types of the patch
        """

        super().__init__(id, pqubit_id_offset, origin)

        self.dx: int = dx
        self.dy: int = dy
        self.rounds: int = r
        self.basis: str = basis
        self.init_state: str = init_state
        self.uniform_error_rate: float = uniform_error_rate

        self.obs_qubits: list[pqubit] = []
        self.z_obs: list[pqubit] = []
        self.x_obs: list[pqubit] = []

        if edge_types is None:
            self.edge_types: dict[str, str] = {'N': 'z', 'S': 'z', 'E': 'x', 'W': 'x'}
        else:
            self.edge_types: dict[str, str] = edge_types

        qubit_id_counter = pqubit_id_offset

        # Init data qubits
        for y in range(0, 2 * dy, 2):
            for x in range(0, 2 * dx, 2):
                pq = pqubit(qubit_id_counter, y + 1 + origin[0], x + 1 + origin[1], 'd', init_state)
                self.qubits.append(pq)
                self.data_qubits.append(pq)
                qubit_id_counter += 1

        # Init weight-4 X parity qubits
        for y in range(2, dy * 2, 2):

            init_x = 4 if y % 4 == 2 else 2

            for x in range(init_x, dx * 2 , 4):
                pq = pqubit(qubit_id_counter, y + origin[0], x + origin[1], 'x', '0')
                self.qubits.append(pq)
                self.parity_qubits.append(pq)

                qubit_id_counter += 1

        # Init weight-4 Z parity qubits
        for y in range(2, (dy * 2), 2):

            init_x = 4 if y % 4 == 0 else 2

            for x in range(init_x, (dx * 2), 4):
                pq = pqubit(qubit_id_counter, y + origin[0], x + origin[1], 'z', '0')
                self.qubits.append(pq)
                self.parity_qubits.append(pq)

                qubit_id_counter += 1

        # Init internal qubit/coordinate maps
        self._reset_qubit_maps()

        # North edge type
        y = (dy*2)
        init_x = 2

        for x in range(init_x, (dx * 2), 2):

            if self.xy_map[(y-2  + origin[0], x + origin[1])].get_type() != self.edge_types['N']:
                pq = pqubit(qubit_id_counter, y + origin[0], x + origin[1], self.edge_types['N'], '0')
                self.qubits.append(pq)
                self.parity_qubits.append(pq)

                qubit_id_counter += 1

        # South edge type
        y = 0
        init_x = 2

        for x in range(init_x, (dx * 2), 2):

            if self.xy_map[(y+2 + origin[0], x + origin[1])].get_type() != self.edge_types['S']:
                pq = pqubit(qubit_id_counter, y + origin[0], x + origin[1], self.edge_types['S'], '0')
                self.qubits.append(pq)
                self.parity_qubits.append(pq)

                qubit_id_counter += 1

        # East edge type
        init_y = 2
        x = (dx*2)

        for y in range(init_y, (dy * 2), 2):

            if self.xy_map[(y + origin[0], x-2 + origin[1])].get_type() != self.edge_types['E']:
                pq = pqubit(qubit_id_counter, y + origin[0], x + origin[1], self.edge_types['E'], '0')
                self.qubits.append(pq)
                self.parity_qubits.append(pq)

                qubit_id_counter += 1

        # West edge type
        init_y = 2
        x = 0

        for y in range(init_y, (dy * 2), 2):

            if self.xy_map[(y + origin[0], x+2 + origin[1])].get_type() != self.edge_types['W']:
                pq = pqubit(qubit_id_counter, y + origin[0], x + origin[1], self.edge_types['W'], '0')
                self.qubits.append(pq)
                self.parity_qubits.append(pq)

                qubit_id_counter += 1

        # Init internal qubit/coordinate maps
        self._reset_qubit_maps()

        # The data qubits are the first qubits initialized, so qubits [0,dx*dy - 1] are the data qubits

        # Go between the two z type edges

        # North to South
        if self.edge_types['N'] == 'z' and self.edge_types['S'] == 'z':
            data_qubit = pqubit_id_offset + int(ceil(dx / 2)) - 1
            for y in range(dy):
                self.z_obs.append(self.id_map[data_qubit + (y * dx)])

        # East to West
        elif self.edge_types['E'] == 'z' and self.edge_types['W'] == 'z':
            data_qubit = pqubit_id_offset + int(dy / 2) * dx
            for x in range(dx):
                self.z_obs.append(self.id_map[data_qubit + (x)])

        # Bad edges
        else:
            """raise ValueError(
                            'Bad edge types; (North, South) must have the same edge types, and (East, West) must have the same edge types.')"""
            warnings.warn(
                "Cannot initialize observable qubits list. Either the North/South or East/West edge types do not match\n"
                "Explicitly declare the observables with <object>.x_obs:list[int] = [...] or <object>.z_obs:list[int] = [...]")

        # Go between the two x type edges

        # North to South
        if self.edge_types['N'] == 'x' and self.edge_types['S'] == 'x':
            data_qubit = pqubit_id_offset + int(ceil(dx / 2)) - 1
            for y in range(dy):
                self.x_obs.append(self.id_map[data_qubit + (y * dx)])

        # East to West
        elif self.edge_types['E'] == 'x' and self.edge_types['W'] == 'x':
            data_qubit = pqubit_id_offset + int(dy / 2) * dx
            for x in range(dx):
                self.x_obs.append(self.id_map[data_qubit + (x)])

        # Bad edges
        else:
            """raise ValueError(
                'Bad edge types; (North, South) must have the same edge types, and (East, West) must have the same edge types.')"""
            warnings.warn(
                "Cannot initialize observable qubits list. Either the North/South or East/West edge types do not match\n"
                "Explicitly declare the observables with <object>.x_obs:list[int] = [...] or <object>.z_obs:list[int] = [...]")

        self.obs_qubits = self.x_obs if self.basis == 'x' else self.z_obs

        """for qubit in self.qubits:
            patch_specific_neighbors = super().get_pqubit_neighbors(qubit.get_id())
            for neighbor in patch_specific_neighbors:
                qubit.add_neighbor(neighbor.get_id(),qubit.get_type() if qubit.get_type() != 'd' else None)"""

    def get_qubit_coords_circuit(self):
        init_circuit = ''

        # INITIALIZE QUBIT COORDINATES
        for qubit in self.qubits:
            init_circuit += f"QUBIT_COORDS({qubit.get_pos_y()},{qubit.get_pos_x()}) {qubit.get_id()}\n"

        return init_circuit

    def get_init_circuit(self):

        init_circuit = ""

        init_circuit += self.get_qubit_coords_circuit()

        # Determine whether to reset the qubit in the X or Z basis
        z_resets = []
        x_resets = []
        for data_qubit in self.data_qubits:
            if data_qubit.get_init_state() == '0' or data_qubit.get_init_state() == '1':
                z_resets.append(data_qubit)
            elif data_qubit.get_init_state() == '+' or data_qubit.get_init_state() == '-':
                x_resets.append(data_qubit)
            else:
                raise ValueError(f'Bad initial state on qubit {data_qubit.get_id()}')

        # Apply excitations to data qubits (if any)
        z_excitations = []
        x_excitations = []
        for data_qubit in self.data_qubits:
            if data_qubit.get_init_state() == '1':
                z_excitations.append(data_qubit)
            if data_qubit.get_init_state() == '-':
                x_excitations.append(data_qubit)

        # TODO - incorporate initialization error on excitations
        if len(z_excitations) > 0:
            init_circuit += 'RZ '
            for data_qubit in z_excitations:
                init_circuit += f'{data_qubit.get_id()} '
            init_circuit += '\n'

        if len(x_excitations) > 0:
            init_circuit += 'RX '
            for data_qubit in x_excitations:
                init_circuit += f'{data_qubit.get_id()} '
            init_circuit += '\n'

        # TICK
        init_circuit += 'TICK\n'
        init_circuit += 'SHIFT_COORDS(0,0,1)\n'

        return init_circuit

    def get_quiescent_circuit(self):

        init_circuit = ''

        # Construct stabilizer circuits
        init_circuit += self.get_stabilizer_circuit(1)

        measurements_reversed = self.measurements[::-1]

        # Apply detectors to relevant stabilizers
        for parity_qubit in self.parity_qubits:

            if parity_qubit.get_exclude_detectors():
                continue

            if parity_qubit.get_type() == self.basis:
                init_circuit += f"DETECTOR({parity_qubit.get_pos_y()},{parity_qubit.get_pos_x()},0) "

                for m_index, m_id in enumerate(measurements_reversed):
                    if m_id == parity_qubit.get_id():
                        init_circuit += f'rec[{-1 * (1 + m_index)}] '

                init_circuit += '\n'

        return init_circuit

    def get_rep_circuit(self, num_rounds):
        """
        Shouldn't be used as a public method. Will do repetitions of patches with shared stabilizers, but this circuit
        will not be a proper merged/surgeried circuit.
        :param without_repeat:
        :return:
        """

        init_circuit = ""

        data_qubits = [qubit for qubit in self.data_qubits]

        parity_qubits = [qubit for qubit in self.parity_qubits]

        # REPEAT BLOCK
        if num_rounds > 1:
            init_circuit += f"REPEAT {num_rounds} " + "{\n"

        # DATA IDLE
        init_circuit += f'DEPOLARIZE1({self.uniform_error_rate}) '
        for data_qubit in data_qubits:
            init_circuit += f'{data_qubit.get_id()} '

        init_circuit += '\n'

        # Construct stabilizer circuits (includes resetting parity qubits prior to CNOTs, performing CNOTs, and measuring
        #   after CNOTs. This method DOES track measurements internally for us in self.measurements.
        init_circuit += self.get_stabilizer_circuit(num_rounds)

        # SHIFT COORDS
        init_circuit += 'SHIFT_COORDS(0,0,1)\n'

        # APPLY DETECTORS BETWEEN ROUNDS
        measurements_reversed = self.measurements[::-1]
        for parity_qubit in parity_qubits:

            # Skip qubits who do not require detectors
            if parity_qubit.get_exclude_detectors():
                continue

            init_circuit += f"DETECTOR({parity_qubit.get_pos_y()},{parity_qubit.get_pos_x()},0) "

            m_count = 0
            for m_index, m_id in enumerate(measurements_reversed):
                if m_id == parity_qubit.get_id():
                    init_circuit += (f'rec[{-1 * (m_index + 1)}] ')
                    m_count += 1

                if m_count >= 2:  # <-------Gets the 2nd latest measurement of the qubit
                    break

            if m_count < 2:
                raise ValueError(f'Couldn\'t find the prior round measurement for q{parity_qubit.get_id()}')

            init_circuit += '\n'

        # END REPEAT BLOCK
        if num_rounds > 1:
            init_circuit += "}\n"

        return init_circuit

    def get_terminating_circuit(self):
        init_circuit = ""

        # Measure all data qubits based on the multipatch/experiment basis
        if self.basis == 'z':
            init_circuit += f"X_ERROR({self.uniform_error_rate}) "
        elif self.basis == 'x':
            init_circuit += f'Z_ERROR({self.uniform_error_rate}) '
        else:
            raise ValueError("Bad basis while creating termination circuit.")

        for data_qubit in self.data_qubits:
            init_circuit += f'{data_qubit.get_id()} '

        init_circuit += '\n'

        if self.basis == 'z':
            init_circuit += 'M '
        elif self.basis == 'x':
            init_circuit += 'MX '
        else:
            raise ValueError("Bad basis while creating termination circuit.")

        for data_qubit in self.data_qubits:
            init_circuit += f'{data_qubit.get_id()} '
            self.measurements.append(data_qubit.get_id())

        init_circuit += '\n'
        measurements_reversed = self.measurements[::-1]

        for parity_qubit in self.get_parity_qubits():

            if parity_qubit.get_exclude_detectors():
                continue

            if parity_qubit.get_type() == self.basis:
                init_circuit += f'DETECTOR({parity_qubit.get_pos_y()},{parity_qubit.get_pos_x()},1) '

                for neighbor_qubit in self.get_pqubit_neighbors(parity_qubit.get_pos()):

                    if neighbor_qubit == -1:  # <---- -1 is the indicator for no neighbor
                        continue

                    init_circuit += f'rec[{-1 * (1 + measurements_reversed.index(neighbor_qubit.get_id()))}] '

                try:
                    init_circuit += f'rec[{-1 * (1 + measurements_reversed.index(parity_qubit.get_id()))}]\n'
                except ValueError:
                    warnings.warn(f"Unable to locate measurement record for {parity_qubit.get_id()}")

        init_circuit += '\n'

        if self.basis == 'x':
            init_circuit += 'OBSERVABLE_INCLUDE(0) '
            for qubit in self.x_obs:
                init_circuit += f'rec[{-1*(measurements_reversed.index(qubit.get_id())+1)}] '

        elif self.basis == 'z':
            init_circuit += 'OBSERVABLE_INCLUDE(0) '
            for qubit in self.z_obs:
                init_circuit += f'rec[{-1 * (measurements_reversed.index(qubit.get_id()) + 1)}] '

        else:
            raise ValueError("Observable unknown.")

        return init_circuit

    def get_full_circuit(self, num_rounds):
        self.measurements = []

        circ = ''

        circ += self.get_qubit_coords_circuit()
        circ += self.get_init_circuit()
        circ += self.get_quiescent_circuit()
        # for step_index in range(len(self.steps)):
        # self.next_step()
        for r in range(num_rounds):
            circ += self.get_rep_circuit(1)
        circ += self.get_terminating_circuit()
        return circ

    def get_west_edge_qubits(self):
        edge_qubits = []
        for qubit in self.qubits:
            if qubit.get_pos_x() <= 1:
                edge_qubits.append(qubit)

        return edge_qubits

    def get_east_edge_qubits(self):
        edge_qubits = []
        for qubit in self.qubits:
            if qubit.get_pos_x() >= (2 * self.dx + 1):
                edge_qubits.append(qubit)

        return edge_qubits

    def get_south_edge_qubits(self):
        edge_qubits = []
        for qubit in self.qubits:
            if qubit.get_pos_y() <= 1:
                edge_qubits.append(qubit)

        return edge_qubits

    def get_north_edge_qubits(self):
        edge_qubits = []
        for qubit in self.qubits:
            if qubit.get_pos_x() >= (2 * self.dx + 1):
                edge_qubits.append(qubit)

        return edge_qubits

    def get_west_edge_type(self):
        return self.edge_types['W']

    def get_east_edge_type(self):
        return self.edge_types['E']

    def get_north_edge_type(self):
        return self.edge_types['N']

    def get_south_edge_type(self):
        return self.edge_types['S']

    def get_metadata(self):
        return {
            'r': self.rounds,
            'b': self.basis,
            'dx': self.dx,
            'dy': self.dy,
            'psi_0': self.init_state,
            'p': self.uniform_error_rate,
            'origin': self.origin
        }

    def draw_patch(self, figsize=(8, 6), ax=None, draw_cx=True, show=True):
        from matplotlib import pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np

        ax_was_none = ax is None
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        x = []
        y = []
        ids = []
        for qubit in self.qubits:
            x.append(qubit.get_pos_x())
            y.append(qubit.get_pos_y())
            ids.append(qubit.get_id())
            ax.text(qubit.get_pos_x(), qubit.get_pos_y(), qubit.get_id())

        ax.scatter(x, y, label=self.id)

        for parity_qubit in self.parity_qubits:

            poly_points = []
            pos = parity_qubit.get_pos()
            try:
                poly_points.append(self.xy_map[(pos[0] - 1, pos[1] - 1)].get_pos()[::-1])
            except KeyError:
                pass
            try:
                poly_points.append(self.xy_map[(pos[0] - 1, pos[1] + 1)].get_pos()[::-1])
            except KeyError:
                pass
            try:
                poly_points.append(self.xy_map[(pos[0] + 1, pos[1] + 1)].get_pos()[::-1])
            except KeyError:
                pass
            try:
                poly_points.append(self.xy_map[(pos[0] + 1, pos[1] - 1)].get_pos()[::-1])
            except KeyError:
                pass

            if len(poly_points) < 3:
                poly_points.append(pos[::-1])
            poly = mpatches.Polygon(np.array(poly_points),
                                    facecolor='blue' if parity_qubit.get_type() == 'z' else 'red', edgecolor='black',
                                    zorder=-2, linewidth=0.15)
            ax.add_patch(poly)

            if draw_cx:
                neighbor_coordinates = []
                for neighbor in [qid for qid in self.get_pqubit_neighbors(parity_qubit.get_id())]:
                    if neighbor == -1:
                        continue
                    neighbor_coordinates.append(self.id_map[neighbor.get_id()].get_pos())
                for ni in range(1, len(neighbor_coordinates)):

                    prev_neighbor = neighbor_coordinates[ni - 1]

                    px_offset = -0.4
                    py_offset = -0.4
                    if prev_neighbor[0] < pos[0]:
                        py_offset = 0.4
                    if prev_neighbor[1] < pos[1]:
                        px_offset = 0.4

                    next_neighbor = neighbor_coordinates[ni]

                    nx_offset = -0.4
                    ny_offset = -0.4
                    if next_neighbor[0] < pos[0]:
                        ny_offset = 0.4
                    if next_neighbor[1] < pos[1]:
                        nx_offset = 0.4

                    ax.arrow(
                        prev_neighbor[1] + px_offset,
                        prev_neighbor[0] + py_offset,
                        (next_neighbor[1] + nx_offset) - (prev_neighbor[1] + px_offset),
                        (next_neighbor[0] + ny_offset) - (prev_neighbor[0] + py_offset),
                        facecolor='white',
                        edgecolor='black',
                        width=0.005,
                        head_width=0.15,
                        length_includes_head=True
                    )

        ax.grid()
        if show and ax_was_none:
            fig.show()
        if ax_was_none:
            return fig, ax

    """def get_pqubit_neighbors(self, coords_or_id: tuple[float,float] | int):
        if isinstance(coords_or_id, int):
            coords_or_id = self.id_map[coords_or_id].get_pos()

        neighbors = []

        for neighbor in self.xy_map[coords_or_id].get_neighbors():
            neighbors.append(neighbor)

        return neighbors"""

    def invert_stabilizers(self):
        for parity_qubit in self.parity_qubits:
            new_type = 'x' if parity_qubit.get_type() == 'z' else 'z'
            parity_qubit.set_type(new_type)
        self.reset_all_stabilizers()

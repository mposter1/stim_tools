from surf.pqubit import pqubit
class qubit_topology:

    x_stab_order: list[list[int]|int] = [
                                    [1,1],
                                    [-1,1],
                                    [1,-1],
                                    [-1,-1],
                                    -1
                                    ]
    z_stab_order: list[list[int]|int] = [
                                    [1,1],
                                    [1,-1],
                                    [-1,1],
                                    [-1,-1],
                                    -1
                                    ]

    stab_order_colors: list[str] = [
        'orange',
        'magenta',
        'cyan',
        'lime',
        'yellow',
        'beige',
        'salmon',
        'olive'
    ]

    def __init__(self, id, pqubit_id_offset, origin):
        self.qubits: list[pqubit] = []
        self.data_qubits: list[pqubit] = []
        self.parity_qubits: list[pqubit] = []
        self.obs_qubits: list[list[pqubit]] = []
        self.id: str = id
        self.pqubit_id_offset: int = pqubit_id_offset
        self.origin: tuple[float,float] = origin
        self.measurements: list = []
        self.uniform_error_rate: float = 0.0


        # Init internal qubit/coordinate maps
        self._reset_qubit_maps()

    def get_num_qubits(self):
        return len(self.xy_map)

    def get_id(self):
        return self.id

    def get_max_id(self):
        return max(self.id_map)

    def get_min_id(self):
        return min(self.id_map)

    def _reset_qubit_maps(self):
        id_map, xy_map = self.get_qubit_maps()
        self.id_map: dict[int, pqubit] = id_map
        self.xy_map: dict[tuple[float, float], pqubit] = xy_map

    def delete_pqubit(self, qubit_id):
        """
        Deletes the pqubit corresponding to the qubit_id from the topology, and removes any reference of the pqubit from
        the neighboring pqubits.
        :param qubit_id:
        :return:
        """
        qubit = self.id_map[qubit_id]
        del self.id_map[qubit_id]
        del self.xy_map[qubit.get_pos()]
        self.qubits.remove(qubit)
        if qubit.get_type() != 'd':
            self.parity_qubits.remove(qubit)
        else:
            self.data_qubits.remove(qubit)

    def add_pqubit(self, id: int, pos: tuple[float,float], qtype: str, is_obs=False):
        new_qubit = pqubit(id, pos[0], pos[1], qtype)
        self.qubits.append(new_qubit)
        if qtype == 'd':
            self.data_qubits.append(new_qubit)
        elif qtype == 'x' or type == 'z':
            self.parity_qubits.append(new_qubit)
        else:
            raise ValueError(f'Bad qubit type: {type}.')

        self.xy_map[new_qubit.get_pos()] = new_qubit
        self.id_map[new_qubit.get_id()] = new_qubit

    def add_existing_pqubit(self, new_qubit: pqubit):
        self.qubits.append(new_qubit)
        if new_qubit.get_type() == 'd':
            self.data_qubits.append(new_qubit)
        elif new_qubit.get_type() == 'x' or new_qubit.get_type() == 'z':
            self.parity_qubits.append(new_qubit)
        else:
            raise ValueError(f'Bad qubit type: {new_qubit.get_type()}.')

        self.xy_map[new_qubit.get_pos()] = new_qubit
        self.id_map[new_qubit.get_id()] = new_qubit

    def replace_pqubit(self, old_qubit_id: int, new_qubit: pqubit, additional_neighbors: list[int] = None, cx_types: list[str] = None):

        if self.id_map[old_qubit_id].get_pos() != new_qubit.get_pos():
            raise ValueError('Cannot replace a pqubit if it does not have the same location as the new pqubit.')

        old_qubit = self.id_map[old_qubit_id]
        del self.id_map[old_qubit_id]
        del self.xy_map[old_qubit.get_pos()]
        self.qubits.remove(old_qubit)

        new_qubit_id = new_qubit.get_id()
        self.id_map[new_qubit_id] = new_qubit
        self.xy_map[old_qubit.get_pos()] = new_qubit
        self.qubits.append(new_qubit)

        if old_qubit.get_type() != 'd':
            self.parity_qubits.remove(old_qubit)
            self.parity_qubits.append(new_qubit)
        else:
            self.data_qubits.remove(old_qubit)
            self.data_qubits.append(new_qubit)

        try:
            self.obs_qubits.remove(old_qubit)
            self.obs_qubits.append(new_qubit)
        except ValueError:
            pass

        if additional_neighbors is not None:

            if len(additional_neighbors) != len(cx_types):
                raise ValueError('Number of additional neighbors and number of CNOT types does not match in qubit replacement.')

            for neighbor, cx_type in zip(additional_neighbors, cx_types):
                new_qubit.add_neighbor(neighbor, cx_type)

    def mirror_across_y(self):

        # Find max x pos
        max_x = self.origin[1]
        for qubit in self.qubits:
            if qubit.get_pos_x() >= max_x:
                max_x = qubit.get_pos_x()

        # Absolute diff between qubit x pos and max x is the new x pos
        for qubit in self.qubits:
            new_x = self.origin[1] + abs(qubit.get_pos_x()-max_x)
            qubit.set_pos_x(new_x)

        self._reset_qubit_maps()
        #self.reset_all_stabilizers()

    def mirror_across_x(self):

        # Find max x pos
        max_y = self.origin[0]
        for qubit in self.qubits:
            if qubit.get_pos_y() >= max_y:
                max_y = qubit.get_pos_y()

        # Absolute diff between qubit x pos and max x is the new x pos
        for qubit in self.qubits:
            new_y = self.origin[0] + abs(qubit.get_pos_y()-max_y)
            qubit.set_pos_y(new_y)

        self._reset_qubit_maps()
        #self.reset_all_stabilizers()

    def get_qubit_maps(self):
        id_map = {}
        xy_map = {}
        for qubit in self.qubits:
            id_map[qubit.get_id()] = qubit
            xy_map[qubit.get_pos()] = qubit
        return id_map, xy_map

    def get_qubit_from_id(self, qubit_id):
        return self.id_map[qubit_id]

    def get_qubit_from_pos(self, qubit_pos):
        return self.xy_map[qubit_pos]

    def get_qubit_coordinates(self, qubit_id):
        return self.id_map[qubit_id].get_pos()

    def get_qubit_id_from_coordinates(self, coordinates):
        return self.xy_map[coordinates].get_id()

    def get_qubits(self):
        return self.qubits

    def get_data_qubits(self):
        return self.data_qubits

    def get_parity_qubits(self):
        return self.parity_qubits


    def get_pqubit_neighbors(self, coords_or_id: tuple[float,float] | int, type_filter=None):

        """
        Returns neighbors in a sequence corresponding to the ordering of the CX performed for stabilizers.
        Data qubits default to Z-type stabilizer ordering for their neighbors.
        """

        if isinstance(coords_or_id, int):
            coords_or_id = self.id_map[coords_or_id].get_pos()

        qtype = self.xy_map[coords_or_id].get_type()

        neighbors = []

        ordering = self.x_stab_order if qtype == 'x' else self.z_stab_order

        if self.xy_map[coords_or_id].custom_offsets is not None:
            ordering = self.xy_map[coords_or_id].custom_offsets

        for offset in ordering:

            if offset == -1:
                neighbors.append(-1)
                continue

            try:
                neighbor = self.xy_map[(coords_or_id[0] + offset[0], coords_or_id[1] + offset[1])]
                if type_filter is None or neighbor.get_type() == type_filter:
                    neighbors.append(neighbor)
                elif neighbor.get_type() != type_filter:
                    neighbors.append(-1)
            except KeyError:
                neighbors.append(-1)

        return neighbors


    def get_stabilizer_circuit(self, reps):
        """
        Returns the circuit that constructs all stabilizers in the topology.
        :return:
        """

        if reps > 1:
            raise NotImplementedError("Keep reps = 1 for now")

        init_circuit = ""

        # RESET THE Z STABILIZERS
        init_circuit += 'RZ '
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'z':
                init_circuit += f'{parity_qubit.get_id()} '

        init_circuit += '\n'

        # APPLY ERROR AFTER RESETTING Z STABILIZERS
        init_circuit += f'DEPOLARIZE1({self.uniform_error_rate}) '
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'z':
                init_circuit += f'{parity_qubit.get_id()} '

        init_circuit += '\n'

        # RESET THE X STABILIZERS
        init_circuit += 'RX '
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'x':
                init_circuit += f'{parity_qubit.get_id()} '

        init_circuit += '\n'

        # APPLY ERROR AFTER RESETTING THE X STABILIZERS
        init_circuit += f'DEPOLARIZE1({self.uniform_error_rate}) '
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'x':
                init_circuit += f'{parity_qubit.get_id()} '

        init_circuit += '\n'

        # TICK
        init_circuit += 'TICK\n'

        cx_rounds: dict[int,list[tuple[int,int]]] = {}

        # PERFORM CNOTS TO "ACTIVATE" STABILIZERS
        for parity_qubit in self.parity_qubits:
            neighbors = [neighbor for neighbor in self.get_pqubit_neighbors(parity_qubit.get_id(),type_filter='d')]
            parity_qubit_id = parity_qubit.get_id()
            for ni, neighbor in enumerate(neighbors):

                # For skipping a round in the timings/neighbor orderings
                if neighbor == -1:
                    continue

                neighbor_id = neighbor.get_id()

                assert(neighbor is not None)

                if ni not in cx_rounds:
                    cx_rounds[ni] = []

                # x: parity targets data
                if parity_qubit.get_type() == 'x':
                    cx_rounds[ni].append((parity_qubit_id, neighbor_id))

                # z: data targets parity
                elif parity_qubit.get_type() == 'z':
                    cx_rounds[ni].append((neighbor_id, parity_qubit_id))

                else:
                    raise ValueError(f'Unknown CX type between qubits {parity_qubit_id} and {neighbor}')

        rounds = cx_rounds.keys()
        sorted_rounds = sorted(rounds)
        for round_num in sorted_rounds:

            init_circuit += 'CX '

            cxs = cx_rounds[round_num]

            for cx in cxs:

                init_circuit += f'{cx[0]} {cx[1]} '

            init_circuit += '\n'

            init_circuit += f'DEPOLARIZE2({self.uniform_error_rate}) '

            for cx in cxs:

                init_circuit += f'{cx[0]} {cx[1]} '

            init_circuit += '\n'

        # TICK
        init_circuit += 'TICK\n'

        # APPLY ERROR PRIOR TO MEASUREMENT OF Z STABILIZERS
        init_circuit += f'DEPOLARIZE1({self.uniform_error_rate}) '
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'z':
                init_circuit += f'{parity_qubit.get_id()} '

        init_circuit += '\n'

        # MEASURE Z STABILIZERS
        init_circuit += 'M '
        z_measurements = []
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'z':
                init_circuit += f'{parity_qubit.get_id()} '
                z_measurements.append(parity_qubit.get_id())

        # TRACK MEASUREMENTS
        for rep in range(reps):
            for measurement in z_measurements:
                self.measurements.append(measurement)

        init_circuit += '\n'

        # APPLY ERROR PRIOR TO MEASUREMENT OF X STABILIZERS
        init_circuit += f'DEPOLARIZE1({self.uniform_error_rate}) '
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'x':
                init_circuit += f'{parity_qubit.get_id()} '

        init_circuit += '\n'

        # MEASURE THE X STABILIZERS (RETURNING FROM SUPERPOSITION)
        init_circuit += 'MX '
        x_measurements = []
        for parity_qubit in self.parity_qubits:
            if parity_qubit.get_type() == 'x':
                init_circuit += f'{parity_qubit.get_id()} '
                x_measurements.append(parity_qubit.get_id())

        # TRACK MEASUREMENTS
        for rep in range(reps):
            for measurement in x_measurements:
                self.measurements.append(measurement)

        init_circuit += '\n'

        return init_circuit

    def get_select_stabilizer_circuit(self, pqubits: list[pqubit]):
        """
        Returns the circuit that constructs all stabilizers in the topology.
        :return:
        """

        raise NotImplementedError
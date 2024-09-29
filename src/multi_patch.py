import warnings

from surf.src.patch import surf_patch
from surf.src.qubit_topology import qubit_topology
from surf.src.pqubit import pqubit
from surf.src.step import surf_step
class multi_surf_patch(qubit_topology):

    def __init__(self, patches: list[surf_patch], uniform_error_rate: float, experiment_basis: str, merge_neighbors: bool = True):
        super().__init__(''.join(f'{patch.get_id()}, ' for patch in patches),
                         0,
                         (0, 0))
        self.patches: list[surf_patch] = patches
        self.edges: dict[str, dict[str, list[pqubit]]] = {}
        self.qubits: list[pqubit] = []
        self.data_qubits: list[pqubit] = []
        self.parity_qubits: list[pqubit] = []
        self.id: str = ''.join(f'{patch.get_id()}, ' for patch in patches)
        self.uniform_error_rate: float = uniform_error_rate
        self.q_id_to_patch: dict[int, surf_patch] = {}

        self.basis = experiment_basis

        for patch in patches:
            self.edges[patch.get_id()] = {
                "N": patch.get_north_edge_qubits(),
                "S": patch.get_south_edge_qubits(),
                "E": patch.get_east_edge_qubits(),
                "W": patch.get_west_edge_qubits(),
            }

        self.xy_map: dict[[float,float], pqubit] = {}
        self.id_map: dict[int, pqubit] = {}
        for patch in patches:
            for qubit in patch.get_qubits():
                if qubit.get_id() in self.id_map and qubit.get_pos() != self.id_map[qubit.get_id()].get_pos():
                    raise ValueError(f'Patches have conflicting physical qubits: q{qubit.get_id()} @ {qubit.get_pos()} and {self.id_map[qubit.get_id()].get_pos()}.')

                if merge_neighbors and qubit.get_pos() in self.xy_map:
                    warnings.warn(f'Overlapping physical qubits in multi-patch:\n'
                                  f' q{qubit.get_id()} in patch {patch.get_id()} with q{self.xy_map[(qubit.get_pos())].get_id()} @ {qubit.get_pos()}.\n'
                                  f'Replacing... q{qubit.get_id()} in {patch.get_id()} with q{self.xy_map[(qubit.get_pos())].get_id()}')

                    patch.replace_pqubit(qubit.get_id(), self.xy_map[(qubit.get_pos())],
                                         #additional_neighbors = qubit.get_neighbors(),
                                         #cx_types = [self.xy_map[(qubit.get_pos())].get_type() for q in self.xy_map[(qubit.get_pos())].get_neighbors()]
                                         )

                    """new_neighbors = qubit.get_neighbors()
                    old_neighbors = self.xy_map[(qubit.get_pos())].get_neighbors()
                    self.neighbor_reset_buffer.add(self.xy_map[(qubit.get_pos())].get_id())
                    for new_neighbor in new_neighbors:

                        if new_neighbor is None:
                            continue

                        self.neighbor_reset_buffer.add(new_neighbor)
                        self.xy_map[(qubit.get_pos())].add_neighbor(new_neighbor, qubit.get_cx_type(new_neighbor))

                    for old_neighbor in old_neighbors:

                        if old_neighbor is None:
                            continue

                        self.neighbor_reset_buffer.add(old_neighbor)"""
                    continue


                self.qubits.append(qubit)
                #if qubit.get_pos() in self.xy_map:
                #    self.xy_map[qubit.get_pos()] = [self.xy_map[qubit.get_pos()].get_id(), qubit.get_id()]
                #else:
                self.xy_map[qubit.get_pos()] = qubit
                self.id_map[qubit.get_id()] = qubit
                if qubit.get_type() == 'd':
                    self.data_qubits.append(qubit)
                elif qubit.get_type() == 'x' or qubit.get_type() == 'z':
                    self.parity_qubits.append(qubit)
                else:
                    raise ValueError('Bad physical qubit type in multi-patch.')

        self.exclude_detectors: dict[int,pqubit] = {}

        for patch in self.patches:
            for qubit in patch.get_qubits():
                self.q_id_to_patch[qubit.get_id()] = patch

        self.steps: list[surf_step] = []
        self.current_step: int = -1

    def draw_patch(self, figsize=(8,6), ax=None, show=True, draw_merged=True, draw_cx=True):
        from matplotlib import pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
        from math import atan2

        drawn_colors = []

        def argsort(seq):
            # http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
            # by unutbu
            # https://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python
            # from Boris Gorelik
            return sorted(range(len(seq)), key=seq.__getitem__)

        def rotational_sort(list_of_yx_coords, centre_of_rotation_yx_coord, clockwise=True):
            # https: // stackoverflow.com / questions / 67735537 / how - to - sort - coordinates - in -python - in -a - clockwise - direction
            # by Thomas Kimber
            cy, cx = centre_of_rotation_yx_coord
            angles = [atan2(y - cy,x - cx) for y, x in list_of_yx_coords]
            indices = argsort(angles)
            if clockwise:
                return [list_of_yx_coords[i] for i in indices]
            else:
                return [list_of_yx_coords[i] for i in indices[::-1]]

        ax_was_none = ax is None
        if ax is None:
            fig,ax = plt.subplots(1,1,figsize=figsize)

        if not draw_merged:
            for patch in self.patches:
                x = []
                y = []
                ids = []
                for qubit in patch.qubits:
                    x.append(qubit.get_pos_x())
                    y.append(qubit.get_pos_y())
                    ids.append(qubit.get_id())
                    ax.text(qubit.get_pos_x(),qubit.get_pos_y(),qubit.get_id())

                ax.scatter(x,y,label=patch.id)

                for parity_qubit in patch.parity_qubits:

                    poly_points = []
                    pos = parity_qubit.get_pos()
                    neighbors = self.get_pqubit_neighbors(pos,type_filter='d')

                    """try:
                        poly_points.append(patch.xy_map[(pos[0]-1,pos[1]-1)].get_pos()[::-1])
                    except KeyError:
                        pass
                    try:
                        poly_points.append(patch.xy_map[(pos[0]-1,pos[1]+1)].get_pos()[::-1])
                    except KeyError:
                        pass
                    try:
                        poly_points.append(patch.xy_map[(pos[0]+1,pos[1]+1)].get_pos()[::-1])
                    except KeyError:
                        pass
                    try:
                        poly_points.append(patch.xy_map[(pos[0]+1,pos[1]-1)].get_pos()[::-1])
                    except KeyError:
                        pass"""

                    for neighbor in neighbors:
                        if neighbor == -1:
                            continue
                        poly_points.append(neighbor.get_pos())

                    if len(poly_points) < 3:
                        poly_points.append(pos)
                    else:
                        poly_points = rotational_sort(poly_points,parity_qubit.get_pos())
                    poly = mpatches.Polygon(np.array([point[::-1] for point in poly_points]), facecolor='blue' if parity_qubit.get_type() == 'z' else 'red', edgecolor='black', zorder=-2, linewidth=0.15)
                    ax.add_patch(poly)
        else:
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
                neighbors = self.get_pqubit_neighbors(pos,type_filter='d')
                """try:
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
                    pass"""

                for neighbor in neighbors:
                    if neighbor == -1:
                        continue
                    poly_points.append(neighbor.get_pos())

                if len(poly_points) < 3:
                    poly_points.append(pos)
                else:
                    poly_points = rotational_sort(poly_points, parity_qubit.get_pos())
                poly = mpatches.Polygon(np.array([point[::-1] for point in poly_points]),
                                        facecolor='blue' if parity_qubit.get_type() == 'z' else 'red',
                                        edgecolor='black', zorder=-2, linewidth=0.15)
                ax.add_patch(poly)

                if draw_cx:

                    neighbors = self.get_pqubit_neighbors(parity_qubit.get_id(),type_filter='d')
                    for ni,neighbor in enumerate(neighbors):

                        if neighbor == -1:
                            # Skip
                            continue

                        x_offset = -0.15
                        y_offset = -0.15
                        neighbor_pos = neighbor.get_pos()
                        if neighbor_pos[0] < pos[0]:
                            y_offset = 0.15
                        if neighbor_pos[1] < pos[1]:
                            x_offset = 0.15

                        label = f'{ni}' if self.stab_order_colors[ni] not in drawn_colors else None

                        if label is not None:
                            drawn_colors.append(self.stab_order_colors[ni])

                        ax.scatter([neighbor_pos[1] + x_offset],
                                   [neighbor_pos[0] + y_offset],
                                   c=self.stab_order_colors[ni],
                                   marker='D',
                                   label=label,
                                   edgecolor='black'
                                   )


                    sanitized_neighbors = [neighbor for neighbor in neighbors if neighbor != -1]
                    neighbor_coordinates = [neighbor.get_pos() for neighbor in sanitized_neighbors]

                    for ni in range(1,len(neighbor_coordinates)):

                        prev_neighbor = neighbor_coordinates[ni-1]

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

                        #cx_orderings[neighbor_indices[ni]].append(arrow)

        ax.grid()
        ax.legend()#cx_orderings,[str(i) for i in range(len(cx_orderings))],handler_map={tuple: HandlerTuple(ndivide=None)})
        if show and ax_was_none:
            fig.show()
        if ax_was_none:
            return fig,ax


    def get_qubit_coords_circuit(self):
        init_circuit = ''

        # INITIALIZE QUBIT COORDINATES
        for qubit in self.qubits:
            init_circuit += f"QUBIT_COORDS({qubit.get_pos_y()},{qubit.get_pos_x()}) {qubit.get_id()}\n"

        for step in self.steps:
            for qubit in step.get_pqubits():
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

        # We DO include step qubits here so we don't have to initialize them later
        for step in self.steps:
            for data_qubit in step.get_pqubits():
                if data_qubit.get_type() == 'd':
                    if data_qubit.get_init_state() == '0' or data_qubit.get_init_state() == '1':
                        z_resets.append(data_qubit)
                    elif data_qubit.get_init_state() == '+' or data_qubit.get_init_state() == '-':
                        x_resets.append(data_qubit)
                    else:
                        raise ValueError(f'Bad initial state on qubit {data_qubit.get_id()}')

        # Apply resets accordingly
        if len(z_resets) > 0:
            init_circuit += 'RZ '
            for data_qubit in z_resets:
                init_circuit += f'{data_qubit.get_id()} '
            init_circuit += '\n'
            init_circuit += f'X_ERROR({self.uniform_error_rate}) '
            for data_qubit in z_resets:
                init_circuit += f'{data_qubit.get_id()} '
            init_circuit += '\n'

        # TICK
        init_circuit += 'TICK\n'
        init_circuit += 'SHIFT_COORDS(0,0,1)\n'

        if len(x_resets) > 0:
            init_circuit += 'RX '
            for data_qubit in x_resets:
                init_circuit += f'{data_qubit.get_id()} '
            init_circuit += '\n'

            init_circuit += f'Z_ERROR({self.uniform_error_rate}) '
            for data_qubit in x_resets:
                init_circuit += f'{data_qubit.get_id()} '
            init_circuit += '\n'

        """# Gather parity qubits across steps
        parity_qubits = [qubit for qubit in self.parity_qubits]
        for step in self.steps:
            for qubit in step.get_pqubits():
                if qubit.get_type() != 'd':
                    parity_qubits.append(qubit)

        # Reset all parity qubits in Z basis
        init_circuit += 'RZ '
        for parity_qubit in parity_qubits:
            init_circuit += f'{parity_qubit.get_id()} '
        init_circuit += '\n'
        init_circuit += f'X_ERROR({self.uniform_error_rate}) '
        for parity_qubit in parity_qubits:
            init_circuit += f'{parity_qubit.get_id()} '
        init_circuit += '\n'

        # TICK
        init_circuit += 'TICK\n'"""

        # Apply excitations to data qubits (if any)
        z_excitations = []
        x_excitations = []
        for data_qubit in self.data_qubits:
            if data_qubit.get_init_state() == '1':
                z_excitations.append(data_qubit)
            if data_qubit.get_init_state() == '-':
                x_excitations.append(data_qubit)

        for step in self.steps:
            for data_qubit in step.get_pqubits():
                if data_qubit.get_type() == 'd':
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

    def get_rep_circuit(self, num_rounds, add_to_measurements_before_detectors=None):
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

        if add_to_measurements_before_detectors is not None:
            init_circuit += '\nMX '
            for qi in self.steps[self.current_step].get_pqubits():
                if qi.get_type() == 'd':
                    init_circuit += f'{qi.get_id()} '
                    self.measurements.append(qi.get_id())
            init_circuit += '\n'

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
            for m_index,m_id in enumerate(measurements_reversed):
                if m_id == parity_qubit.get_id():
                    init_circuit += (f'rec[{-1*(m_index+1)}] ')
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

                for neighbor_qubit in self.get_pqubit_neighbors(parity_qubit.get_pos(),type_filter='d'):

                    if neighbor_qubit == -1:  # <---- -1 is the indicator for no neighbor
                        continue

                    init_circuit += f'rec[{-1 * (1 + measurements_reversed.index(neighbor_qubit.get_id()))}] '

                try:
                    init_circuit += f'rec[{-1 * (1 + measurements_reversed.index(parity_qubit.get_id()))}]\n'
                except ValueError:
                    warnings.warn(f"Unable to locate measurement record for {parity_qubit.get_id()}")

        init_circuit += '\n'

        return init_circuit

    def get_full_circuit(self, num_rounds):
        self.measurements = []

        circ = ''

        circ += self.get_qubit_coords_circuit()
        circ += self.get_init_circuit()
        circ += self.get_quiescent_circuit()
        for r in range(num_rounds-1):
            circ += self.get_rep_circuit(1)
        circ += self.get_terminating_circuit()
        return circ


    def invert_stabilizers(self):
        for patch in self.patches:
            patch.invert_stabilizers()

    def add_step(self, step: surf_step):
        self.steps.append(step)

    def next_step(self):
        """
        Updates the circuit elements to remove any qubits from prior steps and contain only the temporary qubits from
        the next step within the multi-patch sequence.
        """

        if self.current_step > -1:
            for qubit in self.steps[self.current_step].get_pqubits():

                # Don't delete the qubit if it's a part of an existing patch
                for patch in self.patches:
                    if qubit in patch.get_qubits():
                        continue

                self.delete_pqubit(qubit.get_id())

        self.current_step += 1

        for qubit in self.steps[self.current_step].get_pqubits():
            if qubit not in self.qubits:
                self.add_existing_pqubit(qubit)

    def draw_layers(self):
        """
        Draws all layers of the circuit.
        :return:
        """
        import matplotlib.pyplot as plt
        import numpy as np

        cx_count_colors = ['orange','magenta','cyan','lime',]

        fig,ax = plt.subplots(len(self.steps), 4, figsize=(8*4,6*len(self.steps)))
        if isinstance(ax,np.ndarray):
            ax = ax.flat
        for i in range(0,len(ax),4):
            init_circ = self.get_init_circuit()
            quis_circ = self.get_quiescent_circuit(i,self.steps[i].get_basis())
            rep_circ = self.get_rep_circuit(self.steps[i].get_rounds(), i)
            term_circ = self.get_terminating_circuit(i)

            cx_counter = 0
            cx_dict = {}
            axis = ax[0]
            for qubit in self.qubits:
                pos = qubit.get_pos()
                axis.scatter(pos[1],pos[0])
                axis.text(pos[1],pos[0],str(qubit.get_id()))
            for line in init_circ.splitlines():
                items = line.split()
                if len(items) == 0:
                    continue
                if items[0] == 'CX':
                    cxs = []
                    for k in range(1,len(items),2):
                        cxs.append((items[k],items[k+1]))
                        cx_dict[cx_counter] = cxs
                    cx_counter += 1
            for cx_num in cx_dict:
                cxs = cx_dict[cx_num]
                for cx in cxs:
                    q1 = self.id_map[int(cx[0])].get_pos()
                    q2 = self.id_map[int(cx[1])].get_pos()
                    axis.arrow(q1[1],q1[0],q2[1]-q1[1],q2[0]-q1[0],length_includes_head=True,head_width=0.15)
                    axis.text((q1[1] + q2[1])/2,(q1[0]+q2[0])/2,str(cx_num),c=cx_count_colors[cx_num], zorder=4)

            cx_counter = 0
            cx_dict = {}
            axis = ax[2]
            for qubit in self.qubits:
                pos = qubit.get_pos()
                axis.scatter(pos[1], pos[0])
                axis.text(pos[1],pos[0],str(qubit.get_id()))
            for line in rep_circ.splitlines():
                items = line.split()
                if items[0] == 'CX':
                    cxs = []
                    for k in range(1, len(items), 2):
                        cxs.append((items[k], items[k + 1]))
                        cx_dict[cx_counter] = cxs
                    cx_counter += 1
            for cx_num in cx_dict:
                cxs = cx_dict[cx_num]
                for cx in cxs:
                    q1 = self.id_map[int(cx[0])].get_pos()
                    q2 = self.id_map[int(cx[1])].get_pos()
                    axis.arrow(q1[1], q1[0], q2[1] - q1[1], q2[0] - q1[0], length_includes_head=True,head_width=0.15)
                    axis.text((q1[1] + q2[1])/2,(q1[0]+q2[0])/2, str(cx_num), c=cx_count_colors[cx_num], zorder=4)

        fig.show()
        return fig,ax
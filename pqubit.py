class pqubit:
    """
    Physical qubit abstraction for surface-code patches
    """
    def __init__(self, id: int, pos_y: float, pos_x: float, type: str, init_state: str, exclude_detectors: bool = False):
        self.id: int = id
        self.pos_y: float = pos_y
        self.pos_x: float = pos_x
        self.pos: tuple[float,float] = (pos_y, pos_x)
        self.type: str = type
        self.init_state: str = init_state
        self.custom_offsets = None
        self.exclude_detectors: bool = exclude_detectors

    def get_pos_x(self):
        return self.pos_x

    def get_pos_y(self):
        return self.pos_y

    def get_pos(self):
        return self.pos

    def get_id(self):
        return self.id

    def get_type(self):
        return self.type

    def set_pos_x(self, val):
        self.pos_x = val
        self.pos = (self.pos_y,self.pos_x)

    def set_pos_y(self, val):
        self.pos_y = val
        self.pos = (self.pos_y,self.pos_x)

    def set_pos(self, val):
        self.pos = val

    def set_id(self, val):
        self.id = val

    def set_type(self, val):
        self.type = val


    def set_init_state(self, val: str):
        self.init_state = val

    def get_init_state(self):
        return self.init_state

    def get_exclude_detectors(self):
        return self.exclude_detectors

    def set_exclude_detectors(self, val: bool):
        self.exclude_detectors = val

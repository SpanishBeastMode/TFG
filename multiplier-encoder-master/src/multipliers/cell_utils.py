import math

DELTA_k = 4
LEN_K = 12
STEP_z = 3
class CellWrapper(object):
    def __init__(self, cell, z_i):
        self.cell = cell
        self.z_i = z_i
        self.MOVE = {
            45: {
                "delta_k": (DELTA_k, DELTA_k),
                "delta_z": (1, 1),
                "up": (self._plus_k, self._minus_k, self._minus_z, self._plus_z),
                "down": (self._minus_k, self._plus_k, self._plus_z, self._minus_z)
            },
            60: {
                "delta_k": (DELTA_k, LEN_K - DELTA_k),
                "delta_z": (STEP_z - 1, 1),
                "up": (self._minus_k, self._minus_k, self._minus_z, self._minus_z),
                "down": (self._plus_k, self._plus_k, self._plus_z, self._plus_z)
            },
            30: {
                "delta_k": (LEN_K - DELTA_k, DELTA_k),
                "delta_z": (1, STEP_z - 1),
                "up": (self._minus_k, self._minus_k, self._minus_z, self._minus_z),
                "down": (self._plus_k, self._plus_k, self._plus_z, self._plus_z)
            },
            0: {
                "delta_k": (LEN_K, 0),
                "delta_z": (0, STEP_z),
                "up": (self._minus_k, self._plus_k, self._plus_z, self._minus_z),
                "down": (self._plus_k, self._plus_k, self._plus_z, self._plus_z)
            },
            90: {
                "delta_k": (0, LEN_K),
                "delta_z": (STEP_z, 0),
                "up": (self._minus_k, self._minus_k, self._minus_z, self._minus_z),
                "down": (self._plus_k, self._plus_k, self._plus_z, self._plus_z)
            }
        }

    def next(self, dir, up_or_down):
        V, H = self.cell[0:4], self.cell[4:8]
        z_i_v, z_i_h = self.z_i

        k_v_shift, k_h_shift, z_v_shift, z_h_shift = self.MOVE[dir][up_or_down]
        delta_k_v, delta_k_h = self.MOVE[dir]["delta_k"]
        delta_z_v, delta_z_h = self.MOVE[dir]["delta_z"]

        V_up = [tuple([v] + k_v_shift(w, k, delta_k_v) + [z_v_shift(z, z_i_v, delta_z_v)[0]]) for (v, w, k, z) in V]
        H_up = [tuple([v] + k_h_shift(w, k, delta_k_h) + [z_h_shift(z, z_i_h, delta_z_h)[0]]) for (v, w, k, z) in H]
        z_v, z_h = V[0][3], H[0][3]
        z_i_v, z_i_h = z_v_shift(z_v, z_i_v, delta_z_v)[1], z_h_shift(z_h, z_i_h, delta_z_h)[1]
        return CellWrapper(tuple(V_up + H_up), (z_i_v, z_i_h))

    def _plus_k(self, w, k, delta_k):
        return [w + math.floor((k + delta_k) / LEN_K), (k + delta_k) % LEN_K]

    def _minus_k(self, w, k, delta_k):
        return [w - (1 - (math.floor((k + (LEN_K - delta_k)) / LEN_K))), (k + (LEN_K - delta_k)) % LEN_K]
    
    def _minus_z(self, z, z_i, delta_z):
        return [z - (1 - math.floor((z_i + (STEP_z - delta_z)) / STEP_z)), (z_i + (STEP_z - delta_z)) % STEP_z]

    def _plus_z(self, z, z_i, delta_z):
        return [z + math.floor((z_i + delta_z) / STEP_z), (z_i + delta_z) % STEP_z]

def multiplier_origin(m=16):
    cell = tuple([(0, 6, k, 0) for k in range(0, 4)] + [(1, 0, k, 5) for k in range(8, 12)])
    z_i = (0, 2)
    return CellWrapper(cell, z_i)

def multiplier_origin_diagonal_45_shifted_2units(m=16):
    cell = tuple([(0, 6-1, k, 0) for k in range(4, 8)] + [(1, 0+1, k, 5) for k in range(4, 8)])
    z_i = (2, 0)
    return CellWrapper(cell, z_i)

def multiplier_origin_right_shifted_3units(m=16, shift=3):
    cell = tuple([(0, 6 + shift, k, 0) for k in range(0, 4)] + [(1, 0, k, 5 + shift) for k in range(8, 12)])
    z_i = (0, 2)
    return CellWrapper(cell, z_i)

def origin_4_max_rect(m=16):
    cell = tuple([(0, 8, k_h, 0) for k_h in range(0, 4)] + [(1, 0, k_v, 7) for k_v in range(8, 12)])
    z_i = (0, 2)
    return CellWrapper(cell, z_i)
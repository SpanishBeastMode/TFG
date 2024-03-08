import abc

class Construct_subgraph(abc.ABC):

    @abc.abstractmethod
    def prepare_base_tiles(self, base_cellwrpr, shared_direction, pegasus_to_linear):
        pass

    @abc.abstractmethod
    def embed(self, cfa_plmt, shared_from, tn=8):
        pass

@Construct_subgraph.register
class Construct_CFAgraph_for_Multiplier_ver2():
    def __init__(self, base_cellwrpr, shared_direction, pegasus_to_linear):
        self.tile_dict = self.prepare_base_tiles(base_cellwrpr, shared_direction, pegasus_to_linear)

    def prepare_base_tiles(self, base_cellwrpr, shared_direction, pegasus_to_linear):
        base_tile = {'base': base_cellwrpr.cell}
        start = base_cellwrpr
        extra_tiles = {k: start.next(*v).cell for k, v in shared_direction.items()}
        assert start == base_cellwrpr
        tile_dict = {**base_tile, **extra_tiles}

        tile_dict = {key: [pegasus_to_linear(qbt) for qbt in cell] for key, cell in tile_dict.items()}
        return tile_dict

    def embed(self, cfa_plmt, shared_from, tn=8):
        graph = {}
        for v, p in cfa_plmt.items():
            if p < tn:
                graph[v] = self.tile_dict['base'][p]
            else:
                rp = cfa_plmt[shared_from[v]]
                if p == tn:
                    graph[v] = self.tile_dict['out'][rp]
                else:
                    graph[v] = self.tile_dict['carry'][rp]
        return graph

@Construct_subgraph.register
class Construct_CFAgraph_for_Multiplier_ver3(Construct_CFAgraph_for_Multiplier_ver2):
    def __init__(self, base_cellwrpr, shared_direction, pegasus_to_linear):
        self.tile_dict = self.prepare_base_tiles(base_cellwrpr, shared_direction, pegasus_to_linear)

    def embed(self, cfa_plmt, shared_from, tn=8):
        graph = {}
        for v, p in cfa_plmt.items():
            if p < tn:
                graph[v] = self.tile_dict['base'][p]
            else:
                sp = cfa_plmt[shared_from[v]]
                sw = v.split('out')[0]
                if not sw:
                    graph[v] = self.tile_dict[v][sp]
                elif sw == 'in1_':
                    graph[v] = self.tile_dict['A'][sp]
                elif sw in ['c_', 'enable_']:
                    graph[v] = self.tile_dict['carry'][sp]
        return graph

@Construct_subgraph.register
class Construct_chain_for_adder_carry():
    def __init__(self, base_cellwrpr, shared_direction, pegasus_to_linear):
        self.next_tile, self.next_next_tile = self.prepare_base_tiles(base_cellwrpr, shared_direction, pegasus_to_linear)

    def prepare_base_tiles(self, base_cellwrpr, propg_direction, pegasus_to_linear):
        next_cellwrpr =  base_cellwrpr.next(*propg_direction['B'])
        next_tile = [pegasus_to_linear(qbt) for qbt in next_cellwrpr.cell]
        next_next_tile = [pegasus_to_linear(qbt) for qbt in next_cellwrpr.next(*propg_direction['out']).cell]
        return next_tile, next_next_tile

    def embed(self, cfa_plmt, shared_from):
        # p = cfa_plmt[shared_from['c_out']]
        p = cfa_plmt['c_in']
        ch = []
        endp = 7
        if p in range(2, 4):
            path = [(p, endp)]
        elif p in range(4, 6):
            via = 3
            path = [(p, via), (via, endp)]
        ch += [tuple([self.next_tile[x] for x in ft]) for ft in path]

        p = cfa_plmt['in2']
        edg = (self.next_tile[endp], self.next_next_tile[cfa_plmt['in2']])
        ch.append(edg)
        return ch


@Construct_subgraph.register
class Construct_chain_for_adder_carry_bit_for_all_link_version(Construct_chain_for_adder_carry):
    def __init__(self, base_cellwrpr, shared_direction, pegasus_to_linear):
        self.base_tile = [pegasus_to_linear(qbt) for qbt in base_cellwrpr.cell]
        self.next_tile, self.next_next_tile = self.prepare_base_tiles(base_cellwrpr, shared_direction, pegasus_to_linear)

    def embed(self, cfa_plmt, shared_from):
        ch = [(self.base_tile[cfa_plmt['c_out']], self.next_tile[cfa_plmt['c_in']])]

        ch += super().embed(cfa_plmt, shared_from)
        return ch

def construct_horizontal_chains(cfa_graphs, vs=['enable']*2):
    chains = []
    for i, adder_graph in enumerate(cfa_graphs):
        chains.append([])
        for cfa_graph, next_cfa_graph in zip(adder_graph[:-1], adder_graph[1:]):
            edg = (cfa_graph[vs[0]], next_cfa_graph[vs[1]])
            chains[i].append(edg)
    return chains

def construct_vertical_chains(cfa_graphs, v='in1'):
    def __bridge(fi, ti, j, cfa_graphs=cfa_graphs):
        return tuple([cfa_graphs[i][j][v] for i in [fi, ti]])

    l, w = len(cfa_graphs[0]), len(cfa_graphs)
    chains = []
    for j in range(l):
        ch = [__bridge(i, i+1, j) for i in range(w - 1)]
        chains.append(ch)
    return chains

def construct_chains_for_cfa_out(cfa_graphs, vs=('out', 'in2')):
    l, w = len(cfa_graphs[0]), len(cfa_graphs)
    chains = []
    for i in range(w-1):
        chains.append([])
        for j in range(1, l):
            cfas = cfa_graphs[i][j], cfa_graphs[i+1][j-1]
            edg = tuple([cfa[v] for v, cfa in zip(vs, cfas)])
            chains[i].append(edg)
    return chains

def construct_chains_for_cfa_carry(cfa_graphs, vs=('c_out', 'c_in')):
    l, w = len(cfa_graphs[0]), len(cfa_graphs)
    chains = []
    for i in range(w):
        chains.append([])
        for j in range(l-1):
            cfas = cfa_graphs[i][j], cfa_graphs[i][j+1]
            edg = tuple([cfa[v] for v, cfa in zip(vs, cfas)])
            chains[i].append(edg)
    return chains

from __future__       import division, print_function

from itertools        import combinations, combinations_with_replacement
from collections      import deque, namedtuple, defaultdict, Counter

from vgraph.norm      import NormalizedLocus


RefNode = namedtuple('RefNode', 'start stop seq')
AltNode = namedtuple('AltNode', 'start stop seq locus')
VariantGraph = namedtuple('VariantGraph', 'graph start_node zygosity_constraints phase_constraints')


class RefCache(object):
    def __init__(self, ref):
        self.ref = ref
        self.cache = {}

    def __getitem__(self, interval):
        node = self.cache.get(interval)
        if node is None:
            start, stop = interval
            node = self.cache[interval] = RefNode(start, stop, self.ref[start:stop])
        return node


def is_valid_path(vg, path):
    if vg.phase_constraints:
        pathset = set(p for p in path if isinstance(p, AltNode))
        for phaseset in vg.phase_constraints.itervalues():
            n = len(pathset & phaseset)
            if n != 0 and n != len(phaseset):
                return False

    return True


def is_valid_geno(vg, paths):
    observed_zygosity = Counter(p for path in paths for p in path if isinstance(p, AltNode))
    return vg.zygosity_constraints == observed_zygosity


def get_path_seq(path):
    return ''.join(p.seq for p in path)


def make_vargraph(ref, start, stop, loci, name):
    start_node, stop_node = RefNode(start, start, ''), RefNode(stop + 1, stop + 1, '')
    vg = VariantGraph(defaultdict(set), start_node, defaultdict(int), defaultdict(set))
    ref_cache = RefCache(ref)

    pos = start_node.start
    new_layer = [start_node]

    for locus in sorted(loci, key=NormalizedLocus.left_order_key):
        sample = locus.record.samples[name]
        indices = sample.allele_indices
        left = locus.left

        if left.start >= pos:
            pos = left.start
            prev_layer, new_layer = new_layer, []
        elif left.start < pos:
            raise ValueError('unordered locus previous start={}, current start={}'.format(pos, left.start))

        for i in set(indices):
            if i:
                node = AltNode(left.start, left.stop, left.alleles[i], locus)

                # Add zygosity constraint
                vg.zygosity_constraints[node] = indices.count(i)

                # Add phase constraint for hets only
                if sample.phased and len(indices) == 2 and indices[0] != indices[1]:
                    phasename = 'M' if indices[0] == i else 'F'
                    vg.phase_constraints[phasename].add(node)
            else:
                node = ref_cache[left.start, left.stop]

            _extend_graph(vg, prev_layer, ref_cache, node)

            new_layer.append(node)

    _extend_graph(vg, new_layer, ref_cache, stop_node)

    return vg


def _extend_graph(vg, prev_layer, ref_cache, node):
    # Link into previous nodes
    for prev in prev_layer:
        # Do not link back to nodes that overlap current node
        if prev.stop > node.start:
            continue

        # Slice in a reference sequence if nodes do not join cleanly
        if prev.stop < node.start:
            ref_gap = ref_cache[prev.stop, node.start]
            vg.graph[prev].add(ref_gap)
            prev = ref_gap

        # Join new node to previous
        vg.graph[prev].add(node)


def _apply_phase_constrants(nodes, phasemap, phasesets, antiphasesets):
    if phasesets:
        # If any adjacent node belongs to an already selected phase set,
        # then allow only those nodes.  Otherwise, all nodes are valid.
        nodes = [node for node in nodes if phasemap.get(node) in phasesets] or nodes

    if antiphasesets:
        # Remove any node belonging to an anti-phaseset
        nodes = [node for node in nodes if phasemap.get(node) not in antiphasesets]

    return nodes


def generate_paths(vg, order='bfs'):
    if order not in ('bfs', 'dfs'):
        raise ValueError('invalid traversal order specified: {}'.format(order))

    graph, start = vg.graph, vg.start_node

    phasemap = {node: name for name, phaseset in vg.phase_constraints.iteritems() for node in phaseset}

    # queue of (path, pathset, phasesets, antiphasesets)
    path = [start]
    queue = deque([(path, set(path), set(), set())])

    next_node = queue.popleft if order == 'bfs' else queue.pop

    while queue:
        path, pathset, phasesets, antiphasesets = next_node()
        adjacent = [node for node in graph[path[-1]] if node not in pathset]

        # yield complete paths
        if not adjacent:
            yield path
            continue

        # Set of new phase sets being added from this node
        # nb: must occur prior to pruning
        add_phasesets = set(phasemap.get(node) for node in adjacent if node in phasemap and node not in phasesets)

        # prune adjacent nodes based on phase and anti-phase constraints
        adjacent = _apply_phase_constrants(adjacent, phasemap, phasesets, antiphasesets)

        for node in adjacent:
            new_pathset = pathset.copy()
            new_pathset.add(node)
            new_phasesets = _update_phasesets(phasesets, phasemap.get(node))
            new_antiphasesets = _update_antiphasesets(antiphasesets, add_phasesets, phasemap.get(node))
            queue.append((path + [node], new_pathset, new_phasesets, new_antiphasesets))


def _update_phasesets(phasesets, phaseset):
    if phaseset is not None and phaseset not in phasesets:
        # only copy when adding a new phaseset
        phasesets = phasesets.copy()
        phasesets.add(phaseset)
    return phasesets


def _update_antiphasesets(antiphasesets, add_phasesets, phaseset):
    if add_phasesets:
        # Copy if adding a new anti-phaseset, one not equal to the current node's phaseset
        add_anti = set(p for p in add_phasesets if p != phaseset and p not in antiphasesets)
        if add_anti:
            antiphasesets = antiphasesets | add_anti
    return antiphasesets


def generate_genotypes(ref, start, stop, loci, name, debug=False):
    vg = make_vargraph(ref, start, stop, loci, name)

    if debug:
        print('-' * 80)
        print('VG [{:d}, {:d})'.format(start, stop))
        for src in sorted(vg.graph):
            print('  {} -> {}'.format(src, vg.graph[src]))
        print()
        print('ZYGOSITY CONSTRAINTS:', vg.zygosity_constraints)
        print()
        print('PHASE CONSTRAINTS:', vg.phase_constraints)
        print()

    paths = map(tuple, generate_paths(vg))
    valid_paths = [path for path in paths if is_valid_path(vg, path)]

    # Any het constraint remove the need to consider homozygous genotypes
    # FIXME: diploid assumptions
    if 1 in vg.zygosity_constraints.values():
        pairs = combinations(valid_paths, 2)
    else:
        pairs = combinations_with_replacement(valid_paths, 2)

    valid_pairs = [pair for pair in pairs if is_valid_geno(vg, pair)]
    valid_genos = sorted(set(tuple(sorted([get_path_seq(p1), get_path_seq(p2)])) for p1, p2 in valid_pairs))

    if debug:
        print('PATHS:')
        for i, path in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, path))
        print()

        print('VALID PATHS:')
        for i, path in enumerate(valid_paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, path))
        print()

        print('VALID HAPLOTYPES:')
        for i, path in enumerate(valid_paths, 1):
            assert len(path) == len(set(path))
            print('{:4d}: {}'.format(i, get_path_seq(path)))
        print()

        print('VALID GENOTYPES:')

        for i, (allele1, allele2) in enumerate(valid_genos, 1):
            print('{:4d}: {}/{}'.format(i, allele1, allele2))
        print()

    return valid_genos

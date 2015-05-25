from itertools        import combinations, combinations_with_replacement
from collections      import deque, namedtuple, defaultdict, Counter


RefNode = namedtuple('RefNode', 'start stop seq')
AltNode = namedtuple('AltNode', 'start stop seq locus')
VariantGraph = namedtuple('VariantGraph', 'graph start_node zygosity_constraints phase_constraints antiphase_constraints')


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
    for a1, a2 in vg.antiphase_constraints:
        if not ((a1 in path) ^ (a2 in path)):
            return False

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
    vg = VariantGraph(defaultdict(set), start_node, defaultdict(int), defaultdict(set), [])
    ref_cache = RefCache(ref)

    prev_layer = [start_node]

    for locus in loci:
        sample = locus.record.samples[name]
        indices = sample.allele_indices
        left = locus.left

        new_layer = []
        for i in set(indices):
            if i:
                node = AltNode(left.start, left.stop, left.alleles[i], locus)

                # Add zygosity constraint
                vg.zygosity_constraints[node] = indices.count(i)

                # Add phase constraint for hets only
                if sample.phased and len(indices) == 2 and indices[0] != indices[1]:
                    phaseset = 'M' if indices[0] == i else 'F'
                    vg.phase_constraints[phaseset].add(node)
            else:
                node = ref_cache[left.start, left.stop]

            _extend_graph(vg, prev_layer, ref_cache, node)

            new_layer.append(node)

        # Add antiphase (het) constraint for alt/alt
        if 0 not in indices and len(new_layer) == 2:
            vg.antiphase_constraints.append(tuple(new_layer))

        prev_layer = new_layer

    _extend_graph(vg, prev_layer, ref_cache, stop_node)

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


def _apply_phase_constrants(nodes, phasemap, phasesets):
    if phasesets:
        for node in nodes:
            if phasemap.get(node) in phasesets:
                return [node]
    return nodes


def generate_paths(vg):
    graph, start = vg.graph, vg.start_node

    phasemap = {node: name for name, phaseset in vg.phase_constraints.iteritems() for node in phaseset}

    # queue of (path, pathset, phasesets)
    path = [start]
    queue = deque([(path, set(path), set())])
    while queue:
        path, pathset, phasesets = queue.popleft()
        adjacent = [node for node in graph[path[-1]] if node not in pathset]

        # yield at end of path
        if not adjacent:
            yield path
            continue

        # otherwise process adjacent nodes subjects to phase constraints
        for node in _apply_phase_constrants(adjacent, phasemap, phasesets):
            new_pathset = pathset.copy()
            new_pathset.add(node)
            new_phasesets = _update_phasesets(phasesets, phasemap.get(node))
            queue.append((path + [node], new_pathset, new_phasesets))


def _update_phasesets(phasesets, phaseset):
    if phaseset is not None and phaseset not in phasesets:
        phasesets = phasesets.copy()
        phasesets.add(phaseset)
    return phasesets


def generate_genotypes(ref, start, stop, loci, name):
    vg = make_vargraph(ref, start, stop, loci, name)

    if 0:
        print '-' * 80
        print 'VG [{:d}, {:d})'.format(start, stop)
        for src in sorted(vg.graph):
            print '  {} -> {}'.format(src, vg.graph[src])
        print
        print 'ZYGOSITY CONSTRAINTS:', vg.zygosity_constraints
        print
        print 'ANTI-PHASE CONSTRAINTS:', vg.antiphase_constraints
        print
        print 'PHASE CONSTRAINTS:', vg.phase_constraints
        print

    paths = map(tuple, generate_paths(vg))

    valid_paths = [path for path in paths if is_valid_path(vg, path)]

    # Any het constraint precludes looking at homozygous genotypes
    if 1 in vg.zygosity_constraints.values():
        comb = combinations(valid_paths, 2)
    else:
        comb = combinations_with_replacement(valid_paths, 2)

    valid_pairs = [pair for pair in comb if is_valid_geno(vg, pair)]
    valid_genos = sorted(set(tuple(sorted([get_path_seq(p1), get_path_seq(p2)])) for p1, p2 in valid_pairs))

    if 0:
        print 'PATHS:'
        for i, path in enumerate(paths, 1):
            assert len(path) == len(set(path))
            print '{:4d}: {}'.format(i, path)
        print

        print 'VALID PATHS:'
        for i, path in enumerate(valid_paths, 1):
            assert len(path) == len(set(path))
            print '{:4d}: {}'.format(i, path)
        print

        print 'VALID HAPLOTYPES:'
        for i, path in enumerate(valid_paths, 1):
            assert len(path) == len(set(path))
            print '{:4d}: {}'.format(i, get_path_seq(path))
        print

        print 'VALID GENOTYPES:'

        for i, (allele1, allele2) in enumerate(valid_genos, 1):
            print '{:4d}: {}/{}'.format(i, allele1, allele2)
        print

    return valid_genos

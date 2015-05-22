from itertools        import imap, combinations, combinations_with_replacement
from collections      import namedtuple, defaultdict, Counter

from vgraph.graph     import bfs_paths
from vgraph.iterstuff import unique_everseen


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


def is_valid_path(path, vg):
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


def is_valid_geno(paths, vg):
    observed_zygosity = Counter(p for path in paths for p in path if isinstance(p, AltNode))
    return vg.zygosity_constraints == observed_zygosity


def get_path_seq(path):
    return ''.join(p.seq for p in path)


def make_vargraph(ref, start, stop, loci, name):
    start_node = RefNode(start, start, '')
    vg = VariantGraph(defaultdict(set), start_node, defaultdict(int), defaultdict(set), [])
    ref_cache = RefCache(ref)

    last_layer = [start_node]

    for locus in loci:
        left_locus = locus.left
        sample = locus.record.samples[name]
        indices = sample.allele_indices

        new_layer = []
        for i in set(indices):
            if i:
                node = AltNode(left_locus.start, left_locus.stop, left_locus.alleles[i], locus)

                # Add zygosity constraint
                vg.zygosity_constraints[node] = indices.count(i)

                # Add phase constraint for hets only
                if sample.phased and len(indices) == 2 and indices[0] != indices[1]:
                    phaseset = 'M' if indices[0] == i else 'F'
                    vg.phase_constraints[phaseset].add(node)
            else:
                node = ref_cache[left_locus.start, left_locus.stop]

            # Link into previous nodes
            for last in last_layer:
                # Do not link back to nodes that overlap current node
                if last.stop > node.start:
                    continue

                # Slice in a reference sequence if nodes do not join cleanly
                if last.stop < node.start:
                    ref_gap = ref_cache[last.stop, node.start]
                    vg.graph[last].add(ref_gap)
                    last = ref_gap

                # Join new node to previous
                vg.graph[last].add(node)

            new_layer.append(node)

        # Add antiphase (het) constraint for alt/alt
        if 0 not in indices and len(new_layer) == 2:
            vg.antiphase_constraints.append(tuple(new_layer))

        last_layer = new_layer

    for last in last_layer:
        vg.graph[last].add(ref_cache[last.stop, stop + 1])

    return vg


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


    paths = imap(tuple, bfs_paths(vg.graph, vg.start_node))
    paths = list(unique_everseen(paths))

    valid_paths = [path for path in paths if is_valid_path(path, vg)]

    # Any het constraint precludes looking at homozygous genotypes
    if 1 in vg.zygosity_constraints.values():
        comb = combinations(valid_paths, 2)
    else:
        comb = combinations_with_replacement(valid_paths, 2)

    valid_pairs = [pair for pair in comb if is_valid_geno(pair, vg)]
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

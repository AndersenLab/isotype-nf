#!/usr/bin/env python

import sys
import copy

import numpy
import scipy.cluster.hierarchy

w_cutoff = 0.99986

def main():
    gt_fname, previso_fname, out_fname = sys.argv[1:4]

    weights, names = load_gtcheck(gt_fname)
    Z = cluster_discordance(weights)

    old_groups, refstrains = load_prev_isotypes(previso_fname, names)

    N = weights.shape[0]
    start = N - len(old_groups) - 120
    end = N - len(old_groups) + 30

    compare_clustering(Z, old_groups, start, end, out_fname)

def load_gtcheck(fname):
    column = {"discordance": 0, "sites": 0, "i": 0, "j": 0}
    dtypes = {"discordance": numpy.int32, "sites": numpy.int32, "i": "U8", "j": "U8"}
    names = {"discordance": "discordant", "sites": "sites", "i": "strain1", "j": "strain2"}
    fs = open(fname)
    header = fs.readline().rstrip('\n').split("\t")
    for i, name in enumerate(header):
        if name in column:
            column[name] = i
    dtype = []
    cols = []
    for k, v in column.items():
        cols.append(v)
        name = names[k]
        dtype.append((name, dtypes[k]))
    data = numpy.loadtxt(fname, skiprows=1, dtype=numpy.dtype(dtype), usecols=cols, delimiter="\t")

    # Convert discordance to concordance
    disconcordance = (data['discordant']) / data['sites']

    # Create concordance matrix
    names = numpy.unique(numpy.r_[data['strain1'], data['strain2']])
    strain_indices = {str(name): i for i, name in enumerate(names)}
    N = names.shape[0]
    weights = numpy.zeros((N, N), numpy.float32)
    for i in range(data.shape[0]):
        v1 = strain_indices[str(data['strain1'][i])]
        v2 = strain_indices[str(data['strain2'][i])]
        weights[v1, v2] = disconcordance[i]
        weights[v2, v1] = disconcordance[i]
    return weights, names

def load_prev_isotypes(fname, names):
    name_set = set([name for name in names])
    name_indices = {name: i for i, name in enumerate(names)}
    cols = {"strain": 0, "isotype": 0, "isotype_ref_strain": 0}
    fs = open(fname)
    header = fs.readline().rstrip().split("\t")
    for i, name in enumerate(header):
        if name in cols:
            cols[name] = i
    groups = {}
    refstrains = {}
    for line in fs:
        line = line.rstrip().split("\t")
        strain = line[cols['strain']]
        isotype = line[cols['isotype']]
        refstrain = line[cols['isotype_ref_strain']]
        if isotype == "NA":
            continue
        if strain not in name_set:
            continue
        refstrains[isotype] = refstrain
        groups.setdefault(isotype, [])
        groups[isotype].append(name_indices[strain])
    fs.close()
    return groups, refstrains

def cluster_discordance(weights):
    # cluster with complete linkage
    N = weights.shape[0]
    Z = scipy.cluster.hierarchy.linkage(weights[numpy.triu_indices(N, 1)], method='complete')

    # rescore tree using min between discordance - max within discordance
    W = numpy.zeros((N, N, 2), numpy.float32)
    W[:, :, :] = weights[:, :, numpy.newaxis]
    W[numpy.arange(N), numpy.arange(N), 0] = numpy.inf
    W[numpy.arange(N), numpy.arange(N), 1] = -numpy.inf
    groups = {i: [i] for i in range(N)}
    mapto = {i: i for i in range(N)}
    diffs = numpy.nanmin(weights, axis=0)
    for i in range(N - 2):
        x = int(Z[i, 0])
        y = int(Z[i, 1])
        z = mapto[x]
        w = mapto[y]
        W[z, :, 0] = numpy.nanmin(W[(w, z), :, 0], axis=0)
        W[:, z, 0] = numpy.nanmin(W[:, (w, z), 0], axis=1)
        W[z, :, 1] = numpy.nanmax(W[(w, z), :, 1], axis=0)
        W[:, z, 1] = numpy.nanmax(W[:, (w, z), 1], axis=1)
        groups[i + N] = groups[x] + groups[y]
        mapto[i + N] = z
        del groups[x]
        del groups[y]
        del mapto[x]
        del mapto[y]
        valid = list(mapto.values())
        valid.sort()
        valid = numpy.array(valid)
        without = list(set(valid).difference(set([z])))
        diffs[z] = numpy.nanmin(W[z, without, 0]) - W[z, z, 1]
        Z[i, 2] = numpy.amin(diffs[valid])
        Z[i, 3] = 1 - numpy.amax(W[valid, valid, 1])
    output = open("clustering.txt", "w")
    for i in range(Z.shape[0]):
        output.write(f"{Z[i, 0]}\t{Z[i, 1]}\t{Z[i, 2]:0.6f}\t{Z[i, 3]: 0.6f}\n")
    output.close()
    return Z

def compare_clustering(Z, old_groups, start, end, fname):
    N = Z.shape[0] + 1
    new_groups = {i: [i] for i in range(N)}
    for i in range(start):
        x = int(Z[i, 0])
        y = int(Z[i, 1])
        new_groups[N + i] = new_groups[x] + new_groups[y]
        del new_groups[x]
        del new_groups[y]
    output = open(fname, 'w')
    output.write("Score\tCutoff\tNewGroups\tOldGroups\tNewIsotypes\tIdentical\tNewJoins\tNewSplits\n")
    for i in range(start, end):
        stats = compare_groups(old_groups, list(new_groups.values()))
        tmp = "\t".join([f"{stats[name]}" for name in ['new_groups', 'old_groups', 'new', 'identical', 'join', 'split']])
        output.write(f"{Z[i, 2]:0.5e}\t{Z[i, 3]:0.7f}\t{tmp}\n")
        x = int(Z[i, 0])
        y = int(Z[i, 1])
        new_groups[N + i] = new_groups[x] + new_groups[y]
        del new_groups[x]
        del new_groups[y]
    output.close()

def compare_groups(old_groups, new_groups):
    stats = {"new_groups": len(new_groups), "old_groups": len(old_groups), 'new': 0, 'identical': 0, 'split': 0, 'join': 0}
    # Keep track of which isotype group each strain belongs to
    old_strain_iso = {}
    for i, group in old_groups.items():
        old_strain_iso.update({j: i for j in group})
    new_strain_iso = {}
    for i, group in enumerate(new_groups):
        new_strain_iso.update({j: i for j in group})
    used = set()
    for i, group in enumerate(new_groups):
        # If group has already been handled, skip it
        if i in used:
            continue
        new_strains = set()
        novel_strains = set()
        # Parse new strains into those present in old strain set and novel strains
        for strain in group:
            if strain in old_strain_iso:
                new_strains.add(strain)
            else:
                novel_strains.add(strain)
        new_isos = set([i])
        old_strains = set()
        old_isos = set()
        # Iterate until sets of old and new groups contain the same strains
        while len(new_strains.intersection(old_strains)) != len(new_strains):
            for strain in new_strains:
                old_iso = old_strain_iso[strain]
                old_isos.add(old_iso)
                old_strains = old_strains.union(set(list(old_groups[old_iso])))
            for strain in old_strains:
                new_iso = new_strain_iso[strain]
                new_isos.add(new_iso)
                for new_strain in new_groups[new_iso]:
                    if new_strain in old_strain_iso:
                        new_strains.add(new_strain)
                    else:
                        novel_strains.add(new_strain)
        # If there is a single new group and one or zero old groups, no need to worry about compound groups
        if len(new_isos) == 1 and len(old_isos) <= 1:
            # If no old groups matched, new group consists of novel strains and reference will be picked based on coverage
            if len(old_isos) == 0:
                stats['new'] += 1
            # If old group matches, use its isotype reference strain
            else:
                new_iso = list(old_isos)[0]
                stats['identical'] += 1
        # If there are more than one new or old groups, there has been a join/split and it is a compound group
        else:
            # For each new group, check which old groups overlap and if the overlap contains the isotype reference strain
            stats['split'] += len(new_isos) - 1
            stats['join'] += len(old_isos) - 1
        # Don't reexamine already parsed new groups
        for new_iso_num in new_isos:
            used.add(new_iso_num)
    return stats

main()
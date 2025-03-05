#!/usr/bin/env python

import sys
import copy

import numpy
import scipy.cluster.hierarchy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.switch_backend('cairo')


def main():
    gt_fname, previso_fname, coverage_fname, cutoff = sys.argv[1:5]
    cutoff = float(cutoff)

    weights, names = load_gtcheck(gt_fname)

    Z = cluster_discordance(weights)
    new_groups = cut_tree(Z, cutoff, names)

    old_groups, refstrains = load_prev_isotypes(previso_fname, names)
    old_names = []
    for strains in old_groups.values():
        old_names += list(strains)
    coverages = load_coverages(coverage_fname)

    new_groups, compound_groups, stats = compare_groups(old_groups, new_groups, coverages)

    write_new_groups(new_groups, refstrains)
    # write_wi_isotype_sample_sheet(new_groups, refstrains)
    write_summary(old_groups, new_groups, weights.shape[0], stats)
    plot_isotype_comparison(compound_groups, names, weights, cutoff)

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

def cut_tree(Z, cutoff, names):
    # Recapitulate clustering until scoring cutoff is reached
    N = Z.shape[0] + 1
    groups = {i: [i] for i in range(N)}
    for i in range(Z.shape[0]):
        if Z[i, 3] < cutoff:
            break
        x = int(Z[i, 0])
        y = int(Z[i, 1])
        groups[N + i] = groups[x] + groups[y]
        del groups[x]
        del groups[y]
    new_groups = []
    for group in groups.values():
        new_groups.append([names[i] for i in group])
    return new_groups

def load_prev_isotypes(fname, names):
    name_set = set([name for name in names])
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
        groups[isotype].append(strain)
    fs.close()
    return groups, refstrains

def load_coverages(fname):
    coverages = {line.split("\t")[0]: float(line.rstrip().split("\t")[1]) for line in open(fname)}
    return coverages

def compare_groups(old_groups, new_groups, coverages):
    stats = {'new': 0, 'identical': 0, 'split': 0, 'join': 0}
    # Keep track of which isotype group each strain belongs to
    old_strain_iso = {}
    for i, group in old_groups.items():
        old_strain_iso.update({j: i for j in group})
    new_strain_iso = {}
    for i, group in enumerate(new_groups):
        new_strain_iso.update({j: i for j in group})
    used = set()
    new_isotypes = {}
    compound_isos = []
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
            new_group = list(new_strains.union(novel_strains))
            # If no old groups matched, new group consists of novel strains and reference will be picked based on coverage
            if len(old_isos) == 0:
                iso_coverages = [(coverages[strain], strain) for strain in new_group]
                iso_coverages.sort()
                new_iso = iso_coverages[-1][1]
                stats['new'] += 1
            # If old group matches, use its isotype reference strain
            else:
                new_iso = list(old_isos)[0]
                stats['identical'] += 1
            new_isotypes[new_iso] = new_group
        # If there are more than one new or old groups, there has been a join/split and it is a compound group
        else:
            new_group_isos = {}
            # For each new group, check which old groups overlap and if the overlap contains the isotype reference strain
            for new_iso_num in new_isos:
                new_group = set(new_groups[new_iso_num])
                iso_coverages = []
                for old_iso in old_isos:
                    overlap = new_group.intersection(set(old_groups[old_iso]))
                    if old_iso in overlap:
                        # If old isotype reference is in new group, need to keep track of coverages in case more than one old references present
                        iso_coverages.append((coverages[old_iso], old_iso))
                if len(iso_coverages) == 0:
                    # No old isotype references present, need to pick biggest coverage strain in new group
                    iso_coverages = [(coverages[strain], strain) for strain in new_group]
                iso_coverages.sort()
                new_iso = iso_coverages[-1][1]
                new_isotypes[new_iso] = new_group
                new_group_isos.update({strain: new_iso for strain in new_group})
            stats['split'] += len(new_isos) - 1
            stats['join'] += len(old_isos) - 1

            # Need to create compound isotype representations
            # First, determine minimum sets of overlap
            nodes = {}
            new_nodes = {}
            old_nodes = {}
            for n_iso in set(list(new_group_isos.values())):
                for o_iso in old_isos:
                    overlap = set(new_isotypes[n_iso]).intersection(set(old_groups[o_iso]))
                    if len(overlap) > 0:
                        nodes[len(nodes)] = {"strains": list(overlap), "edges": set([]), "new_iso": n_iso, "old_iso": o_iso}
                        nodes[len(nodes) - 1]['strains'].sort()
                        new_nodes.setdefault(n_iso, [])
                        old_nodes.setdefault(o_iso, [])
                        new_nodes[n_iso].append(len(nodes) - 1)
                        old_nodes[o_iso].append(len(nodes) - 1)
            for connections in new_nodes.values():
                for x in connections:
                    for y in connections:
                        if x == y:
                            continue
                        nodes[x]['edges'].add(y)
            for connections in old_nodes.values():
                for x in connections:
                    for y in connections:
                        if x == y:
                            continue
                        nodes[x]['edges'].add(y)
            for node in nodes.values():
                node['edges'] = list(node['edges'])
            # Identify possible starting points
            leaves = []
            for key, node in nodes.items():
                if len(set([nodes[edge]['new_iso'] for edge in node['edges']])) == 1 and len(set([nodes[edge]['old_iso'] for edge in node['edges']])) == 1:
                    leaves.append(key)
            best_path = []
            l_index = 0

            def find_path(nodes, node, path=[]):
                curr_path = list(path) + [node]
                best_path = []
                for edge in nodes[node]['edges']:
                    if edge in path:
                        continue
                    tmp_path = find_path(nodes, edge, curr_path)
                    if len(tmp_path) > len(best_path):
                        best_path = tmp_path
                return best_path
            
            # Find complete path through graph or closest to complete
            while len(best_path) < len(nodes) and l_index < len(leaves):
                path = find_path(nodes, leaves[l_index])
                if len(path) > len(best_path):
                    best_path = path
                l_index += 1
            for n in range(len(nodes)):
                if n not in best_path:
                    best_path.append(n)
            new_compound = []
            old_compound = []
            for n in best_path:
                if len(new_compound) == 0 or new_compound[-1][0] != nodes[n]['new_iso']:
                    if len(new_compound) > 0:
                        new_strains = set(new_isotypes[new_compound[-1][0]])
                        for o_iso in old_isos:
                            new_strains = new_strains.difference(old_groups[o_iso])
                        if len(new_strains) > 0:
                            new_compound[-1][1] += list(new_strains)
                            old_compound.append(["NA", list(new_strains)])
                    new_compound.append([nodes[n]['new_iso'], list(nodes[n]['strains'])])
                else:
                    new_compound[-1][1] += nodes[n]['strains']
                if len(old_compound) == 0 or old_compound[-1][0] != nodes[n]['old_iso']:
                    old_compound.append([nodes[n]['old_iso'], list(nodes[n]['strains'])])
                else:
                    old_compound[-1][1] += nodes[n]['strains']
            new_strains = set(new_isotypes[new_compound[-1][0]])
            for o_iso in old_isos:
                new_strains = new_strains.difference(old_groups[o_iso])
            if len(new_strains) > 0:
                new_compound[-1][1] += list(new_strains)
                old_compound.append(["NA", list(new_strains)])
            compound_isos.append([old_compound, new_compound])
        # Don't reexamine already parsed new groups
        for new_iso_num in new_isos:
            used.add(new_iso_num)
    return new_isotypes, compound_isos, stats


def write_new_groups(new_groups, refstrains):
    groups = []
    for key, values in new_groups.items():
        values = list(values)
        values.sort()
        groups.append((values, key))
    groups.sort()
    output = open("isotype_groups.tsv", "w")
    output.write("group\tstrain\tisotype\tisotype_ref_strain\n")
    for i, (group, isotype) in enumerate(groups):
        if isotype in refstrains:
            ref = refstrains[isotype]
        else:
            ref = isotype
        for name in group:
            output.write(f"{i}\t{name}\t{isotype}\t{ref}\n")
    output.close()

def write_wi_isotype_sample_sheet(new_groups, refstrains):
    output = open("wi_isotype_sample_sheet.txt", "w")
    for name in new_groups.keys():
        if name in refstrains:
            ref = refstrains[name]
        else:
            ref = name
        output.write(f"{ref}\n")
    output.close()

def write_summary(old_groups, new_groups, num_strains, stats):
    output = open('isotype_summary.txt', 'w')
    output.write(f"Old isotype groups\t{len(old_groups)}\n")
    output.write(f"New isotype groups\t{len(new_groups)}\n")
    output.write(f"Number of strains\t{num_strains}\n")
    for key, value in stats.items():
        output.write(f"{key}\t{value}\n")
    output.close()

def plot_isotype_comparison(compound_groups, names, weights, cutoff):
    name2index = {name: i for i, name in enumerate(names)}
    with PdfPages(f"isotype_comparison.pdf") as pdf:
        for i in range(len(compound_groups)):
            old_isogroups = compound_groups[i][0]
            new_isogroups = compound_groups[i][1]
            strains = []
            new_isos = []
            old_isos = []
            old_isogroups_na = []
            for x in new_isogroups:
                strains += x[1]
                new_isos += [x[0]] * len(x[1])
            prev_iso = ""
            for x in old_isogroups:
                old_isos += [x[0]] * len(x[1])
                if x[0] == "NA" or x[0] == prev_iso:
                    old_isogroups_na[-1] += x[1]
                else:
                    old_isogroups_na.append(x[1])
                    prev_iso = x[0]
            indices = [name2index[name] for name in strains]
            N = len(indices)
            fig, ax = plt.subplots(1, 2, figsize=(N / 2 + 5.5, N / 4 + 3.0))
            indices = numpy.array(indices)
            W = weights[indices, :][:, indices]
            W[numpy.arange(indices.shape[0]), numpy.arange(indices.shape[0])] = 0
            ax[0].imshow(1-W, vmin=cutoff, vmax=1)
            ax[0].set_xticks(numpy.arange(N))
            ax[0].set_xticklabels([f"{strains[j]} ({old_isos[j]})" for j in range(N)], rotation=35, ha='right')
            ax[0].set_yticks(numpy.arange(N))
            ax[0].set_yticklabels([f"{strains[j]} ({old_isos[j]})" for j in range(N)])
            wo_min = 0
            wo_max = 1
            wi = 0
            for x in old_isogroups_na:
                within = numpy.array([name2index[name] for name in x])
                wi = max(wi, numpy.amax(weights[within, :][:, within]))
                if len(old_isogroups_na) > 1:
                    without = numpy.array(list(set(indices).difference(set(within))))
                    wo_max = min(wo_max, numpy.amin(weights[within, :][:, without]))
                    wo_min = max(wo_min, numpy.amax(weights[within, :][:, without]))
            if len(old_isogroups_na) > 1:
                wo = f"{1 - wo_min:0.6f}-{1 - wo_max:0.6f}"
            else:
                wo = "NA"
            wi = f"{1 - wi:0.6f}"
            ax[0].set_title(f"old groups\nmin within:{wi}\nbetween:{wo}")
            ax[1].imshow(1-W, vmin=cutoff, vmax=1)
            ax[1].set_xticks(numpy.arange(N))
            ax[1].set_xticklabels([f"{strains[j]} ({new_isos[j]})" for j in range(N)], rotation=35, ha='right')
            ax[1].set_yticks(numpy.arange(N))
            ax[1].set_yticklabels([f"{strains[j]} ({new_isos[j]})" for j in range(N)])
            wo_min = 0
            wo_max = 1
            wi = 0
            for x in new_isogroups:
                within = numpy.array([name2index[name] for name in x[1]])
                wi = max(wi, numpy.amax(weights[within, :][:, within]))
                if len(new_isogroups) > 1:
                    without = numpy.array(list(set(indices).difference(set(within))))
                    wo_max = min(wo_max, numpy.amin(weights[within, :][:, without]))
                    wo_min = max(wo_min, numpy.amax(weights[within, :][:, without]))
            if len(new_isogroups) > 1:
                wo = f"{1 - wo_min:0.6f}-{1 - wo_max:0.6f}"
            else:
                wo = "NA"
            wi = f"{1 - wi:0.6f}"
            ax[1].set_title(f"new groups\nmin within:{wi}\nbetween:{wo}")
            prev_old = old_isos[0]
            for j in range(1, N):
                if old_isos[j] != "NA" and old_isos[j] != prev_old:
                    ax[0].vlines([j - 0.5], ymin=-0.5, ymax=N - 0.5, color='white', lw=5)
                    ax[0].hlines([j - 0.5], xmin=-0.5, xmax=N - 0.5, color='white', lw=5)
                    ax[0].vlines([j - 0.5], ymin=-0.5, ymax=N - 0.5, color='black', lw=2)
                    ax[0].hlines([j - 0.5], xmin=-0.5, xmax=N - 0.5, color='black', lw=2)
                    prev_old = old_isos[j]
                if new_isos[j] != new_isos[j - 1]:
                    ax[1].vlines([j - 0.5], ymin=-0.5, ymax=N - 0.5, color='white', lw=5)
                    ax[1].hlines([j - 0.5], xmin=-0.5, xmax=N - 0.5, color='white', lw=5)
                    ax[1].vlines([j - 0.5], ymin=-0.5, ymax=N - 0.5, color='black', lw=2)
                    ax[1].hlines([j - 0.5], xmin=-0.5, xmax=N - 0.5, color='black', lw=2)
            plt.tight_layout()
            pdf.savefig()
            plt.close()


main()
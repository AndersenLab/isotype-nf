#!/usr/bin/env python

import sys
import multiprocessing
import multiprocessing.managers
from multiprocessing.shared_memory import SharedMemory

import numpy

def main():
    vcf_fname, out_fname, nthreads = sys.argv[1:4]
    nthreads = int(nthreads)
    tmp = numpy.load(vcf_fname, mmap_mode='r')
    strains = tmp.dtype.names[2:]
    triu = numpy.triu_indices(len(strains), 1)
    with multiprocessing.managers.SharedMemoryManager() as smm:
        pool = multiprocessing.Pool(nthreads)
        genotypes_view = smm.SharedMemory(len(strains) * tmp.shape[0] * numpy.dtype(numpy.int8).itemsize)
        genotypes_name = genotypes_view.name
        genotypes = numpy.ndarray((len(strains), tmp.shape[0]), dtype=numpy.int8, buffer=genotypes_view.buf)
        for i, strain in enumerate(strains):
            genotypes[i, :] = tmp[strain]
        results_view = smm.SharedMemory(triu[0].shape[0] * 2 * numpy.dtype(numpy.int32).itemsize)
        results_name = results_view.name
        results = numpy.ndarray((triu[0].shape[0], 2), numpy.int32, buffer=results_view.buf)
        indices = numpy.round(numpy.linspace(0, triu[0].shape[0], nthreads +  1)).astype(numpy.int32)
        futures = []
        for i in range(indices.shape[0] - 1):
            s, e = indices[i:i+2]
            kwargs = {'start': s, 'end': e, 'index1': triu[0][s:e], 'index2': triu[1][s:e], 'geno_shape': genotypes.shape,
                      'geno_name': genotypes_name, 'results_shape': results.shape, 'results_name': results_name}
            futures.append(pool.apply_async(find_discordance, kwds=kwargs))
        for result in futures:
            _ = result.get()
        pool.close()
        pool.terminate()        

        output = open(out_fname, 'w')
        output.write("discordance\tsites\ti\tj\n")
        for i in range(results.shape[0]):
            output.write(f"{results[i, 0]}\t{results[i, 1]}\t{strains[triu[0][i]]}\t{strains[triu[1][i]]}\n")
        output.close()

def find_discordance(start, end, index1, index2, geno_shape, geno_name, results_shape, results_name):
    genotypes_view = SharedMemory(geno_name)
    genotypes = numpy.ndarray(geno_shape, numpy.int8, buffer=genotypes_view.buf)
    results_view = SharedMemory(results_name)
    results = numpy.ndarray(results_shape, numpy.int32, buffer=results_view.buf)
    for i in range(end - start):
        x = index1[i]
        y = index2[i]
        valid = numpy.where(numpy.logical_and(genotypes[x, :] >= 0, genotypes[y, :] >= 0))[0]
        results[i + start, 1] = valid.shape[0]
        results[i + start, 0] = numpy.sum(genotypes[x, valid] != genotypes[y, valid])
    genotypes_view.close()
    results_view.close()
    return

main()

#!/usr/bin/env python

import sys
import gzip

import numpy

def main():
    vcf_fname, out_fname = sys.argv[1:3]
    fs = gzip.open(vcf_fname, 'rb')
    line = fs.readline()
    while not line.startswith(b"#C"):
        line = fs.readline()
    header = line.decode("utf-8").rstrip().split("\t")
    strain_index = header.index("FORMAT") + 1
    strains = header[strain_index:]
    data = []
    chrsize = 0
    snps = []
    for line in fs:
        line = line.rstrip().split(b"\t")
        if len(line[3]) > 1 or len(line[4]) > 1:
            snps.append(False)
        else:
            snps.append(True)
        data.append([line[0].decode('utf8'), int(line[1])])
        chrsize = max(chrsize, len(line[0]))
        for c in range(strain_index, len(line)):
            geno = line[c]
            # if geno.count(b"PASS") == 0:
            #     data[-1].append(-1)
            if geno.startswith(b"0/0"):
                data[-1].append("0")
            elif geno.startswith(b"1/1"):
                data[-1].append("1")
            else:
                data[-1].append("-1")
    dtype = [('chr', f"U{chrsize}"), ('pos', numpy.int32)] + [(name, numpy.int8) for name in strains]
    data = numpy.array([tuple(row) for row in data], numpy.dtype(dtype))
    snps = numpy.array(snps, bool)
    data = data[numpy.where(snps)]
    numpy.save(out_fname, data)

main()
#!/usr/bin/env python3
import sys, os, glob

def main():
    try:
        logfile_directory = sys.argv[1]
    except IndexError:
        print("Usage: <Directory of BLAST nr files>")
        return
    logfiles = glob.glob(logfile_directory+'/*')
    clusters = {}
    all_proteins = {}
    all_proteins['NA'] = 'NA'
    for lf in logfiles:
        logf = open(lf)
        in_section = False
        proteins = {}
        for line in logf:
            if line.strip():
                line = line.strip()
                if line.startswith('Sequences producing'):
                    in_section = True
                    continue
                elif in_section is True:
                    all_proteins[line.split(' ')[0]] = ' '.join(line.strip(' ').split(' ')[2:6])
                    try:
                        proteins[line.split(' ')[0]] += 1
                    except KeyError:
                        proteins[line.split(' ')[0]] = 1
                    in_section = False
                    continue
        cluster = lf.split('/')[-1]
        try:
            clusters[cluster] = max(proteins)
        except ValueError:
            clusters[cluster] = 'NA'
    print('cluster\tprotein_id\tprotein_info')
    for cluster, values in clusters.items():
        print('%s\t%s\t%s' % (cluster, values, all_proteins[values]))

if __name__ == '__main__':
    main()

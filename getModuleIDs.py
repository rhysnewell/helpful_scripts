#!/usr/bin/env python3
import sys, os, glob

def main():
    try:
        module_paths = sys.argv[1]
        group_paths = sys.argv[2]
    except IndexError:
        print("Usage: <Module Paths> <KO IDs> \n",
        "Module Paths: Module Paths tsv from enrichm classify \n",
        "KO IDs: Clusters, Orthologs, and KO IDs from getEnrichMKOs.py")
        return
    modules = open(module_paths)
    groups = open(group_paths)
    ko_dict = {}
    for line in modules:
        mod = line.split("\t")[1]
        KOs = line.split("\t")[-1].strip().split(',')  
        for ko in KOs:
            ko_dict[ko] = mod.strip()

        
    print('cluster\tortholog\tko\tmodule')
    for line in groups:
        line = line.split('\t')
        try:
            print('%s\t%s\t%s\t%s' % (line[0], line[1], line[2].strip(), ko_dict[line[2].strip()]))
        except KeyError:
            pass


if __name__ == '__main__':
	main()        
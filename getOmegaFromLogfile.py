#!/usr/bin/env python3
import sys, os, glob

def main():
    try:
        logfile_directory = sys.argv[1]
    except IndexError:
        print("Usage: <Directory of iqtree log files>")
        return
    logfiles = glob.glob(logfile_directory+'/*.log')
    rates = {}
    for lf in logfiles:
        logf = open(lf)
        omega = -1
        kappa = -1
        kappa2 = -1
        for line in logf:
            if line.startswith('Nonsynonymous'):
                omega = float(line.split(' ')[-1])
            elif line.startswith('Transition'):
                k = line.split(' ')
                if len(k) == 4:
                    kappa = float(k[-1])
                else:
                    kappa2 = float(k[-1])
        cluster = lf.split('/')[-1]
        rates[cluster] = [omega, kappa, kappa2]
    print('cluster\tomega\tkappa\tkappa2')
    for cluster, values in rates.items():
        print('%s\t%.3f\t%.3f\t%.3f' % (cluster.split('.')[0], values[0], values[1], values[2]))
        
if __name__ == '__main__':
	main()        
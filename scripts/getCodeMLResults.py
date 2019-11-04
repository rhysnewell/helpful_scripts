#!/usr/bin/env python3
import sys, os, glob

def main():
    try:
        logfile_directory = sys.argv[1]
    except IndexError:
        print("Usage: <Directory of CodeML log files>")
        return
    logfiles = glob.glob(logfile_directory+'/*')
    rates = {}
    for lf in logfiles:
        logf = open(lf)
        M8_omega = -1
        M7_omega = -1
        p_value = -1
        M8 = False
        M7 = False
        for line in logf:
            line = line.strip()
            if line.startswith('M7~'):
                p_value = float(line.split('|')[-1].strip('*'))
            elif line.startswith('- Model M8'):
                M8 = True
            elif line.startswith('- Model M7'):
                M7 = True
            elif line.startswith('*'):
                if M8:
                    M8_omega = float(line.split(': ')[-1])
                    M8 = False
                if M7:
                    M7_omega = float(line.split(': ')[-1])
                    M7 = False        
        cluster = lf.split('/')[-1]
        rates[cluster] = [M7_omega, M8_omega, p_value]
    print('cluster\tm7_null\tm8_alternative\tp_value')
    for cluster, values in rates.items():
        print('%s\t%.3f\t%.3f\t%.3f' % (cluster.split('.')[0], values[0], values[1], values[2]))
        
if __name__ == '__main__':
	main()        
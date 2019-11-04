#!/usr/bin/env python3
import sys, os, glob
from numpy import *

def main():
    try:
        logfile = sys.argv[1]
    except IndexError:
        print("Usage: <Benchmark log files with times>")
        return
        
    lf = open(logfile)
    runs = {}
    tool = None
    run = None
    data = None
    for line in lf:
        if line.startswith('###'):
            data = ''.join(line.strip().strip('#').strip(" ").split(':'))          
            data = '_'.join(data.split(' '))
            #run = tool + '_' + run
            
        elif line.count('COVERM')==1 or line.count('BAMM')==1 or line.count('METABAT')==1 or line.count('MINIMAP2'):
            
            tool = line.strip().strip('!!').strip('#').strip(' ')        
            tool = '_'.join(tool.split(' ')) 
            run = tool + '_' + data
            #print(tool)
            if run not in runs.keys():
                runs[run] = {}

        elif line.count('user')==1 and line.count('elapsed')==1:
            line = line.strip().split()
            user = float(line[0].split('user')[0])
            syst = float(line[1].split('system')[0])
            real = line[2].split('elapsed')[0]
            real = real.split(':')
            real_seconds = int(real[0])*60 + float(real[1])
            try:
                runs[run]['real'].append(real_seconds)
                runs[run]['user'].append(user)
                runs[run]['sys'].append(syst)
            except KeyError:
                runs[run]['real'] = [real_seconds]
                runs[run]['user'] = [user]
                runs[run]['sys']  = [syst]
                     
        elif line.startswith('real'):
            #print('real')
            time = line.strip().split('\t')
            seconds = float(time[1].split('m')[0])*60+float(time[1].split('m')[1].strip('s'))
            try:
                runs[run]['real'].append(seconds)
            except KeyError:
                runs[run]['real'] = [seconds]
                
        elif line.startswith('user'):
            time = line.strip().split('\t')
            seconds = float(time[1].split('m')[0])*60+float(time[1].split('m')[1].strip('s'))
            try:
                runs[run]['user'].append(seconds)
            except KeyError:
                runs[run]['user'] = [seconds]
                
        elif line.startswith('sys'):
            time = line.strip().split('\t')
            seconds = float(time[1].split('m')[0])*60+float(time[1].split('m')[1].strip('s'))
            try:
                runs[run]['sys'].append(seconds)
            except KeyError:
                runs[run]['sys'] = [seconds]
                

    #print(runs)  
    print('tool\treal_1\tuser_1\tsys_1\treal_2\tuser_2\tsys_2\treal_3\tuser_3\tsys_3\treal_mean\tuser_mean\tsys_mean')
    for tool, run in runs.items():
        #print(tool, run)
        print('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (tool,
        run['real'][0], run['user'][0], run['sys'][0],
        run['real'][1], run['user'][1], run['sys'][1],
        run['real'][2], run['user'][2], run['sys'][2],
        mean(run['real']), mean(run['user']), mean(run['sys'])))

if __name__ == '__main__':
    main()

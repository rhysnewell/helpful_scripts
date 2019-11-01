#!/usr/bin/env python3
import os, sys, csv
import kneed

def main(file):
    dev_strain_cnt = []
    dev_dev = []
    with open(file) as dev:
        dev_csv = csv.reader(dev)
        for idx, row in enumerate(dev_csv):
            if idx == 0:
                continue
            else:
                dev_strain_cnt.append(int(row[0]))
                dev_dev.append(float(row[3]))


    dev_strain_cnt, dev_dev = zip(*sorted(zip(dev_strain_cnt, dev_dev)))
    dev_dev = list(dev_dev)
    devs = []
    devs_mean = []
    cnt = 0
    for dev in dev_dev:
        cnt += 1
        devs.append(dev)
        if cnt == 4:
            devs_mean.append(sum(devs)/len(devs))
            cnt = 0
            devs = []

    x = range(len(devs_mean))
    kn = kneed.KneeLocator(x, devs_mean,
                           curve = 'convex',
                           direction = 'decreasing',
                           interp_method='polynomial')
    print(kn.knee+1)

#main('tests/dev.csv')
if __name__=="__main__":
    dev_file = sys.argv[1]
    main(dev_file)

#!/usr/bin/python3
import os
import sys
import json

def ReadFile(f):
    z = open(f, 'r')
    lines = z.readlines()
    return ''.join(lines)

class Payload(object):
    def __init__(self, f):
        j = ReadFile(f)
        self.__dict__ = json.loads(j)
    def save(self, f):
        with open(f, "w") as text_file:
            text_file.write(json.dumps(self.__dict__))

def main():
    results = []
    orders = [2, 4, 6, 8]
    numCell = [12, 24, 32, 48, 64]
    device = {1:'cpu', 2:'gpu'}
    optL = ['opt', 'dbg']
    envVar = 'HIDXOUTPUT'
    if (envVar not in os.environ):
        print('Please set environment variable \'{}\' before running'.format(envVar))
        sys.exit(199)
    outDir = os.environ[envVar]
    for op in optL:
        for dev in device:
            for num in numCell:
                for ord in orders:
                    outfile  = 'output/result.json'
                    savefile = '{}/run_O{}_N{}_D{}_{}.json'.format(outDir, ord, num, device[dev], op)
                    config = '-Dord={} -Dnum={} -Ddev={} -Doutfile={}'.format(ord, num, device[dev], outfile)
                    cmd = './program.{} {}'.format(op, config)
                    print('running: \"{}\"'.format(cmd))
                    sys.stdout.flush()
                    os.system(cmd)
                    pl = Payload(outfile)
                    pl.save(savefile)
                    results.append(pl)
    
    return 0

sys.exit(main())

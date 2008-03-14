#!/usr/bin/env python
"""
Parameter scans
Usage:
   - adapt variables 'params', 'pargs' and 'results'
   - sample run:
        python paramscan.py './evalfct -N 100 -n 1250' res_artdat_wgn.txt
"""
import sys, os

params = {'L':range(0,11),
          'e':['zero','sym', 'zeror', 'smooth'],
          'f':['sureshrink', 'heursure', 'ti', 'conventional'],
          't':['s', 'h']};
#params = {'h':[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],
#          'g':[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]};
pargs= sys.argv[2]; #'data/artdat_wgn_sigreal_100_samples1250.dat'

results = ['avg(RMSE)', 'avg(SNR)',
           'sd(RMSE)', 'sd(SNR)'];

if len(sys.argv)!=4:
    print "Usage: paramscan.py <prog> <input> <output>"
    sys.exit();
    
prog = sys.argv[1];
out  = sys.argv[3];

keys = params.keys();

f = open(out, 'w');

indices={};
i=1
for key in keys+results:
    indices[key]=i;
    f.write(key+"\t")
    i+=1;
f.write("\n");


def runprog(cmd, w, params, keys, rec):
    global indices, f, results
    
    for p in params[keys[rec]]:
        if rec==len(keys)-1:
            print cmd+' -'+str(keys[rec])+' '+str(p)+' '+str(pargs);
            d = os.popen( cmd+' -'+keys[rec]+' '+str(p)+' '+str(pargs), 'r');
            lines = d.readlines()
            res = []
            for line in lines:
                for r in results:
                    if (line.strip()).startswith(r):
                        tmp = line.split('=');
                        tmp = tmp[len(tmp)-1].strip();
                        res.append(tmp)
            for x in w+[p]+res:
                f.write(str(x)+"\t");
            f.write("\n");
        else:
            runprog(cmd+' -'+keys[rec]+' '+str(p),w+[p], params, keys, rec+1);



runprog(prog, [], params, keys, 0);
    
f.close();

#!/usr/bin/env python
import string
import time
import argparse
from scipy.optimize import curve_fit
import numpy as np
import logging
import random
import math

from doris.primer_util import *
from doris.design import *

from doris.base_conversion import convertQuarnary, convertTernary
from doris.huffman import huffman_encode_byte,rotate_encode

needed_distance = { 4: 1,
                    5: 1,
                    6: 1,
                    7: 2,
                    8: 2,
                    9: 2,
                    10: 2,
                    11: 4,
                    12: 4}


# design a comma free codec
def codec_design(l):
    all = []
    if l <= 5:    
        for i in range(4**l):
            d = convertQuarnary(i,l)
            #print d
            if hasRepeat(d) and l>6:
                continue
            bad_distance = False
            for a in all:
                dd = hamming_distance(a,d)
                if dd < needed_distance[l]:
                    bad_distance=True
                    break

            if bad_distance:
                continue

            all.append(d)

        random.shuffle(all)
        return all[0:256]
    else:
       for i in range(3**l):
            if len(all) > 256:
                break
            d = convertTernary(i,l)
            if d in all:
                continue
            bad_distance = False
            for a in all:
                dd = hamming_distance(a,d)
                if dd < needed_distance[l]:
                    bad_distance=True
                    break

            if bad_distance:
                continue

            all.append(d)

       random.shuffle(all)
       return all[0:256]
           
def montecarlo(trials,L,design_rules):
    pg = UniquePrimerGenerator(chars="AGCT",length=20)
    t = 0
    count = 0
    while (t < trials):
        s = pg.get()
        t += 1
        if design_rules.check(s,primers=L):
            L.append(s)
            pg.append(s)
            count = count+1
    return (t,count)

def get_experiment_header():
    return ['kind','codeword_size','pd','ld','lib_size','simulations','estimated_primers']

    
def get_experiment_info(args):
    info = {}
    info['kind'] = 'primer_density_study'
    info['ld'] = args.strand_distance
    info['pd'] = args.distance
    info['lib_size'] = args.libsize
    return info


if __name__ == "__main__":

    parser = argparse.ArgumentParser(\
        description="Estimate how primer count varies under different encodings.")

    parser.add_argument('--csv',type=str,dest="csv",action="store",\
                        default="", help="filename for csv file output")

    parser.add_argument('--samples',dest="samples", action="append", \
                        type=int, default=[500,1000,2000], \
                        help="# of trials to perform to create curve")

    parser.add_argument('--codeword-size',dest="codeword_size", action="append", \
                        type=int, default=[], \
                        help="width of codeword to use")

    parser.add_argument('--verbose',dest="verbose",default=False,action="store_true",\
                        help="Dump design rule information to stdout.")

    parser.add_argument('--distance',type=int,dest="distance",action="store",\
                        default=6, help="Hamming distance between primers")

    parser.add_argument('--strand-distance',type=int,dest="strand_distance",\
                        action="store",\
                        default=6, help="Hamming distance from primer to strands in the Library")

    parser.add_argument('--library-size',type=int,dest="libsize",action="store",\
                        default=10**5, help="Size of library to compare against")

    args = parser.parse_args()

    if len(args.csv)>0:
        csv_file = open(args.csv,"w")
        h = get_experiment_header()
        csv_file.write(",".join(h)+'\n')
        csv_file.close()
    
    codec_dict = {}
    libsize = args.libsize
    
    Raw_lib = []
    for i in range(0,libsize,20):
        r = [ random.randint(0,255) for _ in range(20) ]
        Raw_lib.append(r)

    #print lib

    for cw in args.codeword_size:
        c = codec_design(cw)
        codec_dict[cw] = c
        print (len(c))
        assert len(c) == 256

        Library = []
        for s in Raw_lib:
            if cw == 6:
                cs = [ huffman_encode_byte(_) for _ in s ]
                strand = "A"*20+ rotate_encode("".join(cs)) + "T"*20
            else:
                cs = [ c[_] for _ in s ]
                strand = "A"*20+ "".join(cs) + "T"*20
            Library.append(strand)

        design_rules = DesignRules("Codeword={} ".format(cw))

        design_rules.add_rule(build_GC_rule(45,55))
        design_rules.add_rule(build_longrun_rule())
        design_rules.add_rule(build_hamming_distance_library_rule(6,[]))

        design_rules.add_rule(build_self_complementarity_rule(4))
        design_rules.add_rule(build_inter_complementarity_rule(10,[]))

        design_rules.add_rule(build_Tm_rule(50,55))
        design_rules.add_rule(build_primer3_hairpin_rule())
        design_rules.add_rule(build_primer3_homodimer_rule())
        design_rules.add_rule(build_primer3_heterodimer_bindings_library_rule([]))
        design_rules.add_rule(build_primer3_check_reverse_compliments([]))
        design_rules.add_rule(build_strand_distance_library_rule(6,Library))

        #continue
        x = args.samples
        y = []

        primers = []
        for sims in args.samples:
            t,count = montecarlo(sims,primers,design_rules)
            if args.verbose:
                print (design_rules)
            
            if len(args.csv)>0:
                csv_file = open(args.csv,"a")
                hdr = get_experiment_header()
                vals = get_experiment_info(args)
                vals['codeword_size'] = cw
                vals['simulations'] = design_rules._total
                vals['estimated_primers'] = design_rules._passed
                y.append(design_rules._passed)
                entry = []
                for h in hdr:
                    entry.append("{}".format(vals[h]))

                csv_file.write(",".join(entry)+'\n')
                csv_file.close()

        def logfunc(x,c0,c1):
            return c0 + c1 * np.log(x+1)

        #def sigmoid(z, c0, c1,c2):
        #    return c0 + 1/(1 + c1*np.exp(-np.log(c2*z)))

        # def sqrtfunc(x,c0,c1):
        #     return c0 + c1 * np.sqrt(x)

        # def linear_func(x, a, b):
        #     return a*x + b

        funcs = [logfunc]

        for f in funcs:
            popt,pcov = curve_fit(f,np.array(x),np.array(y))
            size = f(4**20,*popt)
            print ("Using ({}) Extrapolated number of primers = {}".format(f.__name__,size))
            #print 'c0+c1*log(x+1) with c0=%5.3f, c1=%5.3f ' % tuple(popt)

import time
import random
import numpy as np

from doris import primer_util

import concurrent.futures
import threading
import multiprocessing
#import pd_accel

class RandomPrimerGenerator:
    def __init__(self, chars="AGCT", length=20):
        self.chars = chars
        self.len = length

    def get(self):
        return ''.join(random.choice(self.chars) for _ in range(self.len))

class UniquePrimerGenerator(RandomPrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[]):
        RandomPrimerGenerator.__init__(self,chars,length)
        self.avoid = library[:]

    def append(self, l):
        self.avoid.append(l)

    def _get_helper(self):
        return RandomPrimerGenerator.get(self)

    def get(self):
        while True:
            s = self._get_helper()
            if not (s in self.avoid):
                self.avoid.append(s)
                return s

class SpecificPrimerGenerator:
    def __init__(self,value, chars="AGCT", length=20):
        self.chars = chars
        self.len = length

    def get(self,v):
        l=''
        value=v
        for i in range(19,-1,-1):
            if(value>=3*4**i):
              l += "T"
              value-=3*4**i
            elif(value>=2*4**i):
              l += "C"
              value-=2*4**i
            elif(value>=1*4**i):
              l += "G"
              value-=1*4**i
            else:
              l += "A"
        #print value
        # print l
        return l

class LinearPrimerGenerator(SpecificPrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[]):
        SpecificPrimerGenerator.__init__(self,chars,length)
        self.avoid = library[:]

    def append(self, l):
        self.avoid.append(l)

    def _get_helper(self,i):
        return SpecificPrimerGenerator.get(self,i)

    def get(self,i):
        while True:
            s = self._get_helper(i)
            if not (s in self.avoid):
                self.avoid.append(s)
                return s

class LikelyPrimerGenerator(UniquePrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[], last='G', repl_factor=5):
        UniquePrimerGenerator.__init__(self,chars,length,library)
        self.last = last
        self.repl_factor = repl_factor

    def _get_helper(self):
        l = [random.choice(self.chars) for _ in range(self.len)]
        # force replication
        for j in range(self.len/self.repl_factor+1):
            i = random.randint(-2,3)
            if i+j*self.repl_factor+1 < self.len and i+j*self.repl_factor > 0:
                l[i+j*5+1] = l[i+j*5]

        l[-1] = self.last
        return "".join(l)

class MutatePrimerGenerator(UniquePrimerGenerator):
    def __init__(self, chars="AGCT", length=20, library=[], mutations=10):
        UniquePrimerGenerator.__init__(self,chars,length,library)
        self.mutations = mutations

    def _get_helper(self):
        if len(self.avoid)==0:
            l = [random.choice(self.chars) for _ in range(self.len)]
        else:
            l = self.avoid[np.random.randint(len(self.avoid))][:]
        ll = [ x for x in l ]
        rs = [ np.random.randint(self.len) for x in range(self.mutations) ]
        for r in rs:
            ll[r] = random.choice(self.chars)
        return "".join(ll)


class Rule:
    def __init__(self, r, name="", csv_name="", needs_args=False):
        self.r = r
        self.name = name
        self.passed = 0
        self.total = 0
        self.total_time = 0
        self.needs_args=needs_args
        self._csv_name = csv_name
        
    def run_rule(self, strand, **kwargs):
        if self.needs_args:
            return self.r(strand, **kwargs)
        else:
            return self.r(strand)

    def run(self, strand, **kwargs):
        self.total += 1 
        
        t1 = time.time()
        res = self.run_rule(strand, **kwargs)
        t2 = time.time()
        self.total_time += (t2-t1)

        if res == True:
            self.passed += 1
            return True
        else:
            return False

    @property
    def csv_value(self):
        return self.passed
        
    @property
    def csv_name(self):
        return self._csv_name
        
    def __str__(self):
        return "{:50}{:>5} / {:<7} \t Rate = {:<.2e} s/strand Time = {:<.2e} s".format(self.name,self.passed,self.total, self.total_time/(max(self.total,1)),self.total_time) 


class PrimerRule(Rule):
    def __init__(self, r, name="", csv_name=""):
        Rule.__init__(self,r,name,csv_name,True)
        #self._lock = threading.Lock()
            
    def run_rule(self, strand, **kwargs):
        assert 'primers' in kwargs
        for p in kwargs['primers']:
            if self.r(strand, p) == False:
                return False
        return True

class LibraryRule(Rule):
    def __init__(self, r, name="", csv_name="", Library=[]):        
        Rule.__init__(self,r,name,csv_name,False)
        self.Library = Library[:] # library doesn't change during analysis
        self.no_threads = False
        
    def run_rule(self, strand, **kwargs):
        for l in self.Library:
            if self.r(strand, l) == False:
                return False
        return True
    
def build_last_must_be_g_rule():
    return Rule(lambda s: s[-1] == 'G', "Last base must be G")

def build_lastfew_must_be_g_rule():
    return Rule(lambda s: s[-1] == 'G' or s[-2]=='G' or s[-3]=='G', "G must appear at end of strand")

def build_repetition_rule(fraction):
    return Rule(lambda s: not (primer_util.repetitionScore(s) < fraction),"Limit repetition to {}".format(fraction))

def build_numrepeats_rule(num):
    return Rule(lambda s: primer_util.countRepeats(s) >= num,"Require repetition of >= {}".format(num))


def build_self_complementarity_rule(x):
    return Rule(lambda s: not primer_util.self_complementarity(s,x), "Has no self complementing sequence")

def build_singlerun_rule():
    return Rule(lambda s: not primer_util.hasSingleRun(s), "Has no run")

def build_longrun_rule():
    return Rule(lambda s: not primer_util.hasLongRun(s), "Has no run of 4 or longer")

def build_dimerrun_rule():
    return Rule(lambda s: not primer_util.hasDimerRun(s), "Has no dimer run")

def build_GC_rule(f_lo,f_hi):
    return Rule(lambda s: primer_util.checkGC(s,(f_lo,f_hi)), "GC content between {} and {}".format(f_lo,f_hi))

def build_GC_at_end_rule():
    return Rule(lambda s: primer_util.checkGC(s[-5:],(0,60)), "GC at end does not exceed 60%")

def build_Tm_rule(T_lo, T_hi):
    return Rule(lambda s: primer_util.checkTm(s,(T_lo,T_hi)), "Tm between {} and {}".format(T_lo,T_hi))

def build_check_old_strands_rule():
    return Rule(primer_util.check_old_strands, "Check compatibility with s1/s2/s3")


def build_nextera_comparison_rule():
    return Rule(lambda s: primer_util.nextera_strand_comparison(s,3), "Avoid similarity withe Nextera primers")

def build_correlation_distance_library_rule(L):
    return PrimerRule(lambda s,l: not (primer_util.correlation_distance(s,l) > 4), "Correlation distance <= 4")

def build_hamming_distance_library_rule(distance,L):
    return PrimerRule(lambda s,l: not (primer_util.hamming_distance(s,l) < distance), "Hamming distance >= {}".format(distance))

def build_strand_distance_library_rule(distance,L):
    return LibraryRule(lambda s,l: not (primer_util.strand_primer_min_distance(l,s) <= distance),\
                       "Min distance > {} from library of {} strands".format(distance,len(L)),\
                       "",\
                       L)

def build_inter_complementarity_rule(distance,L):
    return PrimerRule(lambda s,l: not (primer_util.inter_complementarity(l,s,distance)),\
                       "Inter-complentarity less than {}".format(distance))

def build_primer3_heterodimer_bindings_library_rule(L):
    return PrimerRule(lambda s,l: primer_util.primer3_check_heterodimer(s,l),"Check heterodimers with primer3")

def build_primer3_homodimer_rule():
    return Rule(primer_util.primer3_check_homodimer, "Avoid homodimer with primer3")

def build_primer3_hairpin_rule():
    return Rule(primer_util.primer3_check_hairpin, "Avoid hairpin with primer3")

def build_primer3_check_reverse_compliments(L):
    return PrimerRule(lambda s,l: primer_util.primer3_check_reverse_compliments(s,l),"check reverse compliments")

class DesignRules:
    def __init__(self, name="", csv_name=""):
        self.rules = []
        self.name = name
        self.no_threads = True
        self._total = 0
        self._passed = 0
        if len(csv_name)==0:
            if len(name)>0:
                self.csv_name = name
            else:
                self.csv_name = "{}".format(random.randint(2**20))
        else:
            self.csv_name = csv_name
            
    def get_rules(self):
        return self.rules[:]

    def add_rule(self, r):
        self.rules.append(r)
        r.no_threads = self.no_threads

    def check(self, strand, **kwargs):
        self._total += 1
        for r in self.rules:
            if r.run(strand,**kwargs) == False:
                return False
        self._passed += 1
        return True

    def getTime(self):
        t = 0
        for r in self.rules:
            t += r.total_time
        return t

    def csv_header_string(self):
        hdr = "trial,primers_tested,"+",".join(self.rules[i].csv_name for i in range(0,len(self.rules)))
        return "{}\n".format(hdr)
        
    def csv_value_string(self):
        val = "{}".format(self.csv_name)
        val += ",{},".format(self._total)
        val += ",".join(str(self.rules[i].csv_value) for i in range(0,len(self.rules)))
        return "{}\n".format(val)
    
    def __str__(self):
        return "{} Results \n\t{}".format(self.name,"\n\t".join(str(self.rules[i]) for i in range(0,len(self.rules))))

def build_standard_design_rules(Library, with_nupack=True, with_primer3=False):
    dr = DesignRules("Standard Design Rules")
    #dr.add_rule(build_last_must_be_g_rule())
    #dr.add_rule(build_singlerun_rule())
    #dr.add_rule(build_dimerrun_rule())
    dr.add_rule(build_GC_rule(40,60))
    #dr.add_rule(build_repetition_rule(0.99))
    dr.add_rule(build_GC_at_end_rule())
    dr.add_rule(build_Tm_rule(45,60))
    dr.add_rule(build_hamming_distance_library_rule(10,Library))
    if with_nupack:
        dr.add_rule(build_nupack_homodimer_rule())
        dr.add_rule(build_nupack_nonspecific_bindings_library_rule(Library))
        dr.add_rule(build_nupack_heterodimer_bindings_library_rule(Library))
    elif with_primer3:
        dr.add_rule(build_primer3_homodimer_rule())
        dr.add_rule(build_primer3_heterodimer_bindings_library_rule(Library))
        dr.add_rule(build_primer3_check_reverse_compliments(Library))

    return dr


if __name__ == "__main__":
    import string
    import random
    import nupack

    def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    s = id_generator(size=20, chars="ACGT")
    print (s)
    print (primer_util.repetitionScore(s))

    r = build_repetition_rule(.99)
    r.run(s)
    print (str(r))

    dr = build_standard_design_rules([],True)
    dr.add_rule(r)
    dr.check(s)
    print (dr)

#!/usr/bin/python

import random
import csv

import editdistance as ed
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC

import primer3

def getTm(seq):
    return mt.Tm_NN(seq)

def getTmWallace(seq):
    return mt.Tm_Wallace(seq)

def checkTm(seq,trange):
    t = mt.Tm_NN(seq)
    return t >= trange[0] and t <= trange[1]

def checkGC(seq,gcrange):
    r = GC(seq)
    return r >= gcrange[0] and r <= gcrange[1] 

def read_primers(filename):
    f = open(filename,"r")
    tmp = f.readlines()
    f.close()
    L = []
    for l in tmp:
        l = l.strip()
        if l.startswith('%'):
            continue
        if len(l) == 0:
            continue
        L.append(l)
    return L

def read_resumed_runs(filename):
    f = open(filename,"r")
    trash = f.readline(); 
    tmp = f.readlines()
    f.close()
    L = []
    for l in tmp:
        l = l.strip()
        if l.startswith('%'):
            continue
        if len(l) == 0:
            continue
        L.append(l)
    return L

def mutate_sequence(seq, d):
    m = {'A': 'C', 'C' : 'T', 'G':'A', 'T':'G'}
    used = {}
    copy = [s for s in seq]
    assert d < len(seq)
    for i in range(0,d):
        while True:
            r = random.randint(0,len(seq)-2)
            #assert r != 19
            if r not in used:
                break
        used[r] = 1
        c = random.choice('ACGT')
        if copy[r] != c:
            copy[r] = c
        else:
            copy[r] = m[c]
    copy = "".join(copy)
    assert hamming_distance(seq,copy)==d
    return copy

def mutate_sequence_ed(seq, d):
    m = {'A': 'C', 'C' : 'T', 'G':'A', 'T':'G'}
    copy = [s for s in seq]
    assert d < len(seq)
    while True: 
        assert len(copy)==len(seq)
        delta = d - ed.eval(seq,"".join(copy))
        if delta < 0:
            copy = [s for s in seq]
        if delta == 0:
            break
        if delta >= 2:
            r = random.randint(0,3)
            if r==0:
                #print "did this!"
                # remove one
                r = random.randint(0,len(seq)-1)
                copy = copy[0:r] + copy[r+1:]
                r = random.randint(0,len(seq)-1)
                copy.insert(r,random.choice('ACGT'))
                assert len(copy)==len(seq)
            elif r==1:
                #print "did this!"
                # insert one
                r = random.randint(0,len(seq)-1)
                copy = copy[0:r] + [random.choice('ACGT')] + copy[r:-1]
                assert len(copy)==len(seq)
            else:
                # substitute one
                c = random.choice('ACGT')
                r = random.randint(0,len(seq)-1)
                if copy[r] != c:
                    copy[r] = c
                else:
                    copy[r] = m[c]
        else:
            c = random.choice('ACGT')
            r = random.randint(0,len(copy)-1)
            if copy[r] != c:
                copy[r] = c
            else:
                copy[r] = m[c]

    copy = "".join(copy)
    assert ed.eval(seq,copy)==d
    return copy

def create_anti_sequence(seq):
    m = {'A': 'C', 'C' : 'T', 'G':'A', 'T':'G'}
    copy = [m[s] for s in seq]
    copy = "".join(copy)
    return copy


def reverse_complement(seq):
    complement = {'T':'A', 'G':'C', 'C':'G', 'A':'T'}
    r = [complement[x] for x in seq]
    r.reverse()
    return "".join(r)

def reverse(seq):
    r = [x for x in seq]
    r.reverse()
    return "".join(r)


def self_complementarity(seq,l):
    for i in range(len(seq)-l):
        s = seq[i:i+l]
        r = reverse_complement(s)
        if r in seq[i+l:]:
            return True
    return False

def inter_complementarity(s1,s2,d):
    for i in range(len(s1)-d):
        s = s1[i:i+d]
        r = reverse_complement(s)
        if r in s2:
            return True
    return False

def repetitionScore(seq):
    score = 0
    for i in range(0,len(seq)-4):
        ss = seq[i:i+5]
        prev = ss[0]
        count = 0
        for s in ss[1:]:
            if prev==s:
                count = count + 1
            prev = s
        if count > 0:
            score = score + 1
    return float(score)/(len(seq)-4)

def countRepeats(seq):
    prev = seq[0]
    count = 0
    i = 0
    while i < len(seq):
        if seq[i]==seq[i-1]:
            count = count + 1
            if i+1 < len(seq) and seq[i+1] == seq[i]:
                # skip single repeating sequence, count as 1
                while i+1 < len(seq) and seq[i+1]==seq[i]:
                    i += 1
            else:
                i += 1
        else:
            i += 1
    return count

def hasSelfDimer(seq, n):
    r = reverse_complement(seq)
    for i in range(len(r)):
        s = r[i:i+n]
        if len(s)<n:
            break
        for j in range(i+1,len(seq)):
            if s == seq[j:j+n]:
                #print c,seq[j:j+n]
                return True
    return False

def hasSingleRun(seq):
    if seq.find("AAA") >= 0:
        return True
    if seq.find("TTT") >= 0:
        return True
    if seq.find("GGG") >= 0:
        return True
    if seq.find("CCC") >= 0:
        return True
    return False

def hasLongRun(seq):
    if seq.find("AAAA") >= 0:
        return True
    if seq.find("TTTT") >= 0:
        return True
    if seq.find("GGGG") >= 0:
        return True
    if seq.find("CCCC") >= 0:
        return True
    return False

def hasRepeat(seq):
    if seq.find("AA") >= 0:
        return True
    if seq.find("TT") >= 0:
        return True
    if seq.find("GG") >= 0:
        return True
    if seq.find("CC") >= 0:
        return True
    return False

def hasDimerRun(seq):
    if seq.find("ATATATAT") >= 0:
        return True
    if seq.find("TATATATA") >= 0:
        return True
    if seq.find("GCGCGCGC") >= 0:
        return True
    if seq.find("CGCGCGCG") >= 0:
        return True
    if seq.find("AGAGAGAG") >= 0:
        return True
    if seq.find("GAGAGAGA") >= 0:
        return True
    if seq.find("CACACACA") >= 0:
        return True
    if seq.find("ACACACAC") >= 0:
        return True
    if seq.find("TGTGTGTG") >= 0:
        return True
    if seq.find("GTGTGTGT") >= 0:
        return True
    if seq.find("CTCTCTCT") >= 0:
        return True
    if seq.find("TCTCTCTC") >= 0:
        return True
    return False

def hasShortDimerRun(seq):
    if seq.find("ATATAT") >= 0:
        return True
    if seq.find("TATATA") >= 0:
        return True
    if seq.find("GCGCGC") >= 0:
        return True
    if seq.find("CGCGCG") >= 0:
        return True
    if seq.find("AGAGAG") >= 0:
        return True
    if seq.find("GAGAGA") >= 0:
        return True
    if seq.find("CACACA") >= 0:
        return True
    if seq.find("ACACAC") >= 0:
        return True
    if seq.find("TGTGTG") >= 0:
        return True
    if seq.find("GTGTGT") >= 0:
        return True
    if seq.find("CTCTCT") >= 0:
        return True
    if seq.find("TCTCTC") >= 0:
        return True
    return False


def hamming_difference(s1, s2):
    value = {'A':0, 'C':1, 'G':2, 'T':3}
    diff = []
    for ch1, ch2 in zip(s1,s2):
        if ch1 == ch2:
            continue
        diff.append( str(ch1 + "->" + ch2) )
    return ",".join(diff)

def hamming_difference_indexes(s1, s2):
    value = {'A':0, 'C':1, 'G':2, 'T':3}
    diff = []
    i = 0
    for ch1, ch2 in zip(s1,s2):
        if ch1 == ch2:
            i += 1
            continue
        else:
            diff.append( i )
        i += 1
    return diff

def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def similarity(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 == ch2 for ch1, ch2 in zip(s1, s2))


def strand_primer_min_distance(strand, primer):
    assert len(strand) > len(primer)
    l1 = len(strand)
    l2 = len(primer)
    plen = len(primer)
    mn = max(l1,l2)
    if l1 - 2*plen > plen:
        for i in range(plen,l1-2*plen):
            d = hamming_distance(strand[i:i+plen],primer)
            #if d <= 6:
            #    print "{}\n{}".format(strand[i:i+plen],primer)            
            mn = min(d,mn)
            #if mn <= 6:
            #    break
            #print "{}: {} vs {}".format(mn,strand[i:i+l],primer[0:l])
        #print "{}: {} vs {}".format(mn,strand,primer)
        #import sys
        #sys.exit(0) 
        #print mn
        return mn
    else:
        return len(primer)
#
#
# The following functions help determine how similar two primers are
#
#
def show_correlation(s1, s2):
    #print "{} Correlations {}".format('-'*10,'-'*10)
    correlated_helper(s1,s2,True)
    correlated_helper(s2,s1,True)

def correlated_helper(s1, s2, show=False):
    l1 = len(s1)
    l2 = len(s2)
    for i in range(0,l1):
        l = min(l2,l1-i)
        s = sum(ch1 == ch2 for ch1, ch2 in zip(s1[i:i+l], s2[0:l]))
        if s1[i:i+l] == s2[0:l]:
        #if similarity(s1[i:i+l],s2[0:l])>7:
            if show:
                print ("")
                print ("   ",s1)
                t = ""
                for ch1,ch2 in zip(s1[i:i+l],s2[0:l]):
                    if ch1 == ch2:
                        t = t+"|"
                    else:
                        t = t+" "                    
                print ("   ",t.rjust(l+i))
                print ("   ",s2.rjust(len(s2)+i))
                print ("")
            return l
        #else:
            #print i, l1, l1-i, l2, l, s, s1[i:i+l], " ", s2[0:l]
    return 0

def uncorrelated(s1, s2):  
    return correlated_helper(s1,s2)==0 and correlated_helper(s2,s1)==0

def correlated(s1,s2):
    return not uncorrelated(s1,s2)

def correlation_distance(s1,s2):
    d1 = correlated_helper(s1,s2)
    d2 = correlated_helper(s2,s1)
    return max(d1,d2)

def check_old_strands(s):
    old_primers = ['CAGGTACGCAGTTAGCACTC','CGTGGCAATATGACTACGGA']
    for o in old_primers:
        if correlation_distance(s,o) > 3:
            return False
    return True

# should go somewhere else!

def checkComplexes(seqs,Tm=50):
    prefix = create_mfe_input(seqs,2)
    args = complexes_args(prefix)
    args.append("-T")
    args.append(str(Tm))
    complexes(args)
    c = read_mfe_output(prefix,".ocx-mfe")
    #for cc in c: 
        #print cc['deltaG']
    problems = [cc for cc in c if cc['deltaG'] < -10.0]
    return problems

def nupack_check_homodimer(s):
    c = checkComplexes([s,s])
    if len(c) > 0:
        #print "*****Homodimer: {} {}".format(s, c[0]['pattern'])
        return False
    else:
        return True

def nupack_check_heterodimer_complexes(s, l, Tm=50):
    # repeat s to search for homo-dimers
    c = checkComplexes([s,l],Tm)
    if len(c) > 0:
        #print "*****Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern'])
        return False        
    return True

def nupack_check_complexes(s, l, Tm=50):
    c = checkComplexes([l,reverse_complement(s)],Tm)
    if len(c) > 0:
        #print "*****(1) False binding ({}): {} vs {} -- {} {}".format(Tm,l,s, c[0]['pattern'], c[0]['deltaG'])
        return False        
    c = checkComplexes([s,reverse_complement(l)],Tm)
    if len(c) > 0:
        #print "*****(2) False binding ({}): {} vs {} -- {} {}".format(Tm,l,s, c[0]['pattern'],c[0]['deltaG'])
        return False        
    return True

def nupack_check_nextera_bindings(s):
    nextera_primers = get_nextera_primers()
    for l in nextera_primers:
        # repeat s to search for homo-dimers
        c = checkComplexes([s,l])
        if len(c) > 0:
            #print "*****Nextera Heterodimer: {} vs {} -- {}".format(l,s, c[0]['pattern'])
            return False

        c = checkComplexes([l,reverse_complement(s)])
        if len(c) > 0:
            #print "*****Nextera False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
            return False

        c = checkComplexes([s,reverse_complement(l)])
        if len(c) > 0:
            #print "*****Nextera Complement False binding: {} vs {} -- {}".format(l,s, c[0]['pattern'])
            return False
    return True

def nextera_strand_comparison(seq,distance):
    for i in illumina_primers:
        d = correlation_distance(seq,i)
        if d > distance:
            #show_correlation(seq,i)
            return False
    return True

def primer3_check_heterodimer(primer1,primer2):
    thermoResult=primer3.calcHeterodimer(primer1,primer2) 
    dg=thermoResult.dg/1000
    if(dg<-10): 
        return False
    else: 
        return True

def primer3_check_hairpin(primer): 
    thermoResult=primer3.calcHairpin(primer) 
    dg=thermoResult.dg/1000
    if(thermoResult.structure_found==True and dg<-10):
        return False 
    else: 
        return True
    
def primer3_check_homodimer(primer): 
    thermoResult=primer3.calcHomodimer(primer) 
    dg=thermoResult.dg/1000
    if(dg<-10): 
        return False 
    else: 
        return True

def primer3_check_reverse_compliments(primer1, primer2):
    rprimer1= reverse_complement(primer1)
    thermoResult=primer3.calcHeterodimer(primer2, rprimer1) 
    dg=thermoResult.dg/1000
    if(dg<-10): 
        return False 
    else: 
        return True



if __name__ == "__main__":
    from Bio.SeqUtils import GC
    if uncorrelated('ACG','TACG') == True:
        print ("Error")

    if correlated('GTCTCGTGGGCTCGG','TCGTCGGCAGCGTC') == True:
        print ("Error")

    print (reverse_complement("ACG"))
    
    s = repetitionScore('TCAAACATTTACGGGGCGAG')
    print ("Score = ", s)

    #5'-CAGGTACGCAGTTAGCACTC
    #5'-CGTGGCAATATGACTACGGA

    print ('-'*80)
    f = open("refined_ranked.txt","r")
    if f:
        s = f.readlines()
        for seq in s:
            seq.strip()        
            old_primers = ['CAGGTACGCAGTTAGCACTC','CGTGGCAATATGACTACGGA']
            for o in old_primers:
                if correlation_distance(seq,o) > 2:
                    show_correlation(seq,o)

            if GC(seq[-6:]) < .4:
                print ("Good GC ending:", seq)


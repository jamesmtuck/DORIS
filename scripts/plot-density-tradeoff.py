import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import csv
import sys
import math
import matplotlib.lines as mlines

def sqrtfunc(x,c0,c1,c2):
    return c2 + c1 * np.sqrt(c0*x)

def nrootfunc(x,c0,c1):
    return c1*np.power(x,c0)

def expfunc(x,c0,c1):
    return np.exp(-x/c0)+c1

def logfunc(x,c0,c1):
        return c0 + c1 * np.log(x+1)

def sigmoid(z, c0, c1, c2):
    return 1/(1 + np.exp(np.log(z*c0)*c1)) + c2

def linfunc(x,m,b):
    return x*m + b

def constfunc(x,m):
    return x*0+m

def collect_doris_data(pd):
    c = []
    for d in Data:
        if d['kind'] == 'doris' and d['pd']==pd:
            c.append( d['estimated_primers'] )
    return np.array(c)

def collect_primers_vs_samples_data(cwsize,libsize,bounds=None):
    y = []
    x = []
    for d in Data:
        if d['kind'] == 'primer_density_study' and int(d['lib_size'])==libsize and int(d['codeword_size'])==cwsize:
            sims = int(d['simulations'])
            if bounds is None:
                y.append( int(d['estimated_primers']) )
                x.append( int(d['simulations']) )
            elif sims >= bounds[0] and sims <= bounds[1]: 
                y.append( int(d['estimated_primers']) )
                x.append( int(d['simulations']) )

    print ("-"*50)
    print ("cwsize={} libsize={}".format(cwsize,libsize))
    print (x)
    print (y)
    return np.array(x),np.array(y)


def extrapolate_number_primers(x,y,f):
    popt,pcov = curve_fit(f,np.array(x),np.array(y))
    size = f(4**20,*popt)
    print ("Using ({}) Extrapolated number of primers = {}".format(f.__name__,size))
    #print 'c0+c1*log(x+1) with c0=%5.3f, c1=%5.3f ' % tuple(popt)
    print (popt)
    print (pcov)
    #return size
    return np.max(y)

def fit_curve(x,y,f):
    popt,pcov = curve_fit(f,np.array(x),np.array(y))
    print (pcov)
    #print "Using ({}) Extrapolated number of primers = {}".format(f.__name__,size)
    #print 'c0+c1*log(x+1) with c0=%5.3f, c1=%5.3f ' % tuple(popt)
    return popt


def calc_density(num_primers, cw, repl=10, filesize=10**9, payload=160):
    adj_file_size = filesize / repl
    total_strands = adj_file_size * num_primers
    size_of_index = math.ceil(math.log(adj_file_size,256))*cw
    bytes_per_strand = ( payload -  size_of_index ) / float(cw)
    return  total_strands * bytes_per_strand

def get_number_primers(x,y,x_val):
    y_list = []
    for xi,yi in zip(x,y):
        if xi==x_val:
            y_list.append(yi)

    assert len(y_list) > 0
    y_np = np.array(y_list)
    return np.mean(y_np),np.std(y_np)

    
def collect_primers_density_x_y(cwlist, libsize, filesize=10**9, payload=160):
    y_density = []
    x_cw = []
    y_primers = []
    
    for cw in cwlist:
        x,y = collect_primers_vs_samples_data(cw,libsize,(0,1405000))
        #print x,y
        if x.size == 0:
            continue

        p,sigma = get_number_primers(x,y,x_val=403500)

        print ("cw={} p={} sigma={}".format(cw,p,sigma))

        #print "-"*20
        density = calc_density(p,cw,filesize=filesize,payload=payload) / 10**12
        #density = file_size * p * 0.9 * ( 160 - math.ceil(math.log(file_size,256))*cw ) / float(cw) / 10**12 #TB
        #print "Density({}) = {:3.2f} TB/mL".format(cw,density)
        #print "Num. Primers ({}) = {}".format(cw,p)

        y_primers.append(p)
        y_density.append(density)
        x_cw.append(cw)

    return x_cw, y_primers, y_density



def plot_density(data,figname='out.pdf',sims=10**5):
    fig, (p1,p2) = plt.subplots(1, 2, figsize=(9, 4), sharex=False)        

    cw_list = [4,5,6,7,8,9,10,11,12]

    x_cw, y_primers, y_density = collect_primers_density_x_y(cw_list,sims,filesize=10**9)
    x_cw_density = np.add(160, -np.ceil(np.log(x_cw))*x_cw) / x_cw 
    
    p1.plot(np.array(x_cw_density),np.array(y_density),label='PCR')

    best_pcr_index = np.argmax(np.array(y_density))
    best_pcr_y = y_density[best_pcr_index]
    best_pcr_x = x_cw_density[best_pcr_index]
    
    d_x_cw, d_y_primers, d_y_density = collect_primers_density_x_y(cw_list,0,filesize=10**9,payload=157)

    best_overhang_index = np.argmax(d_y_density)
    best_overhang_y = d_y_density[best_overhang_index]
    
    #if len(d_y_density) > 0:
        #p1.plot(np.array(d_x_cw),np.array(d_y_density),label='DORIS (*)')
        #print d_x_cw
        #print d_y_primers
        
    #print "------Overhang----"

    d_x_density = np.add(160, -np.ceil(np.log(d_x_cw))*d_x_cw) / d_x_cw 
    best_overhang_x = d_x_density[best_overhang_index]
    
    p1.plot(np.array(d_x_density),np.array(d_y_density),\
             label='DORIS')
    
    p1.set_ylabel('Database Capacity (TB)')
    p1.set_xlabel('Density (B/strand)')

    speedup = "{:.2f}x".format( best_overhang_y/float(best_pcr_y))
    #p1.annotate(xytext=(best_x,(best_overhang-best_pcr)/2),text=speedup,)

    p1.text(best_overhang_x+0.1,(best_overhang_y+best_pcr_y)/2,speedup)
    p1.annotate("",xy=(best_overhang_x,best_overhang_y),xytext=(best_overhang_x,best_pcr_y),arrowprops=dict(arrowstyle='<->'))

    #p1.annotate("",xy=(best_x,best_pcr_y),xytext=(best_pcr_x,best_pcr_y),arrowprops=dict(arrowstyle='-'))

    #l = mlines.Line2D([best_overhang_x,best_pcr_x], [best_pcr_y,best_pcr_y], linestyle="--", color='black',linewidth=0.7)
    #p1.add_line(l)

    p1.annotate("",xy=(best_overhang_x,best_pcr_y),xytext=(best_pcr_x,best_pcr_y),arrowprops=dict(arrowstyle='<->'))
    speedup = "{:.2f}x".format( best_overhang_x/float(best_pcr_x))
    p1.text((best_overhang_x+best_pcr_x)/2,(best_pcr_y+.1),speedup)

    p2.plot(np.array(x_cw_density),np.array(y_primers))
    p2.plot(np.array(d_x_density),np.array(d_y_primers))

    #p2.plot(np.array(x_cw_density),np.array(y_primers),label="PCR")
    #p2.plot(np.array(d_x_density),np.array(d_y_primers),label="DORIS")

    p2.set_ylabel('Number of primers')
    p2.yaxis.set_label_position("right")
    p2.yaxis.tick_right()
    p2.set_xlabel('Density (B/strand)')

    print (d_y_density)
    print (y_density)
    print (y_primers)
    print (d_y_primers)
    
    fig.legend()
    fig.savefig(figname)
    #plt.show()


def plot_primer_prediction(cw,libsize,bounds,f=logfunc,ax=None):
    x,y = collect_primers_vs_samples_data(cw,libsize,bounds)

    miny = []
    maxy = []
    avgy = []
    newx = []

    x = list(x)
    y = list(y)
    
    for xi in x:
        print (xi, newx)
        if xi in newx:
            continue
        vals = []
        for xii,yi in zip(x,y):
            if xii == xi:
                vals.append(yi)
        miny.append(min(vals))
        maxy.append(max(vals))
        avgy.append(np.mean(np.array(vals)))
        newx.append(xi)

    if ax is None:
        fig, (ax) = plt.subplots(1, 1, figsize=(9, 4), sharex=True)        
        

    def newf(x,c0,c1,c2):
        return c2 + c1 * np.sqrt(c0*x)

    ff = f
    
    v = extrapolate_number_primers(newx,miny,ff)
    #miny.append(v)

    v = extrapolate_number_primers(newx,maxy,ff)
    
    #maxy.append(v)

    v = extrapolate_number_primers(newx,avgy,ff)

    p = fit_curve(newx,avgy, ff)

    storage_size = np.logspace(start=1,stop=6,num=1000,base=10)
    primers = ff(storage_size,*p)
    
    #avgy.append(v)

    #newx.append(4**20)

    #p1.set_xscale('log')
    ax.plot(newx,miny,'x',label="min {}".format(cw))
    ax.plot(newx,maxy,'x',label="max {}".format(cw))
    ax.plot(newx,avgy,'x',label="avg {}".format(cw))
    ax.plot(storage_size,primers,label="{}-{}".format(cw,libsize))

    ax.legend()
    return ax
    
    
if __name__ == "__main__":

    Data = []
    for fname in sys.argv[1:]:
        with open(fname, 'r') as csvfile:
            data = csv.DictReader(csvfile)
            for row in data:
                Data.append(row)

    #plot_density(Data,sims=10**4,figname='CodewordVsDensity4.png')
    #plot_density(Data,sims=10**5,figname='CodewordVsDensity5.png')
    plot_density(Data,sims=10**6,figname='CodewordVsDensity.png')

    


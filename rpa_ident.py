# -*- coding: utf-8 -*-
"""
Created on Thu May 26 18:17:48 2022

@author: JensJ
"""
import re
#import io
#import zipfile
#import pandas as pd
import numpy as np
from  scipy.optimize import leastsq
import pandas as pd
import gzip


def find_tables(txt):#### to do?
    txt=txt.replace('\r','')
    txt=txt.replace('\n','')
    rexh=r'<h2>(.*?)</h2>'
    rext=r'<table.*?/table>'
    headers = re.findall(rexh,txt)#,re.DOTALL)#[:1]
    tables  = re.findall(rext,txt)#[:1]
    rexc=r'Compound: (.*?) Order No: (\d{8}-\d{4}?)(.*?)</td>'
    comps=re.findall(rexc,txt)
    headerss=[sub.replace("'",'dash').replace('*','star') for sub in headers]
    
    return headerss,tables,comps

def get_tdata(t):
    t=t.replace(',','')
    rexd=r'<td>(.*?)</td>'
    tdata = re.findall(rexd,t)
    return tdata

def create_sub(h,t,c):
    filename=h+c[-1].replace('(','').replace(')','').replace(':','_') ##header + (datetime)
    tdata = get_tdata(t)
    tdata[::2] = [td+',' for td in tdata[::2]]
    tdata[1::2] = [td+'\n' for td in tdata[1::2]]
    #tdata[1::2] = tdata[::2]+'\n'
    data = [tdata[0],h+'\n']
    data.extend(tdata[2:])
    return filename,data

def save_subs(headers,tables,comps):
    files=[]
    dd={}
    for h,t,c in zip(headers,tables,comps):
        filename,data = create_sub(h,t,c)
        
        dd[filename]=get_tdata(t) # no comma, no \n
        
        txt = f'# {c[0]} order {c[1]}, date {c[2]}\n'
        txt += ''.join(data)
        if mysystem=='web':
            text_file = anvil.BlobMedia('text/plain', txt.encode('utf-8'), name=filename+'.csv')
            anvil.media.download(text_file)
        else:
            with open(filename+'.csv', 'w+') as ff:
                ff.write('%s'%txt)
        
            
        print(h,c)
        #with open(filename,'w+') as f:
        #  f.write('#%s\n'%(h))
        #zf.write(filename)
        #
        #print ('written %s'%filename)
        #shutil.make_archive('mya','zip','/tmp')
        #zf.close()
            
        files.append(filename)
    return files,dd




        
def synchronize(sdash,nstar,temp,etarate=1):
    # create a data array with [time_s,S',n*,TKinv,etarate]
    # 
    # sdas, nstar , and temp are list of tuples (time, value)
    # etarate is a scalar
    #
    sda = np.array(sdash)
    nst = np.array(nstar)
    tem = np.array(temp)
    time = sda[:,0]
    sda = sda[:,1]
    heatrate = (tem[-1,1]-tem[0,1])/(tem[-1,0]-tem[0,0])

    ndata=len(sda)
    
    
    
    if len(nst)!= ndata:
        nda = np.interp(time, nst[:,0], nst[:,1])
    else:
        nda =nst[:,1]
    # inverse temperature K
    tempc = np.interp(time, tem[:,0],tem[:,1])
    tki = 1./(tempc +273.15)
    
    time*=60 #in sec
    trate = np.ones_like(time)*heatrate/60  # per sec
    eta = np.ones_like(time)*etarate*2*np.pi  # per sec
    dataset =np.vstack((time,sda,nda,tempc,tki,eta,trate))
    return dataset.reshape(7,-1).T

    
    


def read_html(htmlfile=None,btxt=None):
    '''
    read htmlfile data
     with measurements from rpa

    Parameters
    ----------
    htmlfile : str, optional
        filename. The default is None.
    btxt : bytes, optional
        compressed content of htmlfile. The default is None.

    Returns
    -------
    headers : list 
        List of headers per dataset (measure)        
    tables : list
        list of str with html table 
    comps : list
        tuples of metadata for measurements 

    '''
    # measurements are made
    #@etarates(1.25, 10, 5, 2.5) # 1/s  # not recorded, must be in this sequence
    #@heatrates(1/2,1/3,1/6,1/12) K/s # can be computed from temperature
    # gives 16 measurements of time vs. (Temp,nstar and Sdash)
    # makes 16 combinations for 3 curves, yields 48 curves
    if htmlfile!="":
        try:
            with open(htmlfile,'rb') as fd: #encoding='utf-8'
                btxt=fd.read()
            print(len(btxt))
            txt = btxt.decode()    
        except:
            print('seems to be zipped')
            try:
                with gzip.open(htmlfile,mode="rt") as fz:
                    txt=fz.read()
            except:
                print('also zip does not work for %s'%htmlfile)
    else:
        assert btxt is not None, 'No file and no content given!'
        try:
            txt = gzip.decompress(btxt).decode('utf-8')
        except:
            txt = btxt
    
    print(txt[:120])
    #txt = base64.b64encode(enc_str)
    headers,tables,comps=find_tables(txt)  
    return (headers,tables,comps)

def stack_rpa_data(headers,tables,comps,etarates_sequence=[1.25,10,2.5,5]):
    #fnames,dd = save_subs(headers,tables,comps)
    # sequence of rates 1/s
    etarates = np.array(etarates_sequence*12).reshape(12,4).T.ravel()
    ic = 0
    datafield = np.ndarray(shape=(0,7))
    for h,t,c,er in zip(headers,tables,comps,etarates):
        data = get_tdata(t)
        filename,datax = create_sub(h,t,c)
        #print(filename,len(datax))
        #txt = ','.join(data[2:])
        mya = np.array(data[2:],dtype=np.float32).reshape((-1,2))
        if ic==0:
            
            if 'Sdash' not in  filename:
                print (filename,c)
            sda = mya   
            
            #print(c[-1])
        if ic ==1:
            if 'Temp' not in  filename:
                print (filename,c)
            #if c[-1]!= ordertime:
            #    print (h,c)
            temp = mya
        if ic==2:
            #if c[-1]!= ordertime:
            #    print (h,c)
            if 'nstar' not in  filename:
                print (filename,c)

            nstar = mya
            
            newset = synchronize(sda, nstar, temp,er)
            datafield = np.vstack((datafield,newset))
            ic=-1
            
        ic+=1
    df = pd.DataFrame(datafield,columns=['time',
                                          'sdash',
                                          'nstar',
                                          'tempc',
                                          'tki',
                                          'eta',
                                          'hrate'],
                      )
                
    return df


def fit_visco(df,lowert=80,uppert=120):
    dfs=df.loc[(df.tempc<=uppert)&(df.tempc>=lowert)]

          
    p0=[5,1936,0.297]
    x=np.log(dfs.eta.values)
    y=dfs.tki.values
    z=np.log(dfs.nstar.values)
    
    result = leastsq(f_visco,p0,args=(x,y,z))    
    ((A,C,n),pp)=result
    #print(A,C,n,pp)
    return {'A':A,'C':C,'n':n,'pp':pp}

def viscosity(x,y,A,C,n):
    return A + C*y + (n-1)*x

def f_visco(p,x,y,z_data):
    A=p[0]
    C=p[1]
    n=p[2]
    z_model = viscosity(x,y,A,C,n)
    diffs = z_model-z_data
    return diffs.flatten()



if __name__=='__main__':
    import os
    htmlfile=r"Rheology M-870-6 Batch 31.html"
    #htmlfile=r"Rheology M-870-6 Batch 31.zip"
    fullpath = os.path.join('testdata',htmlfile)

    with open(fullpath) as fi:
        txt=fi.read()
    btxt = gzip.compress(bytes(txt, 'utf-8'))
    h,t,c =  read_html("",txt)
    df = stack_rpa_data(h,t,c)
    df.to_csv(htmlfile[:-5] + '.csv', index=False)
    
    print(fit_visco(df))
    
    
    
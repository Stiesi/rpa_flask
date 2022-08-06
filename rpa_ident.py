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
    # check 16 x time in html file
    nsets = len(re.findall('time',txt,re.IGNORECASE))
    assert nsets==48; "html file has not 4 x 4 x [S', S'', n*] time sets"
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
    meta = dict(name=c[0],order=c[1],date=c[2])
    tdata[::2] = [td+',' for td in tdata[::2]]
    tdata[1::2] = [td+'\n' for td in tdata[1::2]]
    #tdata[1::2] = tdata[::2]+'\n'
    data = [tdata[0],h+'\n']
    data.extend(tdata[2:])
    return filename,data,meta

def save_subs(headers,tables,comps):
    files=[]
    dd={}
    for h,t,c in zip(headers,tables,comps):
        filename,data,meta = create_sub(h,t,c)
        
        dd[filename]=get_tdata(t) # no comma, no \n
        
        txt = f'# {c[0]} order {c[1]}, date {c[2]}\n'
        txt += ''.join(data)
        #if mysystem=='web':
        #    text_file = anvil.BlobMedia('text/plain', txt.encode('utf-8'), name=filename+'.csv')
        #    anvil.media.download(text_file)
        #else:
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


        
def synchronize(sdash,nstar,temp,testid=1,trigger=0):
    # create a data array with [time_s,S',n*,TKinv,gammarate]
    # 
    # sdas, nstar , and temp are list of tuples (time, value)
    # gammarate is a scalar
    #
    sda = np.array(sdash)
    nst = np.array(nstar)
    tem = np.array(temp)
    
    if trigger==0:
        time = sda[:,0]
    else:        
        time = np.linspace(sda[0,0],sda[-1,0],num=trigger)
    
    #heatrate = (tem[-1,1]-tem[0,1])/(tem[-1,0]-tem[0,0])
    # or
    heatrate=(np.diff(tem[:,1])/np.diff(tem[:,0])).mean()

    #sda = sda[:,1]
    ndata=len(time)
    sda = np.interp(time, sda[:,0], sda[:,1])
    
    
    
    if len(nst)!= ndata:
        nda = np.interp(time, nst[:,0], nst[:,1])
    else:
        nda =nst[:,1]
    # inverse temperature K
    tempc = np.interp(time, tem[:,0],tem[:,1])
    #tki = 1./(tempc +273.15)
    
    time*=60 #in sec
    trate = np.ones_like(time)*heatrate/60  # per sec
    #trate = np.ones_like(time)*heatrate  # per min
    #gammap = np.ones_like(time)*gammarate*2*np.pi  # per sec
    testnr = np.ones_like(time)*testid
    dataset =np.vstack((time,sda,nda,tempc,trate,testnr))
    return dataset.reshape(6,-1).T

# transform dataframes 
def dfto_json(df):
    dfjson = df.to_json(orient='split')
    return dfjson

def dfread_json(dfjson):
    df = pd.read_json(dfjson,orient='split')
    return df


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
    #@gammarates(1.25, 10, 5, 2.5) # 1/s  # not recorded, must be in this sequence
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
    
    #print(txt[:120])
    #txt = base64.b64encode(enc_str)
    headers,tables,comps=find_tables(txt)  
    return (headers,tables,comps)

def stack_rpa_data(headers,tables,comps,
                   gammarates=[1.25,2.5,5,10],
                   heatrates=[5,10,20,30],
                   trigger=1000):
    #fnames,dd = save_subs(headers,tables,comps)
    # sequence of rates 1/s
    # sorting of sequences is unimportant (see below), program will find
    # gammarates according to nstar at different rates @100°C
    # heatrates from temperature tempc (number values are arbitrary (but must be different))
    #gammarates = np.array(gammarates*12).reshape(12,4).T.ravel()
    ic = 0
    datafield = np.ndarray(shape=(0,6))
    testid=0
    for h,t,c in zip(headers,tables,comps):
        data = get_tdata(t)
        filename,datax,meta = create_sub(h,t,c)
        #print(filename,len(datax))
        #txt = ','.join(data[2:])
        mya = np.array(data[2:],dtype=np.float32).reshape((-1,2))
        if ic==0:
            testid+=1
            
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
            
            
            newset = synchronize(sda, nstar, temp,testid,trigger=trigger)
            datafield = np.vstack((datafield,newset))
            ic=-1
            
        ic+=1
    
    df = pd.DataFrame(datafield,columns=['time',
                                          'sdash',
                                          'nstar',
                                          'tempc',
                                          #'tki',
                                          #'gammap',
                                          'hrate',
                                          'test_no'],
                      )
    #df['gammap']=pd.cut(df.gammap,4,labels=gammarates_sequence)
    #for test_no in df.test_no.unique():
    #    dfs=df[df.test_no==test_no]
        #dfs[]
    # group test by nstar@100°C
    xx=df[(df.tempc>99)&(df.tempc<101)].groupby(df.test_no).nstar.mean()
    df=df.merge(xx,on='test_no')
    # create increasing rates
    gammarates.sort(reverse=True)
    gammaratesa=np.array(gammarates)*2*np.pi
    # associate gammarates with values of nstar at 100°C
    df['gp100']=pd.qcut(df.nstar_y,4,labels=gammaratesa)
    # habe nur 3 categorien !?
    #df['gammap']=pd.qcut(df.nstar_y,4,labels=gammarates)
    df.rename(columns={'nstar_x':'nstar'},inplace=True)
    # set heatrates to given categories (unique values)
    df['hrate']=pd.qcut(df.hrate,4,labels=heatrates)
    df['hrate']=df['hrate'].astype('float')
    df.nstar*=1e5
    df['gammap']=df['gp100'].values.astype('float')
    df['test_no']=df['test_no'].astype('int')
    df.drop(['gp100','nstar_y'],axis=1,inplace=True)
    df = df.round({'nstar':1,'gammap':2,'tempc':2,'sdash':3,'time':2})
    return df


def fit_visco(df,lowert=80,uppert=140):
    dfs=df.loc[(df.tempc<=uppert)&(df.tempc>=lowert)]

          
    p0=[300,1936,0.297]
    x=np.log(dfs.gammap.values)
    y=1./(dfs.tempc.values+273.15)
    z=np.log(dfs.nstar.values)
    #z=dfs.nstar.values
    
    full_output=0    
    result = leastsq(f_visco,p0,args=(x,y,z),full_output=full_output)
    if full_output!=0:
        print(result)
    A,C,nexp=result[0]
    pp=result[-1]
        
    #print(A,C,n,pp)
    return {'A':A,'C':C,'n':nexp,'pp':pp,'lower T[C]':lowert,'upper T[C]':uppert}

def viscosity(loggammap,tki,A,C,n):
    return A + C*tki + (n-1)*loggammap


def viscosity_log(loggammap,tki,A,C,n):
    return np.log(A) + C*tki + (n-1)*loggammap

def f_visco(p,loggammap,tki,log_nstar):
    A=p[0]
    C=p[1]
    n=p[2]
    lognstar_model = viscosity_log(loggammap,tki,A,C,n)
    diffs = lognstar_model-log_nstar
    return diffs.flatten()

def plot(df,para,title='RPA',filename='nstar'):
    A=para['A']
    C=para['C']
    n=para['n']
    lowert = para['lower T[C]']
    uppert = para['upper T[C]']
    df['x']=np.log(df.gammap)
    df['y']=1./(df.tempc+273.15)
    # measured values
    df['z']=np.log(df.nstar)
    #df['z']=df.nstar
    # model solution (analytical)
    df['za']=df.apply(lambda x: viscosity_log(x.x,x.y,A,C,n),axis=1)

    dfs = df[(df.tempc>lowert-10)&(df.tempc<uppert+20)]
    
    ax=dfs.plot(x='tempc',y='z',
                label='RPA',
                grid=True,
                #color='gammap',
                style='.')
    dfs.plot(x='tempc',y='za',
             label=f'model A={A:.2f}, C={C:.1f}, n={n:.4f}',
             style='.',
             grid=True,ax=ax,ylabel='log n*',title=title,xlabel='Temperature in °C',
             )
    fig=ax.get_figure()
    fig.savefig(filename+'.png')
    

if __name__=='__main__':
    import os
    htmlfile=r"Rheology M-870-6 Batch 31.html"
    #htmlfile=r'M-1470-224 Viscosity Data.html'
    #htmlfile=r"Rheology M-870-6 Batch 31.zip"
    fullpath = os.path.join('testdata',htmlfile)

    with open(fullpath) as fi:
        txt=fi.read()
    btxt = gzip.compress(bytes(txt, 'utf-8'))
    h,t,c =  read_html("",txt)
    df = stack_rpa_data(h,t,c,trigger=1000)# gammarates_sequence=[1.25,2.5,5.,10.])
    # to zip
    zipfile = htmlfile[:-5]+'.zip'
    compression_options = dict(method='zip', archive_name=htmlfile[:-5] + '.csv')
    
    df.to_csv(zipfile,compression=compression_options,index=False)
    
    #df,res = fitdata(savename)        

    para=fit_visco(df)
    print(para)
    A=para['A']
    C=para['C']
    n=para['n']
    df['x']=np.log(df.gammap)
    df['y']=1./(df.tempc+273.15)
    df['z']=np.log(df.nstar)
    #df['z']=df.nstar
    df['za']=df.apply(lambda x: viscosity_log(x.x,x.y,A,C,n),axis=1)
    df['nstara']=np.exp(df.za)
    df['za_matlab']=df.apply(lambda x: viscosity_log(x.x,x.y,336,1735,0.2395),axis=1)
    #df['za_matlab']=df.apply(lambda x: viscosity_log(x.x,x.y,0.0048,1882,0.275),axis=1)
    if 0:
        for testno in df.test_no.unique():
            dfs = df[df.test_no==testno]
            #print(len(dfs),dfs.gammap.iloc[0],dfs.hrate.iloc[0])
            gammap= dfs.gammap.iloc[0]
            hrate = dfs.hrate.iloc[0]
            
            ax=dfs.plot(x='tempc',y='z',label='%d: h=%.2f, gp=%.2f'%(testno,hrate,gammap),grid=True)
            dfs.plot(x='tempc',y='za',label='python h=%.2f, gp=%.2f'%(hrate,gammap),grid=True,ax=ax)
            #dfs.plot(x='tempc',y='za_matlab',label='mlab h=%.2f, gp=%.2f'%(hrate,gammap),grid=True,ax=ax)
    
    for gammap in df.gammap.unique():
        dfs = df[(df.gammap==gammap)&(df.tempc>80)&(df.tempc<140)]
        #print(len(dfs),dfs.gammap.iloc[0],dfs.hrate.iloc[0])
        gammastr= dfs.gammap.iloc[0]
        hrate = dfs.hrate.iloc[0]
        
        ax=dfs.plot(x='tempc',y='z',
                    label=' h=%.2f, gp=%.2f'%(hrate,gammastr),
                    grid=True,
                    style='.')
        dfs.plot(x='tempc',y='za',
                 label='python h=%.2f, gp=%.2f'%(hrate,gammastr),
                 style='.',
                 grid=True,ax=ax)
        
        

    
    
    plot(df,para,title=htmlfile)
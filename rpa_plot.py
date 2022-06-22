# -*- coding: utf-8 -*-
"""
Created on Thu May 26 18:17:48 2022

@author: JensJ
"""
import numpy as np
import pandas as pd

import plotly.graph_objects as go

import rpa_ident as ri

DEFAULT_PLOTLY_COLORS=['rgb(31, 119, 180)', 'rgb(255, 127, 14)',
                       'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
                       'rgb(148, 103, 189)', 'rgb(140, 86, 75)',
                       'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
                       'rgb(188, 189, 34)', 'rgb(23, 190, 207)']

def plot(df,sx,sy,sz,title='RPA',fitstring='fit',filename='RPA'):
    traces=[]
    # create colors and index from heatrate
    # equal number of values in categories
    categories = 4
    df['col']=pd.qcut(df.hrate,categories,labels=range(categories))
    for testno in df.test_no.unique():
        dfs = df[df.test_no==testno]
        z=np.log(dfs.nstar.values)        
        x=dfs.tempc.values
        y=dfs.gammap.values
        
    
        hr = dfs.hrate.mean()
        col = dfs.col.values[0]
        #print(col)
        gammap = dfs.gammap.values[0]
        traces.append(go.Scatter3d(z=z, x=x, y=y,
                                   legendgroup='%.1f'%col,
                                   name=f'{testno} h {hr:.2f} gp {gammap:.2f}',
                                   mode='markers',
                                   marker = dict(color=DEFAULT_PLOTLY_COLORS[col],
                                                 #colorscale='Viridis',   # choose a colorscale
                                                 #opacity=0.8
                                                 )
                                   ))
    
    traces.append(go.Surface(x=sx,y=sy,z=sz,name=fitstring,showscale=False))
        
    fig = go.Figure(data=traces)
    fig.update_layout(title=title,
                      scene=dict(
                          xaxis_title='temp in C',
                          yaxis_title='shear rate in 1/s',
                          zaxis_title='log (n*)'),
                      autosize=True,
                      width=1500, height=1500,
                      #margin=dict(l=65, r=50, b=65, t=90)
                      )
    #fig.show()
    fig.write_html(f'{filename}.html', auto_open=True)
    


if __name__=='__main__':
    import os
    #htmlfile=r"Rheology M-870-6 Batch 31.html"
    zipfile=r"Rheology M-870-6 Batch 31.zip"
    compression_options = dict(method='zip', archive_name=zipfile[:-4] + '.csv')
    df = pd.read_csv(zipfile,compression='zip')
    
    print (df.info())
    
    res = ri.fit_visco(df,lowert=80,uppert=120)
    A=res['A']
    C=res['C']
    n=res['n']
    parstr = f'log(n*) = log({A:.3f}) + {C:.1f}/(tempK) +  ({n:.3f}-1)*log(gamma)'
    # create plot
    gammap = np.linspace(df.gammap.min(),df.gammap.max(),num=11)
    tempc= np.linspace(80,120,num=11)
    tki = 1./(tempc+273.15)
    gpa, tkia = np.meshgrid(gammap,tki)
    logz = ri.viscosity_log(np.log(gpa),tkia,A,C,n)
    sr,tc = np.meshgrid(gammap,tempc)
    plot(df,tc,sr,logz,title=parstr,fitstring=parstr,filename=zipfile[:-4])

import os
import json
from io import BytesIO
import uuid
import zipfile

#from deta import Deta
from flask import Flask, render_template, request,send_file,send_from_directory
from werkzeug.utils import secure_filename

import rpa_ident as rid

app = Flask(__name__)
dir,file=os.path.split(__file__)
app.config['UPLOAD_FOLDER']=os.path.join(dir,'upload')
app.config['MAX_CONTENT_PATH']=16000000


#db = deta.Base('test')

@app.route('/')
def index():
   return render_template('upload.html')

@app.route('/upload')
def upload():
   return render_template('upload.html')

@app.route('/uploader', methods = ['GET', 'POST'])
def upload_file():
   if request.method == 'POST':

      unique_dirname = str(uuid.uuid4())
      os.chdir(os.path.join(app.config['UPLOAD_FOLDER']))
      os.mkdir(unique_dirname)
      relative_dir = os.path.join(app.config['UPLOAD_FOLDER'],unique_dirname)

      f = request.files['file']
      savename=os.path.join(relative_dir,secure_filename(f.filename))
      basename=os.path.basename(savename) # filename without path
      base,ext = os.path.splitext(basename) # filename without extension
      #f.save(f.filename)
      f.save(savename)
      print(f'file {f.filename} uploaded successfully as {savename}')

      df,res = fitdata(savename)   
      #figfile = os.path.join(relative_dir,base)
      #rid.plot(df,res['fit_parameter'],title=basename,filename=figfile)
      
      arc_file = get_zip_data(df,res,base,relative_dir)
      
      print(unique_dirname, arc_file)
      return send_from_directory(relative_dir,arc_file,as_attachment=True) 
      
def get_zip_data(df,meta,base,relative_dir):
    
    os.chdir(relative_dir)    
    #outname = base+'.zip'
    # write df to zip
    #fullname=os.path.join(unique_dirname,outname)
    arcfile = f'{base}.zip'
    metafile= base+'.txt'
    #figfile= base+'.png'
    compression_options = dict(method='zip', archive_name=f'{base}.csv')
    df.to_csv(arcfile, compression=compression_options, index=False)
    # write meta dict to file
    with open(metafile,'w') as fm:
        print(meta,file=fm)
    zf = zipfile.ZipFile(arcfile,mode='a')
    zf.write(metafile,arcname=base+'.txt')
    #zf.write(figfile,arcname=base+'.png')
    zf.close()
    return arcfile
    
    
def fitdata(filename):#,gammarates_sequence=[1.25,10,2.5,5]):
    with open(filename) as fh:
        txt=fh.read()
    
    # convert to pandas
    #h,t,c =  rid.read_html(None,htmlfile)
    headers,tables,comps=rid.find_tables(txt)
    df = rid.stack_rpa_data(headers,tables,comps)
    #df.to_csv(htmlfile[:-5] + '.csv', index=False)
    
    # do fitting of data
    try:
        parameters = rid.fit_visco(df)
    
    except:
        parameters={'Error in Fit'}
        print('Error in Fitting')
    #pj = json.dumps(parameters)
    #down(df,filename)
    # here we set our output parameters in the form of a json
    
    return df,{
            'fit_parameter': parameters ,
            'flow function': 'log(n_star) = log(A) + C * (1/temperature[K]) + (n-1)* log(gammad[rad/s])',
            #'lower T[C]' : '80',
            #'upper T[C]' : '120',
            }
            
@app.route('/downloader', methods = ['GET', 'POST'])    
def down(df,filename):
    basename=os.path.basename(filename)
    base,ext = os.path.splitext(basename)
    outname = base+'.csv'
    print(outname)
    #pdstr = df.to_string()
    #buffer = BytesIO()
    #return  send_file(BytesIO(df.to_csv(index=False).encode()),mimetype='text/csv',filename=outname) 
    return  send_file(BytesIO(df.to_csv(index=False).encode()),mimetype='text/csv',attachment_filename=outname) 
      
      
if __name__ == '__main__':
   app.run(debug = True)



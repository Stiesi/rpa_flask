
import os
import json

from flask import Flask, render_template, request
from werkzeug.utils import secure_filename

import rpa_ident as rid

app = Flask(__name__)
dir,file=os.path.split(__file__)
app.config['UPLOAD_FOLDER']=os.path.join(dir,'upload')
app.config['MAX_CONTENT_PATH']=16000000


@app.route('/')
def index():
   return render_template('upload.html')

@app.route('/upload')
def upload():
   return render_template('upload.html')

@app.route('/uploader', methods = ['GET', 'POST'])
def upload_file():
   if request.method == 'POST':
      f = request.files['file']
      savename=os.path.join(app.config['UPLOAD_FOLDER'],secure_filename(f.filename))
      f.save(savename)
      #f.save(f.filename)
      print('file uploaded successfully')
      res = fitdata(savename)
      return res
      
def fitdata(filename):
    
    with open(filename) as fh:
        txt=fh.read()
    eta_sequence=[1.25,10,2.5,5]
    # convert to pandas
    #h,t,c =  rid.read_html(None,htmlfile)
    headers,tables,comps=rid.find_tables(txt)
    df = rid.stack_rpa_data(headers,tables,comps,etarates_sequence=eta_sequence)
    #df.to_csv(htmlfile[:-5] + '.csv', index=False)
    
    # do fitting of data
    try:
        parameters = rid.fit_visco(df)
    
    except:
        parameters={'Error in Fit'}
        print('Error in Fitting')
    pj = json.dumps(parameters)

    # here we set our output parameters in the form of a json
    return {#'rpa_data': df.to_json(), 
            'fit_parameter':pj
            }        
        
      
      
if __name__ == '__main__':
   app.run(debug = True)



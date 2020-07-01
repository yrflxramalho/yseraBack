from flask import Flask, url_for, send_from_directory, request, render_template
import os
import json
from werkzeug.utils import secure_filename
from ysera import myfunction, ysera
from flask_cors import CORS
import time


app: Flask = Flask(__name__)
CORS(app)

PROJECT_HOME = os.path.dirname(os.path.realpath(__file__))


@app.route('/ysera', methods = ['POST'])
def getFile():
    f = request.files['file']
    
    params = {
        "hb" : float(request.form.get("hb")),
        "sb" : float(request.form.get("sb")),
        "db" : float(request.form.get("db")),
        "vdw": float(request.form.get("vdw")),
        "aaanBeg" : float(request.form.get("aaanBeg")),
        "aaanEnd" : float(request.form.get("aaanEnd")),
        "aaspi" : float(request.form.get("aaspi")),
        "aactnBeg" : float(request.form.get("aactnBeg")),
        "aactnEnd" : float(request.form.get("aactnEnd")),

    }
    
    millis = int(round(time.time() * 1000))
    name = 'file_'+str(millis)+'.pdb'
    f.save(PROJECT_HOME+"/temp/"+name)
    
    x = ysera(name, params)
    
    # x = "file_1592291221139.pdb;hb=0;sb=22;db=10;lpi=121,tshaped=0,inter=0,paralel=121;vdw=999;ctn=4;an=0;spi=11"
    
    return json.dumps(x)

@app.route('/record', methods=['POST'])
def download():
    data = request.get_json()
    print(data) 
    
    return send_from_directory(directory='output', filename=data['filename'])
    

app.run(host='0.0.0.0', debug=True)
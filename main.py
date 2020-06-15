from flask import Flask, url_for, send_from_directory, request, render_template
import os
import json
from flask_cors import CORS

app: Flask = Flask(__name__)
CORS(app)

PROJECT_HOME = os.path.dirname(os.path.realpath(__file__))

@app.route('/', methods=['POST'])
def api_root():
   
    return json.dumps('response')

  
@app.route('/result/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    return send_from_directory(directory='RESULT_FOLDER', filename=filename)



app.run(host='0.0.0.0')
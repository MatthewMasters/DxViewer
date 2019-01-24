from flask import Flask
from flask import render_template

app = Flask('DxViewer')

@app.route("/")
def main():
	return render_template('index.html')
from flask import Flask, render_template, request, redirect, url_for, flash
import pandas as pd
#import pickle
import sqlite3 as sql
import numpy as np
from sklearn.externals import joblib
#import csv
import os
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
from werkzeug import secure_filename
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from io import BytesIO
import base64
from PyBioMed.PyGetMol.GetProtein import ReadFasta
from str import *
from Bio import SeqIO
from flask import send_file


UPLOAD_FOLDER = '/home/sanika/proj/'
ALLOWED_EXTENSIONS = set(['fasta'])


app = Flask(__name__)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'many random bytes'

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS



@app.route('/')
def main():
    return render_template('index.html')


@app.route('/about/')
def about():
    return render_template('about.html')


@app.route('/dataset/')
def datas():
   con = sql.connect("/home/sanika/proj/test.db")
   con.row_factory = sql.Row

   cur = con.cursor()
   cur.execute("select * from table3")

   rows = cur.fetchall();
   return render_template("datas.html",rows = rows)


@app.route('/predict/', methods=['GET', 'POST'])
def predict():



    if request.method == 'POST':
	
        seq = request.form['seq']
	if(seq):
		with open("random.fasta", "w") as fp:
			fp.write(seq)

		with open("random.fasta", "r") as fp_seq:
			ids = []
		    	seqs = []
	    	   	for seq_record in SeqIO.parse(fp_seq, 'fasta'):  # (generator)
				ids.append(seq_record.id)
				seqs.append(seq_record.seq)

		sequence_peptide = seqs[0]


		pepdesc = PeptideDescriptor('/home/sanika/proj/random.fasta', 'eisenberg')  # use Eisenberg consensus scale
		globdesc = GlobalDescriptor('/home/sanika/proj/random.fasta')

		# --------------- Peptide Descriptor (AA scales) Calculations ---------------
		pepdesc.calculate_global()  # calculate global Eisenberg hydrophobicity
		pepdesc.calculate_moment(append=True)  # calculate Eisenberg hydrophobic moment

		# load other AA scales
		pepdesc.load_scale('gravy')  # load GRAVY scale
		pepdesc.calculate_global(append=True)  # calculate global GRAVY hydrophobicity
		pepdesc.calculate_moment(append=True)  # calculate GRAVY hydrophobic moment
		pepdesc.load_scale('z3')  # load old Z scale
		pepdesc.calculate_autocorr(1, append=True)  # calculate global Z scale (=window1 autocorrelation)

		# --------------- Global Descriptor Calculations ---------------
		globdesc.length()  # sequence length
		globdesc.boman_index(append=True)  # Boman index
		globdesc.aromaticity(append=True)  # global aromaticity
		globdesc.aliphatic_index(append=True)  # aliphatic index
		globdesc.instability_index(append=True)  # instability index
		globdesc.calculate_charge(ph=7.4, amide=False, append=True)  # net charge
		globdesc.calculate_MW(amide=False, append=True)  # molecular weight

		f1=pepdesc.descriptor
		f2=globdesc.descriptor
		result=np.concatenate((f2,f1),axis=1)

		clf=joblib.load('ml_model.pkl')
		pred=clf.predict(result)
		proba=clf.predict_proba(result).tocoo()
		mc=pred.tocoo()
		out=mc.col
		res=[]
		labels=['antiviral','antibacterial','antifungal']
		values=proba.data
		plt.pie(values,labels=labels,autopct='%.0f%%',shadow=True, radius=0.5)
		plt.savefig('/home/sanika/proj/pie_chart.jpg')
	
		figfile = BytesIO()
		plt.savefig(figfile, format='png')
		figfile.seek(0)
		figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
		plt.close()
	    
		for i in range(len(out)):
		    if out[i]==0:
		       res.append("antiviral")
		    elif out[i]==1:
		       res.append("antibacterial")
		    else:
		       res.append("antifungal")
	

		return render_template('seq.html', result_seq = res, sequence_1 = seq, peptide_seq = sequence_peptide, result=figdata_png)
	else:
		return render_template('noseq.html')
    return render_template('predictor.html')

@app.route('/predict/result', methods=['GET', 'POST'])

def downloadimg():
    		return send_file('pie_chart.jpg',
                     mimetype='image/jpeg',
                     attachment_filename='pie_chart.jpg',
                     as_attachment=True)




@app.route('/predict/upload/', methods=['GET', 'POST'])

def upload():

    if request.method == 'POST':
        # This will be executed on POST request.
        upfile = request.files['file']
        if upfile and allowed_file(upfile.filename):
	    
            filename = secure_filename(upfile.filename)
            upfile.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	    #return render_template('upload.html')
	    #flash("File uploaded", "success")
            #with open("/home/sanika/proj/uploads/aa.fasta") as f:
    		#lines = f.readlines()
    		#lines = [l for l in lines if "ROW" in l]


    	    #with open("/home/sanika/proj/uploads/out.fasta", "w") as f1:
        	#f1.writelines(lines)

	    #f = open(filename)
	    #prot_seq = ReadFasta(f)

	    with open(filename) as fasta_file:  # Will close handle cleanly
    		identifiers = []
    		sequence = []
    	   	for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
        		identifiers.append(seq_record.id)
        		sequence.append(seq_record.seq)
	    
	    pepdesc = PeptideDescriptor(filename, 'eisenberg')  # use Eisenberg consensus scale
            globdesc = GlobalDescriptor(filename)

        # --------------- Peptide Descriptor (AA scales) Calculations ---------------
            pepdesc.calculate_global()  # calculate global Eisenberg hydrophobicity
            pepdesc.calculate_moment(append=True)  # calculate Eisenberg hydrophobic moment

        # load other AA scales
            pepdesc.load_scale('gravy')  # load GRAVY scale
            pepdesc.calculate_global(append=True)  # calculate global GRAVY hydrophobicity
            pepdesc.calculate_moment(append=True)  # calculate GRAVY hydrophobic moment
            pepdesc.load_scale('z3')  # load old Z scale
            pepdesc.calculate_autocorr(1, append=True)  # calculate global Z scale (=window1 autocorrelation)

        # --------------- Global Descriptor Calculations ---------------
            globdesc.length()  # sequence length
            globdesc.boman_index(append=True)  # Boman index
            globdesc.aromaticity(append=True)  # global aromaticity
            globdesc.aliphatic_index(append=True)  # aliphatic index
            globdesc.instability_index(append=True)  # instability index
            globdesc.calculate_charge(ph=7.4, amide=False, append=True)  # net charge
            globdesc.calculate_MW(amide=False, append=True)  # molecular weight

            f1=pepdesc.descriptor
            f2=globdesc.descriptor
            result=np.concatenate((f2,f1),axis=1)
	    rs=[]
	    for i in range(len(result)):
	 	prt=np.reshape(result[i],(-1,14))
		clf=joblib.load('ml_model.pkl')
		pred=clf.predict(prt)
		out=pred.toarray()
	#print(clf.predict_proba(result))
		proba=clf.predict_proba(prt).tocoo()
		mc=pred.tocoo()
		out=mc.col
		res=[]
		for i in range(len(out)):
	    		if out[i]==0:
	       			res.append("antiviral")
	    		elif out[i]==1:
	       			res.append("antibacterial")
	    		else: 
	       			res.append("antifungal")
		rs.append(res)
	    a = []
	    for i in range(len(rs)):
      		a.append('-'.join(rs[i]))

	

	    df = pd.DataFrame(data={"id": identifiers, "sequence": sequence, "activity":a},columns=['id','sequence','activity'])
	    df.to_csv("result.csv", sep=',',index=False)

	    countv=0
	    countb=0
	    countf=0
	    countvb=0
	    countbf=0
	    countvf=0
	    countvbf=0
	    for i in range(len(a)):
      		if a[i]=='antiviral':
           		countv=countv+1
      	   	elif a[i]=='antibacterial':
           		countb=countb+1
      		elif a[i]=='antifungal':
           		countf=countf+1
      		elif a[i]=='antiviral-antibacterial':
           		countvb=countvb+1
      		elif a[i]=='antibacterial-antifungal':
           		countbf=countbf+1
      		elif a[i]=='antiviral-antifungal':
           		countvf=countvf+1
      		else:
         		countvbf=countvbf+1
	    Dict={'antiviral':countv,'antibacterial':countb,'antifungal':countf,'antiviral-antibaterial':countvb,'antibacterial-antifungal':countbf,'antiviral-antifungal':countvf,'antiviral-antibacterial-antifungal':countvbf}

	    		
	    D={'AV-P':countv,'AB-P':countb,'AF-P':countf,'AV-P & AB-P':countvb,'AB-P & AF-P':countbf,'AV-P & AF-P':countvf,'AV-P, AB-P & AF-P':countvbf}


	    plt.bar(range(len(D)), D.values(), align='center',width=0.3)  # python 2.x
	    plt.xticks(range(len(D)), D.keys(),rotation='17',fontsize='small')  # in python 2.x
            plt.savefig('/home/sanika/proj/bar_graph.png')
	    figfile = BytesIO()
            plt.savefig(figfile, format='png')
            figfile.seek(0)
	    figdata_png = base64.b64encode(figfile.getvalue()).decode('ascii')
	    plt.close()
		
	   
	    


            

	  
            os.remove(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	    
	    #return render_template('seq.html', seq = rs)
	    return render_template('up.html', mimetype="text/csv",diction = Dict, bar=figdata_png)
	    
	    
		

    

            #flash("File uploaded: Thanks!", "success")
	else:
		error = "PLEASE CHECK THE FORMAT OF FILE TO UPLOAD"
		return render_template('upload.html', error=error)   

    # This will be executed on GET request.
    return render_template('predictor.html')



@app.route('/predict/download/')
def download():
    		return send_file('result.csv',
                     mimetype='text/csv',
                     attachment_filename='result.csv',
                     as_attachment=True)

	    

@app.route('/help/')
def help():
    return render_template('help.html')


if __name__ == '__main__':
    #clf = joblib.load('final_model.pkl')
    app.run(host='127.0.0.1' , port=1024, debug=True)

# damic100
notebooks for the DAMIC100 analysis:  

requirement: Python 3 


instructions:

- python3 -m venv venv 
- source venv/bin/activate 
- pip install -r requirements.txt
- python -m ipykernel install --user --name=venv
- jupyter notebook efficiency.ipynb
- choose the venv as kernel
- Download the pickle folder in Lyon server:  
/sps/hep/damic/gaior/efficiency/20200515/pkl/  
- change the 'datapath' variable to the folder where you downloaded the /pkl/ folder



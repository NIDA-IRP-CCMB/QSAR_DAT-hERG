The code in this repository is the chemoinformatics part of our project in developing an infrastrature for counter screening small-molecule ligands for a selected target and against hERG. 
For the moment, we focus on the ligands for the dopamine transporter (DAT).  

The regression models are based on the methods found in Wacker and Noskov
[Performance of machine learning algorithms for qualitative and quantitative
prediction drug blockade of hERG1 channel](https://doi.org/10.1016/j.comtox.2017.05.001).
The classifier models were initially derived from the methods of
Siramshetty, et al [The Catch-22 of Predicting hERG Blockade Using Publicly Accessible Bioactivity Data](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00150). 

The counter or synergistic screening platform can be easily adapted to identify novel compounds with desired selectivity for other targets, such as the ligands antagonizing both DAT and sigma1 receptor. 

***Our improvements***
 
Due to the stochastic element in the prediction, for regression models, we create multiple models to be used in the prediction and take the average. The classifier models were modified to use the RDKit topological
descritors rather than fingerprint descriptors derived from the Morgan Topological fingerprints.



***Runtime environment***

These models are coded in python 3.6.4.  They were developed using RDKit version
2018.3.1, scikit-learn 0.21.1, MolVS 0.1.1, and xgboost .  Xgboost,
Scikit-learn, RDKit and prerequisites were installed from the anaconda
repository (using the channel rdkit for rdkit and its prerequisites).  MolVS
was installed under anaconda using pip.  If using a recent version of anaconda
python to build and run the models, it is probably necessary to create a virtual
environment based on python 3.6.  At present, trying to install RDKit into a
python 3.7 environment does not appear to be supported and causes anaconda to
regress the entire installation to python 2.7.  Scikit-learn is installed by
default with the full anaconda distribution (not miniconda). The commands to
install xgboost, RDKit and MolVS into an active python 3.6 virtual environment
are:
```
conda install py-xgboost
conda install -c rdkit rdkit 
#http://rdkit.org/docs/Install.html#cross-platform-under-anaconda-python-fastest-install
#just tried using conda-forge
#conda install -c conda-forge rdkit
pip install MolVS
pip install -U scikit-learn
```

***About Descriptor Generation***

The 2D descriptors that are available in rdkit are referenced in Descriptors.descList.  A code snippet that will list
these descriptors at an interactive python prompt is:
```
from rdkit.Chem import Descriptors
for foo in Descriptors.descList:
    print(foo[0])
```

***Strucure of the repository***

    .
    ├── core            <-- major scripts 
    ├── example         <-- example of running filters/buildmodels/prediction
    └── README.md

***unittest***

    module load python
    # unittest for run_filters
    python -m unittest -v core/unittest_filters.py
    # use these following commands
    python -m unittest -v core/unittest_buildmodel.py

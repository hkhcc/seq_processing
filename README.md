# seq_processing
This repository contains a collection of Python-based tools developed by Chi-Chun Ho, initially intended for internal use at Division of Chemical Pathology, Department of Clinical Pathology, Pamela Youde Nethersole Eastern Hospital (PYNEH), Hong Kong, China.
The scripts are provided “as is”, without being actively supported or maintained.
1.	Installation
  * The scripts are written in Python and depend on a number of modules. The easiest way to install the required modules is by installing Anaconda.
  * Anaconda may be downloaded from https://www.anaconda.com/distribution/
  * During installation, the Windows user is suggested to select installation of “Just Me” (as recommended), “Add Anaconda to my PATH environment variable” (despite not recommended) and “Register Anaconda as my default Python 3.x”.
  * Most scripts require a working Internet connection to function properly. For example, the core program autoprimer.py downloads transcript information from the Ensembl REST server, and other programs (e.g. web_search.py) may require other web services. 
2.	Running
  * Check that Python is properly installed. In the Windows command processor (cmd.exe), type
```
python
```
and the Anaconda Python interpreter should be invoked:
```
Python 3.7.0 (default, Aug 14 2018, 19:12:50) [MSC v.1900 32 bit (Intel)] :: Anaconda, Inc. on win32
Type "help", "copyright", "credits" or "license" for more information.
>>>
```


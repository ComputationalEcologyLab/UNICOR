======
README
======

------------------  
UNICOR 3.0 release
------------------

Welcome to the UNICOR v 3.0 release!  UNICOR is a network simulation model designed to help researchers and managers identify vulnerable species and landscapes and evaluate the relative merits of conservation, restoration, or assisted mitigation.  UNICOR is specifically designed to enable quantification of future changes to habitat fragmentation and connectivity by comparing predicted landscape changes stemming from climate regime shifts, urban development, and their interactions.  This release includes installation instructions, version notes, some examples, and technical documentation.

Program Contributors: Erin Landguth, Brian Hand, Joe Glassy, Casey Day, Patrick Jantz
Link: http://github.com/computationalecologylab/UNICOR
Version: 3.0
Python: 3.8.3
Release Date: 2022.01.10
README Update: 2022.01.10 (ell)

--------
Contents
--------

Included in this release are the following:

README.pdf 			this file

UNICOR Source:
UNICOR.py - Python driver code
UNICORFun.py - Python main module/function code
UNICORMaps.py - Python display module code
UNICOROutputs.py - Python file output module code
unicor_script.py - Python files for running batch files
  
Additional Packages:
RipMgr - Installation file for input package by Glassy (2011)
priodict.py - Priority queue by Eppstein (2002)

UNICOR Example Files: 
small_test.rsg - Sample resistance grid
small_test.xy - Sample point locations
small_test.rip - Sample run input parameters			

---------------------------------------
Requirements and Pre-requisite Software
---------------------------------------

1. Baseline Requirements.  UNICOR requires the Python3 interpreter and the NumPy and SciPy package.  Several optional Python module packages, if enabled, facilitate additional UNICOR functionality. Remember that Python modules usually require  particular Python interpreters, so be sure the version ID for any external Python module or package (e.g. NumPy or others) matches the version of your Python interpreter (normally v3.x.x).

-------------------	
UNICOR Installation
-------------------

1. Unpack the UNICOR Archive.  Navigate to the directory on your PC where you wish to install UNICOR, and unpack the supplied zip archive file using a free archive tool like 7Zip (7z.exe), Pkunzip, Unzip, or an equivalent.  Seven-Zip (7Z.exe) is highly recommended since it can handle all common formats on Windows, MAC OS X and Linux. On Windows, it is best to setup a project specific modeling subdirectory to perform your modeling outside of any folder that has spaces in its name (like "My Documents").

2. Install UNICOR.  Next, install the UNICOR software itself by unpacking the zip archive supplied. At this point you should be able to execute the supplied test inputs (eh8k.txt and MT10points.txt with user_input.ip).

------------------
Example UNICOR Run
------------------

The example run is for 10 points representing individuals on a Euclidean distance resistance surface. To run the following example, follow these steps:

1. Double check the UNICOR source, UNICOR additional packages, and UNICOR example files are in the same directory.

2. The included .rip file specifies the parameters that can be changed and used in the sample UNICOR run.  Open example.rip in your editor of choice (e.g., notepad or wordpad).

3. This file is the stanza format following RipMgr documentation.  All '#' signs are comments followed by variable names with a tab to the parameter specified.  The parameter can be changed for running UNICOR, but downloaded parameters will run as is.  See table 1 for more details on each parameter along with its dependency.

4. Start the program with a graphical interface or at the command line: For example, if you downloaded Python 3 from www.python.org, then you are provided with a graphical interface, IDLE.  In Windows you can find IDLE from your Start menu > All Programs > Python 3.x > IDLE (Python GUI).  Alternatively, if you use Python from the command line, then open a terminal window and change your shell directory to the UNICOR home directory.

5. Run the program: There are a number of ways to run this program.  For example, if you are using a command shell you can run the program by typing “python UNICOR.py example.rip”.

6. Check for successful simulation run completion: The program will provide a log file in your UNICOR home directory.  Once completed, output files will be created in UNICOR home directory.

Happy Simulations!

Erin.

Contact Information:

Erin Landguth
Computational Ecology Laboratory
The University of Montana
32 Campus Drive
Missoula MT, 59812-1002
(406) 210-9332
Erin.landguth@mso.umt.edu
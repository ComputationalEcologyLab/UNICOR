======
README
======

------------------  
UNICOR 1.0 release
------------------

Welcome to the UNICOR v 2.0 release!  UNICOR is a network simulation model designed to help researchers and managers identify vulnerable species and landscapes and evaluate the relative merits of conservation, restoration, or assisted mitigation.  UNICOR is specifically designed to enable quantification of future changes to habitat fragmentation and connectivity by comparing predicted landscape changes stemming from climate regime shifts, urban development, and their interactions.  This release includes installation instructions, version notes, some examples, and technical documentation.

Program Contributors: Erin Landguth, Brian Hand, Joe Glassy
Link: http://cel.dbs.umt.edu/software/UNICOR/
Version: 2.0
Python: 2.7.2
Release Date: 2012.03.01
README Update: 2012.03.12 (ell)

--------
Contents
--------

Included in this release are the following:

README_UNICOR.txt - this file

UNICOR Source Folder.
UNICOR.py - Python driver code
UNICORFun.py - Python main module/function code
UNICORMaps.py - Python display module code
UNICOROutputs.py - Python file output module code

UNICOR Utilities Folder.
  
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

1. Baseline Requirements.  UNICOR requires the Python2.6.x interpreter and the NumPy and SciPy package.  Several optional Python module packages, if enabled, facilitate additional UNICOR functionality. Remember that Python modules usually require  particular Python interpreters, so be sure the version ID for any external Python module or package (e.g. NumPy or others) matches the version of your Python interpreter (normally v2.6.x).

-------------------	
UNICOR Installation
-------------------

1. Navigate to the directory on your PC where you wish to install UNICOR, and unpack the supplied zip archive file using a free archive tool like 7Zip (7z.exe), Pkunzip, Unzip, or an equivalent.  Seven-Zip (7Z.exe) is highly recommended since it can handle all common formats on Windows, MAC OS X and Linux. On Windows, it is best to setup a project specific modeling subdirectory to perform your modeling outside of any folder that has spaces in its name (like "My Documents").

2. Install UNICOR.  Next, install the UNICOR software itself by unpacking the zip archive supplied. At this point you should be able to execute the supplied test inputs (small_test.rsg and small_test.xy with small_test.ip). Note with v2.0 you will install 2 folders (utilities and unicor).

------------------
Example UNICOR Run
------------------

The example run is for 10 points representing individuals on an example landscpae resistance surface. To run the following example, follow these steps:

1. Double check the 2 folders (utilities and unicor) are installed in the directory location of your choice.  

2. The included .rip file specifies the parameters that can be changed and used in the sample UNICOR run.  Open small_test.rip in your editor of choice (e.g., notepad or wordpad).

3. This file is the stanza format following RipMgr documentation.  All '#' signs are comments followed by variable names with a tab to the parameter specified.  The parameter can be changed for running UNICOR, but downloaded parameters will run as is.  See table 1 in the complete user manual for more details on each parameter along with its dependency.

4. Start the program at the command line: If you use Python from the command line, then open a terminal window and change your shell directory to the UNICOR home directory. For example in Windows you can go to the Start Menu -> and open cmd.exe by searching for it or locating it in Accessories. After ">" type cd C:/LOCATIONOFYOURINSTALLATION.

5. Run the program: There are a number of ways to run this program.  If you are using a command shell you can run the program by typing “python UNICOR.py example.rip”.

6. Check for successful simulation run completion: The program will provide a log file in your UNICOR home directory.  Once completed, output files will be created in UNICOR home directory.

Happy Simulations!

Erin.

Contact Information:

Erin Landguth
Computational Ecology Laboratory
Division of Biological Sciences
The University of Montana
32 Campus Drive
Missoula MT, 59812-1002
(406) 210-9332
Erin.landguth@mso.umt.edu
<cel.dbs.umt.edu>


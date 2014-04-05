   --  Description:  --

Moltemplate is a cross-platform text-based molecule builder for LAMMPS. 

   --  Typical usage:  --

moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd]  system.lt

   --  Web page:  --

Documentation, examples, and supporting code can be downloaded at:

    http://www.moltemplate.org
The most up-to-date version of moltemplate can be downloaded here.
(After download, you can unpack the archive using:
    tar xzf moltemplate_2012-3-31.tar.gz
The date will vary from version to version.)

----------------------------------------------------
----------  INSTALLATION INSTRUCTIONS:  ------------
----------------------------------------------------

This directory should contain two folders:
    src/       <-- location of all python and bash scripts
    common/    <-- location of shared force fields and molecules 

The ``moltemplate.sh'' script and the python scripts that it invokes are 
located in the ``src/'' subdirectory.  You should update your PATH environment 
variable to include this directory.  

If you do not know what a PATH environment variable is, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)

It is also a good idea to set your MOLTEMPLATE_PATH environment variable to 
point to the ``common/'' subdirectory.
(Force fields and commonly used molecules will eventually be located here.)

   --  Installation example ---

Suppose the directory with this README.TXT file is located at ~/moltemplate.

If you use the bash shell, typically you would edit your 
~/.profile, ~/.bash_profile or ~/.bashrc files to contain the following lines:
    export PATH="$PATH:$HOME/moltemplate/src"
    export MOLTEMPLATE_PATH="$HOME/moltemplate/common"
If you use the tcsh shell, typically you would edit your 
~/.login, ~/.cshrc, or ~/.tcshrc files to contain the following lines:
    setenv PATH "$PATH:$HOME/moltemplate/src"
    setenv MOLTEMPLATE_PATH "$HOME/moltemplate/common"

   --  Requirements:  --

Moltemplate requires the Bourne-shell, and a recent version of python 
(2.7, 3.0 or higher), and can run on OS X, linux, or windows (if a 
suitable shell environment has been installed).

   --  License:  --

Moltemplate is available under the terms of the open-source 3-clause BSD 
license.  (See LICENSE.TXT.)

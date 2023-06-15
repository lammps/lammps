
Moltemplate
===========

##  Description

This folder used to contain a distribution of Moltemplate, a general purpose,
cross-platform, text-based molecule and topology builder for LAMMPS.
Moltemplate was originally conceived for building custom coarse-grained
molecular models, but it has since been generalized for all-atom simulations
as well.  It currently supports the OPLS, COMPASS, AMBER(GAFF,GAFF2),
MARTINI, SPICA(SDK), LOPLS(2015), and TraPPE(1998) force fields, and includes
(New force fields and examples are added continually through user
contributions).

Moltemplate is now distributed via the pip python package manager,
or can be downloaded from http://www.moltemplate.org/download.html
or https://github.com/jewettaij/moltemplate/releases
The most up-to-date version is usually available through GitHub

## Typical usage

    moltemplate.sh [-atomstyle style] [-pdb/-xyz coord_file] [-vmd] system.lt

## Web page

Documentation, examples, and supporting code can be downloaded at:

http://www.moltemplate.org

## Requirements

Moltemplate requires the Bourne-shell, and a recent version of python
(2.7, 3.0 or higher), and can run on OS X, linux, or windows. (...if a
suitable shell environment has been installed.  See below.)


## INSTALLATION INSTRUCTIONS

## Download / install via pip
If you have a working installation of "pip" and access to the internet,
you can download and install it into your home directory ($HOME/.local/bin
and the corresponding python package folders) simply by:

   pip install moltemplate --user

If instead, you want a system-wide installation, you need superuser access,
e.g. via:

   sudo pip install moltemplate

Updates to this distribution method are less frequent, than others, so if
you need a more recent version, you can download it as a .tar.gz or .zip
archive from the moltemplate home page or GitHub (see link above). After
downloading an archive and unpacking it, you should have 3 folders.

    moltemplate/                  <-- source code and force fields
    doc/                          <-- the moltemplate reference manual
    examples/                     <-- examples built with moltemplate

From here on, you can use pip to install this downloaded version as well
by issuing the following command from within the top-level directory of
the unpacked source archive:

    pip install . --user

If you want to install into a system folder, then you need to run pip with
superuser privileges. e.g. with:

    sudo pip install .

Make sure that your default pip install bin directory is in your PATH.
This is usually something like ~/.local/bin/ or ~/anaconda3/bin/.
If you have installed anaconda, this will be done for you automatically.
As you can see, this is mostly identical to the direct install with pip.

Later, you can uninstall moltemplate using:

    pip uninstall moltemplate

If you continue to run into difficulty, try installing moltemplate into
a temporary virtual environment by installing virtualenv, downloading
moltemplate (to "~/moltemplate" in the example below), and running
these commands:

    cd ~/moltemplate
    virtualenv venv
    source venv/bin/activate
    pip install .
    #(now do something useful with moltemplate...)

(You will have to "run source ~/moltemplate/venv/bin/activate" beforehand
every time you want to run moltemplate.
The *virtualenv* tool is explained in detail here:

http://docs.python-guide.org/en/latest/dev/virtualenvs/

If all this fails, then try installing moltemplate by manually updating
your \$PATH environment variable.  Instructions for doing that are included below.

## Manual installation:

Alternatively, after downloading and unpacking the moltemplate source archive,
you can edit your $PATH environment variable manually to include the subdirectory
where the "moltemplate.sh" script is located, as well as the subdirectory where
most of the python scripts are located.
Suppose the top-level directory of the unpacked archive "moltemplate"
and is located in your home directory:

If you use the bash shell, typically you would edit your
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files
to contain the following lines:

    export PATH="$PATH:$HOME/moltemplate/moltemplate"
    export PATH="$PATH:$HOME/moltemplate/moltemplate/scripts"

If you use the tcsh shell, typically you would edit your
`~/.login`, `~/.cshrc`, or `~/.tcshrc` files to contain the following lines:

    setenv PATH "$PATH:$HOME/moltemplate/moltemplate"
    setenv PATH "$PATH:$HOME/moltemplate/moltemplate/scripts"

After making these changes, you may need to start a new terminal (shell) for the changes to take effect.  If you do not know what a `PATH` environment variable is and are curious, read:
    http://www.linfo.org/path_env_var.html
(I receive this question often.)


### WINDOWS installation suggestions

You can install both moltemplate and LAMMPS in windows, but you will
first need to install the BASH shell environment on your computer.
If you are using Windows 10 or later, try installing the "Windows Subsystem for Linux (WSL)"

https://solarianprogrammer.com/2017/04/15/install-wsl-windows-subsystem-for-linux/
https://msdn.microsoft.com/en-us/commandline/wsl/faq

If you are using an older version of windows, try following the tutorial written by Yanqing Fu instead:

https://sourceforge.net/p/lammps/mailman/message/32599824/

To use LAMMPS and moltemplate, You will also need to install (and learn how to use) a text editor.  (Word, Wordpad, and Notepad will not work.)  Popular free text editors which you can safely install and run from within the WSL terminal include: **nano**, **ne**, **emacs**, **vim**, and **jove**.  (Unfortunately, as of 2017-5-17, [graphical unix-friendly text editors such as Atom, VSCode, Notepad++, and sublime won't work with WSL, and may cause file system corruption.  Avoid these editors for now.](https://www.reddit.com/r/bashonubuntuonwindows/comments/6bu1d1/since_we_shouldnt_edit_files_stored_in_wsl_with/))


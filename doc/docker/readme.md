**To use Docker to build the documentation:**

*One-time setup*
- Get and install a version of Docker appropriate for your host system. 
Get started here: https://docs.docker.com/get-started/
- Open a command prompt on your host system and make this folder your 
current working directory.
- Build a container with the commandline: 
    docker build --tag=lammpsdoc .
will build or update a container with the name "lammpsdoc". 
Yes, there is a trailing . on that line.
- To test/explore your new container use the commandline:
    docker run -it
which will run the container and open up an interactive bash shell you can use to 
look around.

*Building the documentation*
- Unless you love vim and all tools commandline, don't use your new container for 
git or editing files. 
- Instead, use the container only for the build and do all your git and editing 
work on your host machine.
- When you are ready to build, run the container with your host lammps folder
as a volume:
    docker run -it -v <local repo folder>:/<container foldername> <container tag>
example from my Windows box:
    docker run -it -v C:/openSource/lammps:/masterlammps lammpsdoc
  - C:/openSource/lammps is my local lammps git repo
  - /masterlammps is the work folder in the container for the build
  - lammpsdoc is the tag used in the build step above

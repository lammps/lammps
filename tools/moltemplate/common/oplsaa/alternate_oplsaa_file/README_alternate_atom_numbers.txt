MOST USERS CAN IGNORE THIS DIRECTORY.

This is the original version of oplsaa.txt released with 
moltemplate version 1.20.  Later on I decided to direct 
users to the (hopefully current) version of the "oplsaa.prm" 
file distributed with TINKER's force-field-explorer program.

http://dasher.wustl.edu/ffe/distribution/params/oplsaa.prm
http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm

The original oplsaa.txt from moltemplate version 1.20 file is located here.
Feel free to pick whichever one works well for your application, but 
if you publish a paper using these parameters please cite TINKER,
since that appears to be the source of all of these files.
I think these are the relevant citations:

Ponder, J. W., & Richards, F. M. (1987). "An efficient newton‐like method for molecular mechanics energy minimization of large molecules. Journal of Computational Chemistry", 8(7), 1016-1024.

Ponder, J. W, (2004) "TINKER: Software tools for molecular design", http://dasher.wustl.edu/tinker/


Here is a discussion of a difference between these two files:

>> On Fri, Apr 11, 2014 at 6:44 PM, Andrew Jewett <jewett.aij@gmail.com> wrote:
>> Jason
>> Lambert graciously wrote the "oplsaa_moltemplate.py" script and
>> provided the "oplsaa.txt" file which it reads as an input file (after
>> some minor manual editing on the part of the user).  I did some
>> googling, and I noticed that this file is very similar to the
>> "oplsaa.prm" file located at this web site:
>>
>> http://dasher.wustl.edu/ffe/distribution/params/oplsaa.prm
>>
>> I suspect it uses the same parameter conventions.
>>
>> However I just realized that the two files are not identical.  (I am
>> emailing Jason now to ask about the difference.)
>>
>> Where did you get the data for the "oplsaa.txt" file?  More
>> specifically, where do the atom numbers come from?

>On Tue, Apr 15, 2014 at 6:06 PM, jason Lambert <jlamber9@gmail.com> wrote:
>To be honest I do not know the original source of the OPLSA force
> field file. The file was provided by a graduate student from a group
> I worked with and I did not question it other than spot checking with
> values that were published.  I am surprised that the files are different.
> Anyways, I am not aware of any logic to the observed switches in
> the parameters. I did a diff -y on the two files and the difference is
> primarily atom switching. The code works with either one though
> so it is not a problem.  We can provide a link to the prm file as well.
> I know I added some dihedrals explicitly but the ones I added were
> already encompassed by the dihedrals containing wild cards.

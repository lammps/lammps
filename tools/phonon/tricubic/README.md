# libtricubic

This folder contains a slightly refactored version of version 1.0 of
the _libtricubic_ library, developed by François Lekien in 2004.

The following paper explaining the method was published in 2005:

> Lekien, F., & Marsden, J. (2005). _Tricubic interpolation in three
> dimensions._ International Journal for Numerical Methods in Engineering,
> 63(3), 455–471. doi:10.1002/nme.1296

Some additional notes and the full matrix can be found in the technical notes:

> Lekien, F., Coulliette, C., & Marsden, J. (2004). _Tricubic Engine - Technical
> Notes and Full Matrix.

The main article refers to the author's website
(http://gyre.cds.caltech.edu/pub/software/tricubic/) to download the code and
documentation. Unfortunately, this website no longer exists. Even the
archive.org snapshot is useless; [A single snapshot from December 3rd 2009 is
available](https://web.archive.org/web/20091203115835/http://gyre.cds.caltech.edu/pub/software/tricubic)
which seems like an FTP listing. No files are accessible.

The source code was obtained from https://github.com/nbigaouette/libtricubic/
Only the sources for the library were retained. No functional changes were
made, but some common programming conventions were applied, source files
merged, and build support for CMake added.


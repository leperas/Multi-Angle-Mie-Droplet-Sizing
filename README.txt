Welcome to the Multi-Angle Mie Sizing project.

These Matlab routines are an important part of an optical droplet measurement system based on planar multi-angle Mie scattering.  Sizing information consists of a mean droplet diameter and droplet distribution estimates for every individual point within a planar area of interest.  The planar method makes possible the fast acquisition of data within a large field of interest, and uses relatively inexpensive instrumentation. 

Currently, the method has demonstrated the ability to measure water droplets from a typical simplex spray nozzle, across the range of 5-50$\mu$m within +/-10\% of known values, and in addition return an estimate of the shape and width of the size distribution at each location within the planar region of interest.  Measurements were successfully completed for three separate flow rate spray nozzles; 1-gallon-per-hour, 2.5-gallon-per-hour, and 4.5-gallon-per-hour.

Additionally the limits of the technique have been explored with simulated data.  Conclusions from these exercises show that the multi-angle planar Mie scatter method is capable of measuring droplet distribution characteristics and means within a nominal range of 0.3$\mu$m up to 150$\mu$m.

A quick disclimer: I am not a programmer; suggestions and improvements are welcome.
----------------
VERSIONS and BRANCHES

-- master: This is the repository you want; there may be many other branches but Master should contain the latest working version of everything, EXCEPT you must download the scattering coeffcient database files separately and unzip them into the ./multi_angle_mie_sizing/database_files directory.  I will keep a copy of any databases I generate updated on my Google+ Drive.  These are NOT part of the Repo.  Installation instructions are found in README.install.

-- Version 1 - Tag: dissertation

The master repository TAG "dissertation" should take you back to the original commit of the dissertation version of all the routines, if for some reason you would want that (you don't).  

-- The compressed download (old)

There are no updates, bug fixes, or new features...but perhaps the simplest way to get started is to download the tar.bz2 package.  All the files referenced within the dissertation are contained within multi-angle-mie-sizing_v1.tar.bz2, including the appropriate cross section database files.  This v1 package will not change or be updated in order to maintain consistancy with the original dissertation document.  This package is available on my Google+ Drive at the following URL:

https://docs.google.com/open?id=0By0Xr9D8LwH2aGZ5YVpabi1ra2M

To install, put this file into it's own directory and run:

tar -jxf multi-angle-mie-sizing_v1.tar.bz2

Additional installation instructions are found within the README.install included in the multi-angle-mie-sizing_v1.tar.bz2 package.  Do NOT use the Repo install instructions with the old v1 package, use the install instructions included in the zipped package!!!

DOCUMENTATION

The MASTER branch in the repository now has limited documentation.  Documentation is still mostly my dissertation (available online from http://www.lib.vt.edu/).  However, there is now also limited documentation in the form of a User Manual PDF in the documentation directory.  Currently this consists only of a complete HOWTO step-by-step example in the Appendix.  There is no in-depth explanation of options or use, but the example will work with the most recent repository and most options are documented in the code itself, so all is not lost.  The long term goal is for this documentation to gradually improve and eventually replace the disertation.

If you are unsure about what to do, consider starting out by downloading the Version 1 tar.bz2 package and installing it as suggested above.  Get my dissertation from http://www.lib.vt.edu/.  Read it and go through the example in the appendix.  At that point you will have the basic idea about things and can move on to the latest version and features.

If you grab the Repo and want to really modify things and make improvements, feel free to contact me so I can explain what is going on, and we can work together constructively. 

Thanks for your iterest.
Stephen LePera
2013

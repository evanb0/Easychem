# VERSION 2 OF CODE

Contains everything I have coded up to now (Version 2).

The code has changed quite a lot.

Inside each file I have tried to implement comments to explain what's going on, so it will probably be worth reading those if you want to take a look. Just to summarise what should be possible now, and will probably be worth testing:

- App should take SDF, XYZ, .com and gjf files (I have only tested SDF and XYZ, I have not tested any of the gaussian functionality at this point and it will probably just result in an error)
- Interface of software changes from blue to red depending on whether you choose ORCA or Gaussian (useless for what we want to do, just thought it'd be fun to add)
- I tested two runs, one with an SDF file and one with an XYZ file. Both ended up coming back working.
- Right now, the OPT keyword is forced by default for Job types, I mainly did this just so I could test multi-job running and it did indeed work, I was able to automate the creation of an OPT/FREQ calc.

I don't really know what to test, but I imagine there will be quite a few errors, and I think it'd be worth you taking a look at if you have the time. In terms of what to do next, I think I will wait on your feedback first before I decide what else to add, it definietly is not finished but we are making good progress.

On Monday I plan to run through all the GAUSSIAN tutorials on the research group github page anyway, I've never used it before so it'll probabaly be worth getting some practice with what it can do/how it creates its imput files anyway, since so far i just guessed from the documentation of gaussain.

Hopefully you are satisfied with the changes that have been made, however anything you think needs to be improved on/is just bad please let me know in as much detail as possible.



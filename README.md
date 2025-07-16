# summerproj
Updated solvation models for ORCA and GAUSSIAN.
Orca documentation was pretty clear, I have included every method other than the cluster methods (Idk if we are going to run with these so I haven't included them yet).
Every solvent option works for ORCA every solvation method we have, I didn't find a definite list for GAUSSIAN solvents so I think it might be worth just keeping them the same, and then if the user generates an imput file and they want a specific solvent they might have to modify it themselves.

CHARGE/MULTIPLICITY SHOULD No longer be hardcoded, I tested it with modified values and it worked, so I think this is good now.

SLURM Button generates 2 things depending on the software chosen:
if software == ORCA -> a submitorca.txt file is made with the necessary things for submission. NOTE THAT THIS AS WELL AS THE INPUT FILE e.g ethane.inp need to be copied (Via scp or other means) to csf3 before it can be submitted properly.

if software == GAUSSIAN -> Theoretically, the same should be the case. A submitgaussian.txt is made with then necessary things for submission. IDK IF The submission file needs the GAUSSIAN file in .com or .inp, however if it is .inp this currently needs to be done manually mv file_name.com file_name.inp. I did not manage to get the GAUSSIAN submission script to work on my CSF3, however I think this is just due to me not having access to GAUSSIAN and not actually an error with the script, If nothing else, I think it'd be worth for you to test this!

I can't think of any other changes I've implemented, but if you can please send over an .SDF file with several conformers, because I don't really understand the idea behind it. May need to talk in detail about this more when im actually in bc idk if its just me being dumb with chemistry or what.

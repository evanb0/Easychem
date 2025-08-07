# summerproj

MOLECULARVIEWER Still does not work on windows (and i still dont know why)

Completely removed SMILES FUNCTIONALITY where initially we made new conformers from a file.
- When a .sdf file is inserted, a 3D Molecular viewer opens at the bottom of the screen, allowing each molecule within that file to be selected and viewed accordingly.
- Option to export 10 lowest energy molecules within the file (this does work for the 6k molecule file I tried it, although it does make the app freeze)
- Changed the UI a little bit (added a zoom button to 3D viewer, changed colour scheme)
- Currently the molecular viewer stuff only works for .sdf files, a lil buggy with .xyz so haven't implemented it yet
- All other features haven't really been touched.

>[!Important] 
>Possible things to do:

- [ ] Ability to upload multiple files at once and treat them all at once
- [ ] Amend the methodology section for both Orca and Gaussian to include more (and more relevant) methods
- [ ] Find a way to work with XYZ files (and many of them at once)
- [ ] Add unit tests (I hear every software endeavour has them idk we're both amateurs here)
- [ ] *For dealing with multiple files later on* there needs to be a way to automatically assign mult/charge and not impose it to all mols at once

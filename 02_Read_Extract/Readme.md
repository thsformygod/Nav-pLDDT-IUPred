This is description of basic scripts used to read, extract features such as sequence, pLDDT and IUPred scores.

Read information from PDB.py:
	This is used to extract many features, such as PDB ids, entries, reference ID, being and end of reference, reference entry ids and so on.
	
Read sequence from PDB.py:
	This is basically to check sequence from PDB chains, to check if there is missing CA and return.
	
Read pLDDT.py:
	This is to read sequence and pLDDT from given alphafold .cif. these .cif could be downloaded from alphafoldDB or calculated on local.
	
Read IUPred.py:
	This is to send sequence for installed iupred.py, and return scores. we need to have access to iupred.py (download and install from  https://iupred3.elte.hu/download_new)


Examples of inputs are either given in scripts, of included in inputs/ folder.
Feel free to contact us if anything wrong.

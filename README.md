# SuCOS select

This script attempts to compare docked conformers with fragment hits using SuCOS

See SuCOS explained in ChemRxiv [manuscript](https://chemrxiv.org/articles/SuCOS_is_Better_than_RMSD_for_Evaluating_Fragment_Elaboration_and_Docking_Poses/8100203) 
Code edited from original SuCOS [github](https://github.com/susanhleung/SuCOS)

This script: 

- Gets SuCOS score of docked conformer of the Moonshot designs and all the Mpro fragment hits
- Checks if highest SuCOS scoring fragment is in list of inspiration fragments 
- Computes average of all the SuCOS scores
 

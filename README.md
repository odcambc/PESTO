# PESTO: Protein Evolutionary Similarity and Topological Organization

## Overview
Some scripts to fetch and align protein structures from the AlphaFold database.

## Instructions

1. Find targets with Foldseek
   * Submit a target PDB to [Foldseek](https://search.foldseek.com/)
   * Download hit tables (in M8 format)
     * You'll get a series of files giving hits across different databases
2. Download model PDBs from AlphaFold database
   * Run `fetch_af_pdb.py` with hit table to download models to a folder
     * `python fetch_af_pdb.py -i alis_afdb50.m8 -o Output/folder` 
       * Get more options using `python fetch_af_pdb.py -h`
   * Alternatively, supply a list of IDs to download
     *  Extract them from the m8 output using awk:
        *  `awk '{print $2}' alis_afdb50.m8 > alis_afdb50_ids.txt`
     *  Run `fetch_af_pdb.py` with output list to download models to a folder
        * `python fetch_af_pdb.py -i alis_afdb50.m8 -o Output/folder` 
3. Run `align_all.py` to perform pairwise alignment and TM-score calculation for all downloaded files
   * `python align_all.py -i Output/folder -o afdb50_tm_scores.csv -j 64`
   * This will generate an output csv in upper triangular format containing the TM-scores between all pairs of structures. For large datasets, this will take a while. Set threads with `-j` flag to speed this up.
   * To also generate a set of PDB files aligned to a reference structure, supply the reference PDB with the `-r` flag. You can also specify a directory to save the aligned PDBs to with the `-d` flag: by default it will save them to a folder named `aligned` in the current directory.
     * `python align_all.py -i Output/folder -o afdb50_tm_scores.csv -j 64 -r AF-O25823-F1-model_v4.pdb`
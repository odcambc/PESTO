# Sequence-Structure Similarity Network

## Overview
Some scripts to fetch and align protein structures from the AlphaFold database.

## Instructions

1. Find targets with Foldseek
   * Submit a target PDB to [Foldseek](https://search.foldseek.com/)
   * Download hit tables (in M8 format)
     * You'll get a series of files giving hits across different databases
2. Extract AlphaFold IDs from hit tables
   * Using `awk`: run the following
     *  `awk '{print $2}' alis_afdb50.m8 > alis_afdb50_ids.txt`
        *  Change inpur and output file names as needed
3. Run `fetch_af_pdb.py` to download PDB files from AlphaFold database
   * `python fetch_af_pdb.py -i alis_afdb50_ids.txt -o Output/folder` 
4. Run `align_all.py` to perform pairwise alignment and TM-score calculation for all downloaded files
   * `python align_all.py -i Output/folder -o afdb50_tm_scores.csv -j 64`
   * This will generate an output csv in upper triangular format containing the TM-scores between all pairs of structures. For large datasets, this will take a while. Set threads with `-j` flag.
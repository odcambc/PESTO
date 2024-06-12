# PESTO: Protein Evolutionary Similarity and Topological Organization

## Overview
Some scripts to fetch and align protein structures from the AlphaFold database.

## Requirements
   * [USAlign standalone](https://zhanggroup.org/US-align/)
   * python libraries (install with `pip install`)
     * requests
     * MDAnalysis
     * pandas
     * pdb-tools

## Instructions

1. Install requirements
2. Find targets with Foldseek
   * Submit a target PDB to [Foldseek](https://search.foldseek.com/)
   * Download hit tables (in M8 format)
3. Download model PDBs from AlphaFold database
   * Run `fetch_af_pdb.py` with hit table to download models to a folder
     * `python fetch_af_pdb.py -i alis_afdb50.m8 -o Output/folder` 
     * Specify minimum tlen or qlen using `--min_tlen` and `--min_qlen` flags 
       * See description of foldseek M8 fields here: https://github.com/steineggerlab/foldseek/issues/59
     * Get more options using `python fetch_af_pdb.py -h`
   * Alternatively, supply a list of IDs to download
     *  Extract them from the m8 output using awk:
        *  `awk '{print $2}' alis_afdb50.m8 > alis_afdb50_ids.txt`
     *  Run `fetch_af_pdb.py` with output list to download models to a folder
        * `python fetch_af_pdb.py -i alis_afdb50.m8 -o Output/folder` 
4. Run `align_all.py` to perform pairwise alignment and TM-score calculation for all downloaded files
   * `python align_all.py -i Output/folder -o afdb50_tm_scores.csv -j 64`
   * This will generate an output csv in upper triangular format containing the TM-scores between all pairs of structures. For large datasets, this will take a while. Set threads with `-j` flag to speed this up.
   * To also generate a set of PDB files aligned to a reference structure, supply the reference PDB with the `-r` flag. You can also specify a directory to save the aligned PDBs to with the `-d` flag: by default it will save them to a folder named `aligned` in the current directory.
     * `python align_all.py -i output/folder -o afdb50_tm_scores.csv -j 64 -r AF-O25823-F1-model_v4.pdb`
   * This expects the USAlign executable to be in the PATH, or at `~/bin/USalign`. If it is not, you must specify the location of the executable with the `-x` flag.
5. Run `merge_pdbs.py` to combine all the aligned PDB files into a single file
   * `python merge_pdbs.py -i output/aligned -o aligned_models.pdb -r output/folder/AF-O25823-F1-model_v4.pdb`
   * This will generate a single PDB file containing all the aligned structures. This is required to generate the pseudo-density map.
6. Open the output merged PDB file in PyMOL and save it as a *single* PDB file containing all models
7. Run `pseudo_density.py` to generate a pseudo-density map from the combined PDB file
   * `python pseudo_density.py -i aligned_models.pdb -o pseudo_density.dx`
   * This will generate a pseudo-density map in DX format that can be opened in, e.g., ChimeraX.
   * You can set the grid spacing (in Angstroms) for sampling with the `-d` flag. The default is 0.5.
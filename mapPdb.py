from Bio import SeqIO
from urllib.request import urlretrieve
import pandas as pd
import numpy as np
import json
import sys
import os
import yaml
import subprocess

def align_pdb(pdbid,
              keep,
              pdb_fasta,
              gen_fasta,
              rdict,
              datapath,):
    '''
    Maps chain indices from a PDB file to genomic indices

    Params
    ______

    pdbid: str
        PDB ID for structure
        
    keep: list
        names of chromosomes in genomic FASTA file to keep
        
    pdb_fasta: str
        Path to a FASTA file containing chains from the PBD structure

    gen_fasta: str
        Path to the genomic FASTA used for alignment

    rdict: dict
        Dictionary renaming chromosomes in pdb_fasta to "pbd_{nameInGenomicFasta}"

    datapath: str
        destination for output files 

    Returns
    _______

    output : one json file per chain
        A json file with genomic indices as keys and PDB indices as values
        
    '''
    # Gather chromosomes from genome file
    hg_dict = {}

    with open(gen_fasta,'r') as f:
        for record in SeqIO.parse(f,'fasta'):
            if record.id in keep:
                hg_dict['gen_'+record.id] = str(record.seq
                                                ).replace('T','U')
            else:
                pass

    # Import PDB fasta and rename chromosomes
    with open(pdb_fasta,'r') as f:
        ls = f.readlines()

    pdb_dict = {}
    
    for i,v in enumerate(ls):
        newname = ''
        if v[0] == '>':
            newname = v.split('|')[2].replace(' ','-')
            if newname in rdict:
                pdb_dict[rdict[newname]] = ls[i+1]
            else:
                pass
        else:
            pass

    # Produce fastas for muscle alignment
    f_ls = []
    
    for k in keep:
        fname = datapath + '_'.join([pdbid,k,'toAlign.fa'])
        f_ls.append(fname)
        with open(fname,'w') as f:
            h_name = 'gen_'+k
            p_name = 'pdb_'+k
            f.write('>'+h_name+'\n')
            f.write(hg_dict[h_name]+'\n')
            f.write('>'+p_name+'\n')
            f.write(pdb_dict[p_name])

    # Write bash scripts for muscle alignment
    bash = []
    afa = []
    for f in f_ls:
        sh = f.split('.')[0]+'.sh'
        bash.append(sh)
        out = f.split('.')[0]+'.afa'
        afa.append(out)
        cmd = 'muscle -align '+f+' -output '+out
        with open(sh,'w') as s:
            s.write('#!/bin/bash'+'\n')
            s.write(cmd)

    # Run muscle alignment
    for b in bash:
        cmd = ['sh',b]
        subprocess.run(cmd)

    # Parse alignment files
    afa_dict = {}
    
    for a in afa:
        with open(a,'r') as f:
            for record in SeqIO.parse(f,'fasta'):
                afa_dict[record.id] = str(record.seq)

    # Export json files
    for k in keep:
        # Get sequences
        gen_seq = list(afa_dict['gen_'+k])
        pdb_seq = list(afa_dict['pdb_'+k])

        # Determine if base or space
        gen_inds = []
        for i,v in enumerate(gen_seq):
            if v == '-':
                gen_inds.append(v)
            else:
                gen_inds.append(i+1)

        pdb_inds = []
        for i,v in enumerate(pdb_seq):
            if v == '-':
                pdb_inds.append(v)
            else:
                pdb_inds.append(i+1) 

        # Match indices
        align_ls = list(zip(gen_inds,pdb_inds))
        align_ls = [i for i in align_ls if '-' not in i]
        align_dict = dict(align_ls)

        # Make json file
        fname = datapath+pdbid+'_'+k+'_mapped.json'
        with open(fname,'w') as jfile:
            json.dump(align_dict,jfile,indent=1)
    
    return print(f'Created alignment files in {datapath}.')

def hex_to_rgb(hex_color):
    """
    Written by Claude.ai
    Convert hex color code to RGB tuple.
    
    Args:
        hex_color: String hex code (e.g., '#FF0000' or 'FF0000' or '#ff0000')
    
    Returns:
        Tuple of (r,g,b) values as integers
    """
    # Remove '#' if present
    hex_color = hex_color.lstrip('#')
    
    # Convert hex to RGB
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    
    return (r, g, b)

def interpolate_color(color1, color2, value):
    """
    Written by Claude.ai
    Interpolate between two RGB colors based on a value between 0 and 1
    
    Args:
        color1: Hex code for minimum value
        color2: Hex code maximum value
        value: Float between 0 and 1
    
    Returns:
        String with hex color code
    """
    r1, g1, b1 = hex_to_rgb(color1)
    r2, g2, b2 = hex_to_rgb(color2)
    
    # Linear interpolation between the two colors
    r = int(r1 + (r2 - r1) * value)
    g = int(g1 + (g2 - g1) * value)
    b = int(b1 + (b2 - b1) * value)
    
    # Convert to hex string
    return f"#{r:02x}{g:02x}{b:02x}"

def export_chains(pdbid,
                  rna_info,
                  protein_info,
                  datapath):
    '''
    Export chain info to csvs for coloring with ChimeraX

    Params
    ______

    pdbid : str
        PDB ID, used only for naming files

    rna_info : dict
        Dictionary in which keys are hg38 ncRNA chromosome
        names and values are a two-object tuple where the 
        first value is the PDB chain ID and the second is
        the hex code used to color the chain

    protein_info : tuple
        A tuple of tuples. Each tuple should have the following structure:
        (SPIDR Target Name, PDB Chain ID, RGB Hex Code, RNA Name from hg38 
        ncRNA)

    datapath : str
        Path to directory where csvs will be deposited

    Returns
    _______

    output : a csv
        Outputs a csv containing color information that ChimeraX will use
        to style chains
    '''
    # Export chain info
    chain_ls = [i[1:3] for i in protein_info if i[1] != '']
    chain_ls = chain_ls + list(rna_info.values())
    
    fname = datapath+pdbid+'_chainInfo.csv'
    pd.DataFrame(chain_ls).to_csv(fname,header=False,index=False)

    return fname

def write_chimerax(pdbid,
                   rna_info,
                   chain_file,
                   cmap_file,
                   datapath,
                   template='chimeraxTemplate.py'):
    '''
    Writes a ChimeraX-readable script based on a template.
    Uses outputs from bedtocmap() and export_chains() as inputs.

    Params
    ______

    pdbid : str
        PDB ID. Tells ChimeraX which file to open

    rna_info : dict
        Dictionary in which keys are hg38 ncRNA chromosome
        names and values are a two-object tuple where the 
        first value is the PDB chain ID and the second is
        the hex code used to color the chain

    chain_file : str
        Relative path to the output of export_chains()

    cmap_file : str
        Relative path to the output of bedtocmap()

    datapath : str
        Path to directory where .py file will be deposited

    template : str
        Path to a text file for the constant parts of the 
        ChimeraX script. Appended to the end of the variables.
        Default: chimeraxTemplate.py

    Returns
    _______
    
    output : a .py file
        A ChimeraX executable Python script for styling the 
        indicated PDB strucutre
    
    '''
    # Get text of template file
    with open(template,'r') as f:
        templs = f.readlines()

    # Grab the RNA chains
    rna_chains = [i[0] for i in rna_info.values()]

    # Get absolute paths for chain and cmap files
    chain_path = os.path.abspath(chain_file)
    cmap_path = os.path.abspath(cmap_file)

    # Get path to save cxs file
    cxs_dir = os.path.dirname(cmap_path)+'/'
    cxs_file = cxs_dir+pdbid+'_sesh.cxs'

    # Make list of lines for top of file
    top_ls = ["pdbid = '"+pdbid+"'\n",
              "rna_chains = "+str(rna_chains)+"\n",
              "chain_file = '"+chain_path+"'\n",
              "cmap_file = '"+cmap_path+"'\n",
              "cxs_file = '"+cxs_file+"'\n"
              "\n"]

    # Merge the lists
    ls_all = top_ls + templs
    # Write out file
    fname = datapath+pdbid+'_chimeraX.py'
    
    with open(fname,'w') as f:
        f.writelines(ls_all)
    return print(f'ChimeraX script written to {fname}')

def bedtocmap(pdbid,
              rna_info,
              protein_info,
              spidrpath,
              datapath,
              floor=2,
              mask=0.2,
              binsize=None):
    '''
    Using jsons exported from align_pdb(), map nucleotides in SPIDR
    data to those in PDB structures and assign colors based on SPIDR
    enrichment values

    Params
    ______
    
    pdbid : str
        PDB ID, used only for naming files

    rna_info : dict
        Dictionary in which keys are hg38 ncRNA chromosome
        names and values are a two-object tuple where the 
        first value is the PDB chain ID and the second is
        the hex code used to color the chain

    protein_info : tuple
        A tuple of tuples. Each tuple should have the following structure:
        (SPIDR Target Name, PDB Chain ID, RGB Hex Code, RNA Name from hg38 
        ncRNA)

    spidrpath : str
        Path to a tidy bedgraph of SPIDR data. Columns should be chr, start, end,
        enrichment, target.

    datapath : str
        Path to directory where csvs will be deposited

    floor : int or float 
        Lowest enrichment value to be considered for normalization and 
        colorization. Default 2

    mask : False or float
        Lowest fraction of maximum enrichment value to be assigned 
        a color. If False, nothing above the floor is masked. Default 
        0.2

    binsize : None or int
        Size of genomic bins for data. If None, 1 nt bins are used.
        If int, bins are generated centered at each nucleotide. 
        Surrounding nucleotides in the bin are set to the maximum
        value of that bin. Default None

    Returns
    _______

    output : two csv files
        The first csv file is generated by export_chains() and contains 
        color information for chains to be included in ChimeraX. The second
        csv contains nucleotide by nucleotide color codes for enriched SPIDR
        bins
    '''
    # Export color info for included chains
    chain_file = export_chains(pdbid,
                               rna_info,
                               protein_info,
                               datapath)

    # Import alignments
    align_dict = {}
    
    for r in rna_info.keys():
        fname = datapath+'_'.join([pdbid,r,'mapped.json'])
        with open(fname,'r') as jfile:
            align_dict[r] = json.load(jfile)
    
    # Load tidy SPIDR data
    spidr = pd.read_csv(spidrpath,header=None)
    spidr[4] = [i.upper() for i in spidr[4]]

    # Get PDB indices and add colors
    df_ls = []
    for p in protein_info:
        # Filter by target
        subdf = spidr[spidr[4] == p[0]].copy()
    
        # Remove low enrichments
        subdf = subdf[subdf[3] >= floor]
    
        # Normalize
        subdf['norm'] = subdf[3] / np.max(subdf[3])

        # Mask low regions
        if mask:
            subdf = subdf[subdf['norm'] > mask]
        else:
            pass
            
        # Filter to rna and assign colors
        subdf = subdf[subdf[0] == p[3]]
        pdb_dict = align_dict[p[3]]
        color1 = '#FFFFFF' #rna_info[p[3]][1]
        color2 = p[2]
        for i in subdf.index:
            subdf.loc[i,'pdb_ind'] = str(pdb_dict[str(subdf.loc[i,1])])
            subdf.loc[i,'color'] = interpolate_color(color1, 
                                                     color2, 
                                                     subdf.loc[i,'norm'])
    
        df_ls.append(subdf)
    
    df_all = pd.concat(df_ls).reset_index(drop=True)
    
    # Filter duplicates, keeping only maximum value
    for r in rna_info:
        subdf = df_all[df_all[0] == r].copy()
        subdf = subdf.loc[subdf.duplicated('pdb_ind',keep=False)]
        if len(subdf) == 0:
            pass
        else:
            for i in np.unique(subdf['pdb_ind']):
                dupes = subdf[subdf['pdb_ind'] == i]['norm']
                maxind = dupes.idxmax()
                df_all.drop(dupes.drop(maxind).index,inplace=True)

    # Binning if desired
    if binsize:
        # Generate range for binning
        if binsize % 2 == 0:
            binhalf = int(binsize/2)
            binrange = (binhalf-1,binhalf+1)
        else:
            binhalf = int((binsize-1)/2)
            binrange = (binhalf, binhalf)

        # Make ChimeraX readable ranges
        starts = [str(int(i)-binrange[0]) for i in df_all['pdb_ind']]
        ends = [str(int(i)+binrange[1]) for i in df_all['pdb_ind']]
        indrange = ['-'.join(i) for i in zip(starts,ends)]

        # Sort by norm so that highest values are set last
        df_all['pdb_ind'] = indrange
        df_all = df_all.sort_values('norm')
    else:
        pass

    # Add chain to nucleotides
    for i in df_all.index:
        df_all.loc[i,'rna_chain'] = rna_info[df_all.loc[i,0]][0]
    
    df_all['chainind'] = df_all['rna_chain']+':'+df_all['pdb_ind']

    # Export chain IDs and hex codes to a csv
    outls = list(zip(df_all['chainind'],df_all['color']))
    
    cmap_file = datapath+pdbid+'_cmap.csv'
    pd.DataFrame(outls).to_csv(cmap_file,header=False,index=False)

    # Write ChimeraX script
    write_chimerax(pdbid,
                   rna_info,
                   chain_file,
                   cmap_file,
                   datapath,
                   template='chimeraxTemplate.py')

    return print(f'Nucleotide colormapping deposited in {fname}')

if __name__ == "__main__":
    # Written with assistance of Claude.ai
    # Check if a YAML file is provided as an argument
    if len(sys.argv) < 2:
        print("Please provide a YAML file path")
        sys.exit(1)
    
    # Get the YAML file path from command-line arguments
    yaml_file_path = sys.argv[1]

    with open(yaml_file_path, 'r') as file:
        params = yaml.safe_load(file)

        # Make data path if it doesn't exist
        pdbid = params['shared_params']['pdbid']
        if not os.path.isdir(pdbid):
            os.makedirs(pdbid)

        # Download FASTA from PDB
        datapath = pdbid+'/'
        pdb_fasta  = ''.join([datapath,'rcsb_pdb_',pdbid,'.fasta'])
        url = 'https://www.rcsb.org/fasta/entry/'+pdbid
        urlretrieve(url, pdb_fasta)

        params['align_pdb']['datapath'] = datapath
        params['align_pdb']['pdb_fasta'] = pdb_fasta

        params['bedtocmap']['datapath'] = datapath        
        
        align_pdb(**params['align_pdb'])
        bedtocmap(**params['bedtocmap'])


                

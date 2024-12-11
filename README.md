# mapPdb: A script for mapping CLIP, CLAP, or SPIDR data to PDB structures

### Overview

This repo contains a script for coloring structural data from PDB based on RNA sequencing enrichment values. The basic workflow is as follows:

    1. The sequences of the desired RNAs are pulled from the specific PDB structure as well as from the genomic FASTA file used for alignment.
    2. The PDB and genomic sequences are aligned using Muscle. The alignment files are used to match nucleotide positions in the genomic file to those in the structure.
    3. The new nucleotide positions are used to update a modified 1nt bedgraph containing the sequencing enrichment values and the identities of protein targets. Additionally, the bedgraph is binned and normalized according to user specifications. 
    4. A custom colormap is generated for each protein target and hex codes are generated for each enrichment value. A new file is generated containing the chain, residue, and color information that will style the PDB structure.
    5. A ChimeraX readable script is generated based on a template and the colormapping information.
    6. The user opens ChimeraX and executes the script. This will open the desired PDB structure, apply the styling specified in the template and colormapping files, and save a ChimeraX session file.

All of these steps are wrapped in a single script: mapPdb.py. mapPdb.py takes a single argument: the seqInfo.yaml file. The yaml file supplies all necessary information for finding the correct PDB structure and applying colors as desired, and also points to the sequencing data.

Sequencing data for mapping should be stored in a modified 1nt bedgraph file referred to here as a "tidy bedgraph". The tidy bedgraph has five columns:

    - Column 1: RNA name 
    - Column 2: Start nucleotide of bin
    - Column 3: End nucleotide of bin
    - Column 4: Enrichment 
    - Column 5: Experimental target

The tidy bedgraph format allows for multiple targets to be mapped to the same structure. For example, if you have CLIP data for RPS5 and RPL36, two tidy bedgraphs can be concatenated without losing the information of which enrichment values belong to each protein.

mapPdb is a relatively simple tool for comparing sequencing and structural data. The ChimeraX script is usually generated in a few seconds, making it easy to adjust styling preferences in the seqInfo.yaml or chimeraxTemplate.py files and rerun the pipeline. From there, publication quality figures can be generated using standard ChimeraX tools.

### Getting started

The final ChimeraX script generated is compatible with ChimeraX 1.8 but has not been tested on other versions. ChimeraX can be downloaded [here](https://www.cgl.ucsf.edu/chimerax/download.html). 

After cloning this repo, create the mapPdb Conda environment. I recommend using [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for package management. On most systems, run:

```mamba env create -f envs/mapPdb.yaml```

Then, before running the script, run:

```conda activate mapPdb```

### Running the example

To make sure that the installation is working, run the example pipeline using:

```python mapPdb.py example/4ug0.yaml```

This should output a folder called `4ug0/`. In `4ug0/` is a script called `4ug0_chimeraX.py`. Open ChimeraX and run the script from the ChimeraX command line:

`run /absolute/path/to/mapPdb/4ug0/4ug0_chimeraX.py`

ChimeraX requires absolute paths for executing scripts. If everything has run correctly, the 4UG0 structure should appear with protein and RNA coloration as specified in the `4ug0.yaml` config file. In the `4ug0/` directory the ChimeraX session is now stored as `4ug0_sesh.cxs`.

Note that all chains other than those specified in the config file have been deleted. In contrast, all RNA residues that did not receive a color have been rendered invisible. To reveal the surfaces of those bases, run:

`show /S2 /L5 surfaces`

in the ChimeraX command line. The surfaces can be hidden using:

`hide /S2 /L5 surfaces`

For a guide on other commands for styling and saving images, consult the [ChimeraX documentation](https://www.cgl.ucsf.edu/chimerax/docs/user/index.html).

### Customizing the config file

The file `seqInfo.yaml` contains instructions on how to adjust the config file for a new visualization. 
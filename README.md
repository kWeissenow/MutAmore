# MutAmore

[![Protein mutation movie](http://img.youtube.com/vi/1XgiFXg-Xrs/0.jpg)](https://www.youtube.com/watch?v=1XgiFXg-Xrs&list=PL0QUUE_zWBuJ6Y5NWtDoY93FUweUUGVuf)

MutAmore visualizes the predicted structural changes of single amino-acid substitutions on a protein sequence as a protein mutation movie (PMM).

Examples rendered in 4K resolution are available on [YouTube](https://www.youtube.com/watch?v=1XgiFXg-Xrs&list=PL0QUUE_zWBuJ6Y5NWtDoY93FUweUUGVuf).

# Running MutAmore on Google Colab
You can try MutAmore on Google Colab in [this notebook](https://colab.research.google.com/drive/1tOzIw2JLnj_jWOHcGQ-kpbSw790O_5OY?usp=sharing).
The Colab notebook is sufficient for movie rendering of (shorter) protein sequences. If you want to create movies for a larger set of sequences or particularly large proteins, consider installing MutAmore locally as described below.

# Local installation
## Requirements

MutAmore is written in Python and requires the Python packages Pillow, Scipy and Biopython:
```
pip install Pillow scipy biopython
```
For movie rendering, you will also need to have [PyMOL](https://pymol.org/2/) and ffmpeg installed. These packages are best installed using the package manager of your Linux distribution, e.g. for Debian/Ubuntu systems:
```
apt-get install pymol
apt-get install ffmpeg
```
Alternatively, you can use Conda:
```
conda install -c conda-forge pymol-open-source
conda install -c conda-forge ffmpeg
```

You will need a structure prediction system for the prediction of mutated structures. We recommend using [ESMFold](https://github.com/facebookresearch/esm) or [ColabFold](https://github.com/sokrypton/ColabFold), but any alternative structure predictor can be used (see below).

## Running MutAmore

MutAmore expects one or more amino-acid sequences in FASTA format as input.

```
python run_movie_rendering.py -i <INPUT_FASTA> -s predictor_scripts/esmfold.sh
```

The predictor script tells MutAmore how to call the structure prediction system. We provide two scripts for ESMFold and ColabFold.

### Using ESMFold

Before running MutAmore, edit the script `predictor_scripts/esmfold.sh` and replace `ESMFOLD_LOCATION` with the path to your ESMFold installation.

### Using ColabFold

Running ColabFold in single-sequence mode works out-of-the-box with the provided script `predictor_scripts/colabfold.sh`. If you want to run ColabFold with locally generated multiple-sequence alignments, refer to the [documentation](https://github.com/sokrypton/ColabFold#generating-msas-for-large-scale-structurecomplex-predictions) and include calls to `colabfold_search` in the predictor script.

### Using an alternative structure prediction system

Any alternative structure prediction system can be used, as long as it expects FASTA files as inputs and generates PDB files.
Create a new script in `predictor_scripts/` and fill in the command-line instruction to run the structure prediction system. Use the placeholders MUTAMORE_INPUT and MUTAMORE_OUTPUT as the input FASTA file and output directory respectively. MutAmore will replace these placeholders at runtime with the appropriate parameters.

## Running MutAmore on two seperate systems

Structure prediction systems usually require powerful GPUs and are often run in server environments. Some users might not be able to install graphical packages such as PyMOL and ffmpeg in their server environment. In this case, you can run MutAmore in two seperate steps on different machines.
You can run only the structure prediction part of the pipeline by specifying an output directory for the structure predictions:
```
python run_movie_rendering.py -i <INPUT_FASTA> -s predictor_scripts/esmfold.sh --only-predict -p <PREDICTION_DIRECTORY>
```

You can then transfer the prediction directory to another machine, e.g. a desktop computer, to run the rendering:
```
python run_movie_rendering.py -i <INPUT_FASTA> --only-render -p <PREDICTION_DIRECTORY>
```

# Additional parameters
### Movie resolution
MutAmore by default renders movies in the resolution 1280x720. You can set a different resolution with the parameters `width` and `height`, e.g. to render in 4K resolution add the following to the MutAmore call:
```
--width 3840 --height 2160
```

### Temporary directory
MutAmore creates a temporary directory to store intermediate files, which is deleted after finishing. By default, it will create the directory `./tmp/` where you call MutAmore. To specify a different directory, use the parameter `-t`:
```
-t <TEMPORARY_DIRECTORY>
```

### Prediction output
MutAmore stores structure predictions in a temporary directory (see above), which is deleted after rendering the movie. In case you want to keep the structure prediction output, you can specify a different directory with the parameter `-p`:
```
-p <PREDICTION_DIRECTORY>
```

### Output directory
By default, MutAmore outputs the finished movies in the directory in which you call MutAmore. To output to a different location, use the parameter `-o`:
```
-o <MOVIE_OUTPUT_DIRECTORY>
```

# Cite
The pre-print is available on [bioRxiv](https://biorxiv.org/cgi/content/short/2023.09.15.557870v1).
Cite with
```
@article {Weissenow2023,
	author = {Weissenow, Konstantin and Rost, Burkhard},
	title = {Rendering protein mutation movies with MutAmore},
	elocation-id = {2023.09.15.557870},
	year = {2023},
	doi = {10.1101/2023.09.15.557870},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/10.1101/2023.09.15.557870v1},
	eprint = {https://www.biorxiv.org/content/10.1101/2023.09.15.557870v1.full.pdf},
	journal = {bioRxiv}
}
```
or
```
K. Weissenow and B. Rost, Rendering protein mutation movies with MutAmore, bioRxiv (2023)
```

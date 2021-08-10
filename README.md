# ZpAnomalon Analysis - RJR to Plots

This is the repo of the third, fourth, and plotting  steps of the ZpAnomalon Analysis. In this analysis, We are using Recrusive Jigsaw Reconstruction to generate the mass estimators of each daughter particle in the Z' decay chain. This calculcation is done with the Restframes package, whose installation details are outlined below. A c++ TMakeClass analyzer is used to calculate the RJR quantiies, and flatten the output trees from [trimmerShed](https://github.com/gracecummings/trimmerShed). The flattened trees are passed to a python native dataframe based analyzer, which creates the analysis object hisotgrams.

## Setting up the analysis environment

This environment was designed to run locally, on an Ubuntu 18.04 machine. There are a lot of packages that this requries, and the below works to set this up locally. At this time, the procedures have not been vetted on the LPC, but the installation of RestFrames, and maybe the Uproot 3 versus 4 distinction, would probably be the only thing that would not work out-of-the-box.

The following python packages are required:

+ pyroot
+ uproot3
+ pandas
+ numpy
+ boost-histogram
+ configparser

if you already have Python 3 working with ROOT, with the above dependences, skip to "Download Repository and Setup Restframes." If you do not, you can follow the conda environment instructions below.

### Conda-forge: Setting up all your needs

We will use a conda environment to keep all versions in sync, and save us the headache of trying to keep ROOT up-to-date with the rest of the universe. These directions are based on Bob Hirosky's directions for the University of Virginia's 56XX computational physics courses. The reasoning behind using conda, what it is, and all that jazz is found [there](https://raw.githubusercontent.com/UVaCompPhys/PHYS56xx-setup/master/notebooks/FAQ/WorkingEnv.ipynb).

To get started, we will install minconda.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

```

This edits your ```.bashrc``` file, so either restart your shell or ```source ~/.bashrc```. The defailt is to start a ```(base)``` conda envoriment upon startup, as indicated from your now altered PS1. I despise this, so to turn off this startup feature, execute

```conda config --set auto_activate_base false```

Now we can begin setting up the useful environment. We will be using the conda-forge repository as the base.

```bash
conda config --add channels conda-forge
conda config --show channels
conda config --set channel_priority strict
```

Setup the conda environment to run this analysis.

```bash
conda create -n ZpAnomalon
conda env list
conda activate ZpAnomalon
```

To deactivate the environment, simply use ```conda deactivate ZpAnomalon```.

Install the packages.

```bash
conda install --yes scipy matplotlib seaborn sympy scikit-learn
conda install --yes root
conda install -c conda-forge pandas
conda install -c conda-forge uproot
conda install -c conda-forge uproot3
conda install -c conda-forge configparser
conda install -c conda-forge boost-histogram
```
To see what you have installed, ```conda list```. The python environment is ready. Proceed to the next section to implement Restframes.

### Download Repository and Setup RestFrames

Download this analysis repository. I would recommend forking it, with the icon in the top left. Detailed directions on how to fork, and sync with this repo, can be found in the [Github Docs](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo). **The current branch of interest in recluster_cluster.** There is no special setup script or compiling. Once you have the repository setup locally,

```bash
cd ZpAnomalonAnalysisUproot
```

We use the [RestFrames](http://restframes.com/) package to calculate the mass esitmators of the Z' cascade decay. Make sure you have you conda environment activated, to compile Restframes with the appropriate ROOT version. I install RestFrames in the repo directory, to make it easy to keep track, and to avoid possible conflicts on my larger system.

```bash
git clone https://github.com/crogan/RestFrames RestFrames
cd RestFrames
./configure --prefix=$PWD
make
make install
```

## Making Topiary (Analysis Third Step)

As stated in the intro, the third step in our analysis chain is c++ TMakeClass analyzer used to trim the skims produced by the trimmer in trimmerShed, do the RJR calculation, and flatten the trees for easier Pandas dataframe use. This is done with a ROOT macro, exectued via a python script. To call RestFrames libraries from a ROOT macro, you must first source the setup script. From the ```ZpAnomalonAnalysisUproot``` directory,

```bash
source RestFrames/setup_RestFrames.sh
```

Currently, these scripts are designed to run locally, with the [trimmerShed](https://github.com/gracecummings/trimmerShed) produced skims stored somewhere on a local machine. Ideally, this step will be combined with the trimmerShed step, but that is for the future. The path can be set by changing

```python
inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
```

in ```runTopiary.py``` to the correct path on your machine. If the source ROOT file is in the executable directory, and no TChaining is needed, you can also just pass the entire name of the ROOT file directly. An example executable command to produce a topiary is:

```bash
python runTopiary.py -s Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8 -c mumu
```

Executables in this framework are naming sensitive, so be sure to include the [ZpAnomalon_TreeMaker](https://github.com/gracecummings/ZpAnomalon_TreeMaker) naming conventions in the sample strings. This naming convention automatically handles the year (for triggers) and the  sample type (for k-factors and triggers). The channel must be explicitly indicated, with either "mumu", "ee" or "emu" selected. Currently only the mumu channel is fully supported.

## Doing selections (Analysis Fourth Step)

The fourth step in our analysis chain is a Python-based analyzer that does our cuts. The flattened topiaries are built for analysis using a dataframe, which is where Pandas and company come into play. The fourth step produces histograms that pass all Pre-selection cuts, impose the binning, and calculate the stats error bars using Sumw2. At this step, the k-factors are applied to DY+Jets samples. No xs, lumi, or stitching is applied at this step. An example excutable is:

```bash
python doSelections.py -f Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8 -zpt 150.0 -hpt 300.0 -met 200.0 -sdm 30.0 -b DeepMassDecorrelTagZHbbvsQCD -wp 0.8 -date 2021-08-10
```

Let's go through some options:

+ `-f` is the sample string
+ `-zpt` is the Z candidate pT cut
+ `-hpt` is the Higgs candidate pT cut
+ `-met` is the MET cut
+ `-sdm` is the lower cut on soft drop mass
+ `-b` is th btagger name to use
+ `-wp` is the btagger working point to use
+ `-date` is actually the directory in the "analysis_output_ZpAnomalon" folder where the input to this selections step is held. Step 3 (Topiary) automatically creates the aforementioned folder, and outputs are nested in directories named for the date the file was created, so this takes a date, classically.
+ `-sr True` is a flag for is signal region selections should be used. Default are sidebands
+ `-c True` is a flag for if the entire region (no signal/sideband differentiation) should be used.
+ `-v True` is a flag to use validation region sidebands and fake signal regions. Must be combined with region flags.

Output for this step is stored in "analysis_output_ZpAnomalon/DATEOFEXECUTION". Three files are created, with different names and file types. The three files produced are identified as:

+ sampleName+`upout`+region+cutdetails+`.root` - ROOT file with weighted preselection histograms
+ sampleName+`selected_errors`+region+cutdetails+`.pkl` - Pickle file with bin-by-bin stats errors
+ sampleName+`totalevents`+region+cutdetails+`.npy` - numpy file with the original event counts

These three files are used in the plotting steps of the analysis.

## Other manipulations and facillitated respins

In many cases, one wants to do go from new selections to plots all in one go. There is a script for that. 

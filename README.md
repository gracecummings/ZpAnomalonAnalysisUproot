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

Download this analysis repository. I would recommend forking it, with the icon in the top left. Detailed directions on how to fork, and sync with this repo, can be found in the [Github Docs](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo). Ther eis no special setup script or compiling. Once you have the repository setup locally,

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

in ```runTopiary.py``` to the correct path on your machine. An example executable command to produce a topiary is:

```bash
python runTopiary.py -s Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8 -y 2017
```

Executables inthis framework are naming sensitive, so be sure to include the [ZpAnomalon_TreeMaker](https://github.com/gracecummings/ZpAnomalon_TreeMaker) naming conventions in the sample strings.

## Doing selections (Analysis Fourth Step)

## Other manipulations and facillitated respins

In many cases, one wants to do go from new selections to plots all in one go. There is a script for that. 

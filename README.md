# ZpAnomalon Analysis - RJR to Plots

This is the repo of the third, fourth, and plotting  steps of the ZpAnomalon Analysis. In this analysis, We are using Recrusive Jigsaw Reconstruction to generate the mass estimators of each daughter particle in the Z' decay chain. This calculcation is done with the Restframes package, whose installation details are outlined below. A c++ TMakeClass analyzer is used to calculate the RJR quantiies, and flatten the output trees from trimmerShed. The flattened trees are passed to a python native dataframe based analyzer, which creates the analysis object hisotgrams.

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

### Conda forge: Setting up all your needs

### Download Repository and Setup RestFrames.

## Making Topiary (Analysis Thrid Step)

As stated in the intro, the third step in our analysis chain is c++ TMakeClass analyzer used to trim the skims produced by the trimmer in trimmerShed, do the RJR calculation, and flatten the trees for easier Pandas dataframe use.

## Doing selections (Analysis Fourth Step)

## Other manipulations and facillitated respins

In many cases, one wants to do go from new selections to plots all in one go. There is a script for that. 

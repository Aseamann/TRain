# TRain
TRain: T-cell Receptor automated immunoinformatics v1.0
Copyright (c) 2024 Austin Seamann.

Thank you for your interest in TRain. For more information, visit the
publication:

URL HERE

If you use this program, please cite the paper.

--------------------------------------------------------------------------------

## Visual Guide

![Flowchart of the programs of TRain](/images/Figure3_TRain_Full_6x.png)

--------------------------------------------------------------------------------

## User Guide

The User's Guide can be found in this GitHub repository as TRain\_User\_Guide.pdf

--------------------------------------------------------------------------------

## Installation
The code for  **TRain**  is available through GitHub at: https://github.com/Aseamann/TRain. Download  
or clone the latest version to your location of choice.

Using git command line:

    git clone https://github.com/Aseamann/TRain
 
 Recommended Setup:
 

    # Enter TRain directory
    cd TRain
    # Setup conda environment with python3.9
    conda create -n TRain python=3.9
    # Modify "TRain/data/config.ini" to point to Rosetta install (see below)
    # Use preferred text editor (Vi, Nano, etc)
    open data/config.ini
    # Build python package
    python -m pip install -e .
    # Activate the conda environment - needed everytime TRain is in use
    conda activate TRain
For modeling and docking, a local version of Rosetta needs to be installed. For licensing details  
visit https://www.rosettacommons.org/software/license-and-download. Installation instructions can  
be found at https://new.rosettacommons.org/docs/latest/getting_started/Getting-Started. Rosetta  
should be built with MPI executables and instructions can be found at  
https://new.rosettacommons.org/docs/latest/build_documentation/Build-Documentation (or follow  
the instructions below). While there is the option to run without MPI, it is not recommended. For more  
instructions on how to run  ModelEngine  and for a list of command line options, see Chapter 5 in the User's Guide.  

If MPI is not installed on your computer, you can install it using the following command (for a  
system with an apt installer):

    # Linux
    sudo apt install libopenmpi-dev
    # Or visit: https://www.open-mpi.org/faq/?category=building
    
    # MacOS
    xcode-select --install
    brew install openmpi
    # Or visit: https://www-lb.open-mpi.org/software/ompi
 
Please read below for general build options of Rosetta for Linux and MacOS systems. RosettaCommons forums are a good place to resolve compilation issues.

    # Navigate to rosetta download
    cd rosetta_src_xx/main/source
    # Recommended build command
    ./scons.py bin mode=release extras=mpi -j<number_of_processors_to_use>

    # MPI not working? Potentially need to point to mpicxx and mpicc
    # Print location of mpicxx and mpicc
    which mpicxx
    which mpicc

    # Use preferred text editor to open site.settings
    open tools/build/site.settings
    
    # Uncomment cxx and cc override xml options,
    # replace with location of your devices executables

The configuration file for  TCRcoupler  is inside directory  /TRain/data/config.ini. Modify the  
configuration file with your preferred text editor specifying the location of Rosetta installation folder  
in line 2 inside the quotation marks.  

Each step of this Quick Start guide contains output files, so you can compare what you obtain  
and make sure each step runs correctly.

## Input
For all the steps, we assume that your working directory is "/TRain". The first step is to construct a full TCR protein sequence using the gene segments and CDR3 regions listed.

    cd tr01input
    ls
    SeqConductor Sample_Input/Test_Table_Single.xlsx -f -a
 The files *alpha.fasta* and *beta.fasta* should appear in your working directory. The contents should be identical to the fasta files location in "tr01input/Sample_Output/Quick_Start/" and "tr01input/Sample_Input/Quick_Start/".

## Model
Now we can submit our fasta files for modeling. For this step we will need a working installation of Rosetta.

    cd ../tr02model
    ls
    ModelEngine -a Sample_Input/Quick_Start/alpha_single.fasta -b Sample_Input/Quick_Start/beta_single.fasta
The resulting model will be under the directory “Modeled”. You can open this PDB file with your  
PDB file viewer of choice. Only the variable region of the TCR will be modeled. It should match the  
model found in "tr02model/Sample_Output/Quick_Start/" and "tr03pair/Sample_Input/Quick_Start/".

## Pair
The TCR we chose binds to the  Influenza  M1 antigen. Therefore, we will pair it with an M1 antigen  
and MHC from PDB model  *5ISZ*. This step will prepare the two structures for the subsequent  
docking step.

    cd ../tr03pair
    TurnTable -t Sample_Input/Quick_Start/ES179M1-01.pdb -p Sample_Input/5isz.pdb
The paired TCR and pMHC will now be in the PDB file under the new directory "Paired". The new name of the PDB will start with the TCR name and then the pMHC name.

## Dock
For this Quick Start guide, we will perform a very short docking run to ensure that everything is working properly. There are two sample runs listed so you can determine if the MPI binaries are working on your system, or else you can use Rosetta without MPI. These examples are not recommended for actual biological interpretation of the interaction, as they are very short. To ensure you're utilizing **TRain** properly, visit chapter 6. 

We will begin by creating a new directory to carry out docking. This is done so our output can be organized, as several sub-directories will be created for each docking step.

    cd ../tr04dock
    mkdir QS_Dock
    cd QS_Dock
    cp ../Sample_Input/Quick_Start/ES179M1-01_5isz.pdb .

    # Example 1 (w/o MPI)
    TCRcoupler ES179M1-01_5isz.pdb -q -r -a 10 -d 100 -e 10
    
    # Example 2 (with MPI)
    TCRcoupler ES179M1-01_5isz.pdb -r -a 10 -d 100 -e 10
Since we are testing with the rigid docking protocol, you will find three folders within the  
“/output_files” directory. Within the folder there will be the “relax”, “dock”, and “refine” directories.  
The best scoring refined PDB will be found in the “refine” directory.  
We will now take the final structure and renumber the pMHC to what was present in 5isz.pdb.  
Your input to PostCoupler.py will be unique, so please follow the instructions below.

    cd ../
    PostCoupler QS_Dock/output_files/refine/ES<tab> -c ACDE -r Sample_input/Quick_Start/5isz.pdb

## Analysis
To conduct an analysis of a single docked TCRpMHC complex, we will begin with both the Rosetta Energy Breakdown Table and Heatmap. Energy Breakdown will provide us a summary of which amino acids are interacting between the TCR and the pMHC.

    cd ../tr05analysis
    
    # Table (output = output.csv)
    DataDepot -t Sample_Input/Quick_Start/ES179M1-01_5isz_Docked.pdb
    
    DataDepot -t Sample_Input/Quick_Start/ES179M1-01_5isz_Docked.pdb -m
    
    # Heatmap
    DataDepot -e Sample_Input/Quick_Start/ES179M1-01_5isz_Docked.pdb
    
    DataDepot -e Sample_Input/Quick_Start/ES179M1-01_5isz_Docked.pdb -m
ABusage calculates interface scores for alpha chain to pMHC and beta chain to pMHC, utilizing the Rosetta application InterfaceAnalyzer.

    DataDepot -a Sample_Input/Quick_Start/ES179M1-01_5isz_Docked.pdb


--------------------------------------------------------------------------------

Please send an email to Austin Seamann or Dario Ghersi with any questions
regarding TRain. If having difficulties with installing dependencies, please
visit their installation/support pages before reaching out.
(aseamann or dghersi [at] unomaha.edu)


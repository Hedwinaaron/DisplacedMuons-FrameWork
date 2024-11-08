# Displaced Muons Framework

Documentation for the Ntuplizer is from https://github.com/CeliaFernandez/standard-Ntuplizer/tree/main
## How to install

Recommended release for this analyzer is CMSSW_13_0_13 or later. Commands to setup the analyzer are:

```
cmsrel CMSSW_13_0_13

cd CMSSW_12_4_0/src

cmsenv

git clone git@github.com:CeliaFernandez/standard-Ntuplizer.git

scram b -j 8
```
## How to Run 

The Ntuplizer can be run with the setup configuration through the cfg file:

```
cmsRun test/runNtuplizer_cfg.py
```

This requires a path to a valid MiniAOD sample in `runNtuplizer_cfg.py`:
https://github.com/Hedwinaaron/DisplacedMuons-FrameWork/blob/6085a81c85cee6c76bd7ba265f1cf9847bed2a27/Ntuplizer/test/Cosmics_runNtuplizer_MiniAOD_cfg.py#L32

Alternatively, the code can be run using CRAB 3. For this, make sure that the `crab_CosmicsAnalysis_Run2023_MiniAOD.py` file is configured properly:
1. Set a valid working area:
https://github.com/Hedwinaaron/DisplacedMuons-FrameWork/blob/6085a81c85cee6c76bd7ba265f1cf9847bed2a27/Ntuplizer/test/crab_CosmicsAnalysis_Run2023_MiniAOD.py#L9
2. Set a valid dataset:
https://github.com/Hedwinaaron/DisplacedMuons-FrameWork/blob/6085a81c85cee6c76bd7ba265f1cf9847bed2a27/Ntuplizer/test/crab_CosmicsAnalysis_Run2023_MiniAOD.py#L9
   The listed dataset is valid in the GRID and can be used to run the Ntuplizer.
3. Set a valid storage address:
https://github.com/Hedwinaaron/DisplacedMuons-FrameWork/blob/6085a81c85cee6c76bd7ba265f1cf9847bed2a27/Ntuplizer/test/crab_CosmicsAnalysis_Run2023_MiniAOD.py#L27
Make sure that you have a valid proxy before running and do at least once:

```
voms-proxy-init --voms cms
```
Then run  `crab_CosmicsAnalysis_Run2023_MiniAOD.py` as:
```
crab submit crab_CosmicsAnalysis_Run2023_MiniAOD
```
## Plots
The resulting ntuples can be used to generate plots using a basic tag and probe selection by running the `Plots/Dmuons_plots_tnp.C` script. First, provide a valid path to the root files:
https://github.com/Hedwinaaron/DisplacedMuons-FrameWork/blob/6085a81c85cee6c76bd7ba265f1cf9847bed2a27/Plots/Dmuons_plots_tnp.C#L21

Then run the script using:
```
root -q -l Dmuons_plots_tnp.C
```

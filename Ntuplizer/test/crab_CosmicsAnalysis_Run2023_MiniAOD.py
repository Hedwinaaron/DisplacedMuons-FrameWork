import CRABClient
from CRABClient.UserUtilities import config, getLumiListInValidFiles
#from FWCore.DataStructs.LumiList import LumiList
from FWCore.PythonUtilities.LumiList import LumiList

config = config()

# General
config.General.workArea = '/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023'
config.General.requestName = 'CosmicsAnalysis_Run2023_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_9_22_test5'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Cosmics_runNtuplizer_MiniAOD_cfg.py'
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['ntuples.root']

# Data
config.Data.inputDataset = '/NoBPTX/Run2023C-22Sep2023_v4-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023'
config.Data.publication = False
config.Data.outputDatasetTag = 'CosmicsAnalysis_Run2022F_MiniAOD-Ntuples_TnP_CMSSW_13_2_0_pre1_9_23_test5'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'

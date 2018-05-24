from WMCore.Configuration import Configuration
config = Configuration()
#General
config.section_('General')
config.General.requestName  = 'job_Run2016B-23Sep2016_ReReco-v3'
config.General.workArea     = 'ReReco-BCDEFG_PromptReco-H'
config.General.transferLogs = True
#JobType
config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles  = ['Spring16_23Sep2016AllV2_DATA.db', 'Spring16_25nsV10p2_DATA_L2Relative_AK8PFchs.txt', 'Spring16_25nsV10p2_DATA_L3Absolute_AK8PFchs.txt', 'Spring16_25nsV10p2_DATA_L2L3Residual_AK8PFchs.txt']
config.JobType.psetName    = 'run_data_Run2016B-23Sep2016_ReReco-v3.py'
config.JobType.outputFiles = ['Data_Run2016B-23Sep2016_ReReco-v3.root']
config.JobType.sendExternalFolder = True
#Data
config.section_('Data')
config.Data.inputDataset = '/SinglePhoton/Run2016B-23Sep2016-v3/MINIAOD'
config.Data.lumiMask     = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 10
config.Data.publication  = False
config.Data.totalUnits   = -1
#config.Data.ignoreLocality = True
config.Data.outLFNDirBase  ='/store/user/rgarg/13TeV/Ntuples/80X/Data/ReReco-BCDEFG_PromptReco-H/Run2016B-23Sep2016_ReReco-v3'
#User
config.section_('User')
#Site
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'

#!/usr/bin/python
#------------------------------------------------------------------------------
# rmci1 family: regeneration of the rmci0 with the fixed generator module
#------------------------------------------------------------------------------
from local_classes import *
from mixing_inputs import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name]            = Stage(name,self);
        return self.fStage[name]

    def __init__(self):
        project                      = 'Mu2eII'
        familyID                     = 'rmci1'

        self.fProjectName            = project
        self.fDsid                   = familyID
        self.fStage                  = {}
#------------------------------------------------------------------------------
# s4:sim : zero luminosity dataset
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator', 'rmci1s00b0', 'local')
        job.fNInputFiles             = 100                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  100000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s4:sim_b1 : 1-batch mode dataset
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_b1_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','rmci1s00b1','local');
        job.fNInputFiles             = 500
        define_mixing_inputs(job);

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       =  20000                   # for 1-batch mode signal+minbias events
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# init s4:concat concatenate s41b0
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_'+familyID+'.fcl'

        dsn                          = 'dig.mu2e.rmci1s41b0.Mu2eII.art'   # dataset: mu2e.rmci1s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'rmci1s41b0','local')  # dataset: mu2e.rmci1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  10                                 # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s5:reco_stn  : reconstruct and stntuple zero-luminosity dataset
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.rmci1s41b0.Mu2eII.art'   # dataset: mu2e.rmci1s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(defname,'','local')        # dataset: mu2e.rmci1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                            # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s5:reco_stn_b1  : reconstruct and stntuple 1-batch mode dataset, don't write output, only stntuple's
#------------------------------------------------------------------------------
        job                          = s.new_job('reco_stn_b1');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_b1_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.rmci1s41b1.Mu2eII.art'  # files in PNFS, not on tape yet
        job.fInputDataset            = Dataset(defname,'','local')       # dataset: mu2e.rmci1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  2
        job.fNEventsPerSegment       =  -1                               # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'        ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

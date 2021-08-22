#!/usr/bin/python

from local_classes import *

class Project:
    #------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
    def new_stage(self,name):
        self.fStage[name] = Stage(name);
        return self.fStage[name]

    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'bmup2'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s1:sim
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s1');
        job                        = s.new_job('sim');

        job.fRunNumber             = 1000;
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDataset          = Dataset('generator','bmup2s00b0','local');
        job.fNInputFiles           = 500                      # for an event generator - N(job segments) 

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment     =  200000
        job.fResample              = 'no'                     # yes/no
        job.fRequestedTime         = '15h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath            = ['trigmubeam' ]
        job.fOutputStream          = ['mubeamout'  ]
        job.fOutputDsID            = ['bmup2s11b0' ]
        job.fOutputFnPattern       = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat          = [ 'art'       ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s2:sim
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s2');
        job                        = s.new_job('sim');

        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDataset          = Dataset('sim.mu2e.bmup2s11b0.Mu2eII.art','bmup2s11b0','local');
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 50
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath            = ['trigmubeam' ]
        job.fOutputStream          = ['mubeamout'  ]
        job.fOutputDsID            = ['bmup2s21b0' ]
        job.fOutputFnPattern       = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat          = [ 'art'       ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s3:sim : of interest: stopped muons
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s3');
        job                        = s.new_job('sim');

        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDataset          = Dataset('sim.mu2e.bmup2s21b0.Mu2eII.art','bmup2s21b0','local');
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 100
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '12h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath            = ['tgtStopOutput', 'ootStopOutput']
        job.fOutputStream          = ['tgtstops'     , 'ootstops'     ]
        job.fOutputDsID            = ['bmup2s31b0'   , 'bmup2s32b0'   ]
        job.fOutputFnPattern       = [ 'sim.mu2e.'+job.fOutputDsID[0], 'sim.mu2e.'+job.fOutputDsID[1] ]
        job.fOutputFormat          = [ 'art'                         , 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4: sim_digi
#------------------------------------------------------------------------------        
        s                          = Stage('s4');

        job                        = Job('sim_digi');
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_sim_digi_'+dsid+'.fcl'
        job.fInputStage            = 's3'
        job.fInputStream           = 'mothers'

        job.fInputDsID             = 'bmup2s33b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup2s33b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream          = ['mothers'      ]
        job.fOutputFnPattern       = ['mothers'      ]
        job.fOutputPath            = ['defaultOutput']
        job.fOutputDsID            = ['bmup2s41b0'   ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s5: reco_stn
#------------------------------------------------------------------------------        
        s                          = Stage('s5');

        job                        = Job('reco_stn');
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
        job.fInputStage            = 's4'
        job.fInputStream           = 'sim_digi'  # this doesn't seem to be needed any more

        job.fInputDsID             = 'bmup2s41b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup2s41b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream          = ['mothers'      ]   # I think this is needed only for output integrity checks
        job.fOutputFnPattern       = ['mothers'      ]
        job.fOutputPath            = ['defaultOutput']
        job.fOutputDsID            = ['bmup2s51b0'   ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

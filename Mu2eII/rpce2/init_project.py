#!/usr/bin/python

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name]           = Stage(name,self);
        return self.fStage[name]

    def __init__(self):

        project                      = 'Mu2eII'
        familyID                     = 'rpce2'

        self.fProjectName            = project
        self.fDsid                   = familyID
        self.fStage                  = {}
#------------------------------------------------------------------------------
# init s3 - stage3 job (corresponds to pbar stage 4)
#------------------------------------------------------------------------------
        s                            = self.new_stage('s3');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s3_sim_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.rpce2s21b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  -1
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'tgtStopOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'root'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim (corresponds to pbar stage 5)
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_'+familyID+'.fcl'
#        job.fInputStage              = 's0'

        job.fInputDataset            = Dataset('generator','rpce2s31b0','local');
        job.fNInputFiles             = 50                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  1000000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s5:reco_stn
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'
#        job.fInputStage              = 's4'

        defname                      = 'dig.mu2e.rpce2s41b0.Mu2eII.art'  # 
        job.fInputDataset            = Dataset(defname,'rpce2s41b0','local') # dataset: mu2e.rpce2s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

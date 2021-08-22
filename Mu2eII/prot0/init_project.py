#!/usr/bin/python

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name]                = Stage(name,self);
        return self.fStage[name]

    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'prot0'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:sim 
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_gen_sim_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('generator','prot0s00b0','local');
        job.fNInputFiles             = 500                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  1000000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:concat
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'
        job.fInputStage              = 's4'

        job.fInputDsID               = 'prot0s41b0'                       # concatenation
        dsn                          = project+'.'+job.fInputDsID+'.art'  # dataset: mu2e.prot0s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'prot0s41b0','local')  # dataset: mu2e.prot0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  5                                 # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

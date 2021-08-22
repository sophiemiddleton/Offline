#!/usr/bin/python

from local_classes import *
from mixing_inputs import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name]                = Stage(name,self);
        return self.fStage[name]

    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'fele1'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# init s4; stage can have one or several jobs associated with it
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
#------------------------------------------------------------------------------
# s4:sim 
#------------------------------------------------------------------------------
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_gen_sim_digi_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('generator','fele1s00b0','local');
        job.fNInputFiles             = 50                  # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  20000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim_b1 : 1-batch mode electron dataset, no concatenation needed, 1M = 500 (segments)x2000(events per segment)
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_gen_sim_digi_b1_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('generator','fele1s00b1','local');
        job.fNInputFiles             = 500                      # number of the job segments for 1M events total

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       =  2000                    # for 1-batch mode signal+minbias events generate 2000 events/segment
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b1'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
 
        define_mixing_inputs(job)
#------------------------------------------------------------------------------
# s4:sim_b2 : 2-batch mode electron dataset, no concatenation needed, 1M = 1000 (segments)x1000(events per segment)
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b2');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_gen_sim_digi_b2_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('generator','fele1s00b2','local');
        job.fNInputFiles             = 500                      # number of the job segments for 1M events total

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       =  2000                    # for 1-batch mode signal+minbias events generate 2000 events/segment
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b2'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
 
        define_mixing_inputs(job)
#------------------------------------------------------------------------------
# s4:concat
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'
        job.fInputStage              = 's4'

        defname                      = 'sim.mu2e.fele1s41b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'fele1s41b0','local') # dataset: mu2e.fele1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                    # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  5                                    # 
        job.fNEventsPerSegment       =  20000                                # not used for concatenation
        job.fResample                = 'no'                                  # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s5: reco_stn
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        defname                      = 'sim.mu2e.fele1s41b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'fele1s41b0','local') # dataset: mu2e.fele1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                    # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

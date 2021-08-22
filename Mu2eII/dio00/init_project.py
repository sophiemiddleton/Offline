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
        familyID                             = 'dio00'

        self.fProjectName                = project
        self.fDsid                       = familyID
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:sim 
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','dio00s00b0','local');
        job.fNInputFiles             = 500                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  500000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim_b1 : 1-batch mode DIO dataset, no concatenation needed, 1M = 500 (segments)x2000(events per segment)
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_b1_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','dio00s00b1','local');
        job.fNInputFiles             = 500                      # number of the job segments for 1M events total

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       =  2000                    # for 1-batch mode signal+minbias events generate 2000 events/segment
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
 
        define_mixing_inputs(job)
#------------------------------------------------------------------------------
# s4:sim_b2 : 1-batch mode DIO dataset, no concatenation needed, 1M = 500 (segments)x2000(events per segment)
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b2');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_b2_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','dio00s00b2','local');
        job.fNInputFiles             = 500                      # number of the job segments for 1M events total

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       =  2000                    # for 1-batch mode signal+minbias events generate 2000 events/segment
        job.fResample                = 'no'                     # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b2'            ]
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
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_'+familyID+'.fcl'
        job.fInputStage              = 's4'

        job.fInputDataset            = Dataset('sim.mu2e.dio00s41b0.Mu2eII.art','dio00s41b0','local') 
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  5                                 # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:gen_sim_crv : input for CRV mixing, dsid=dio00s42b0
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('gen_sim_crv');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_crv_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','dio00s00b0','local');
        job.fNInputFiles             = 200                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  500000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'2b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:concat_crv  - for dio00s42b0, concatenate by x2
#------------------------------------------------------------------------------
        job                          = s.new_job('concat_crv');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_crv_'+familyID+'.fcl'
        job.fInputStage              = 's4'

        job.fInputDataset            = Dataset('sim.mu2e.dio00s42b0.Mu2eII.art','dio00s42b0','local') 
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  2                                 # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'2b0'        ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

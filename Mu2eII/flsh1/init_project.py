#!/usr/bin/python

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name] = Stage(name);
        return self.fStage[name]

    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'flsh1'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s1:sim
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s1');
        job                        = s.new_job('sim');

        job.fRunNumber             = 1000;
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_beam_'+dsid+'.fcl'

        job.fInputDataset          = Dataset('generator','flsh1s00b0','local');
        job.fNInputFiles           = 25000

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment     =  10000;
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '15h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath            = ['trigmubeam'                    ]
        job.fOutputStream          = ['mubeamout'                     ]
        job.fOutputDsID            = ['flsh1s10b0'                    ]          # 
        job.fOutputFnPattern       = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat          = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s1:concat
#------------------------------------------------------------------------------        
        job                          = s.new_job('concat');

        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'
#        job.fInputStage              = 's0'

        job.fInputDataset            = Dataset('sim.mu2e.flsh1s10b0.Mu2eII.art','flsh1s10b0','local');
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment = 50
        job.fNEventsPerSegment       =  -1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '3h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = ['out'           ]
        job.fOutputStream            = ['defaultOutput' ]
        job.fOutputDsID              = ['flsh1s11b0'    ]          # 
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = ['art'          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s2:sim
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s2');
        job                          = s.new_job('sim');

        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_beam_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('sim.mu2e.flsh1s11b0.Mu2eII.art','flsh1s11b0','local');
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment =  1
        job.fResample                = 'no'                   # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputPath              = ['trigmubeam' ]
        job.fOutputStream            = ['mubeamout'  ]
        job.fOutputDsID              = ['flsh1s21b0' ]
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = ['art'       ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s3:sim
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s3');
        job                          = s.new_job('sim');

        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_resample_'+dsid+'.fcl'
        job.fRunNumber               = 1000;

        job.fInputDataset            = Dataset('sim.mu2e.flsh1s21b0.Mu2eII.art','flsh1s21b0','local');
        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment =  1
        job.fResample                = 'yes'   # yes/no
        job.fNEventsPerSegment       =  200000
        job.fRequestedTime           = '15h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputPath              = [ 'mothersOutput' ]
        job.fOutputStream            = [ 'mothers'       ]
        job.fOutputDsID              = [ 'flsh1s33b0'    ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'           ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_sim_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('sim.mu2e.flsh1s33b0.Mu2eII.art','flsh1s33b0','local');
        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment =  10
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputPath              = [ 'mothersOutput' ]
        job.fOutputStream            = [ 'mothers'       ]
        job.fOutputDsID              = [ 'flsh1s41b0'    ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'           ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s5:sim
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s5');
        job                          = s.new_job('sim');

        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_add_proton_time_map_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('sim.mu2e.flsh1s41b0.Mu2eII.art','flsh1s41b0','local');
        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment =  100
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath              = [ 'mothersOutput' ]
        job.fOutputStream            = [ 'mothers'       ]
        job.fOutputDsID              = [ 'flsh1s51b0'    ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'           ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

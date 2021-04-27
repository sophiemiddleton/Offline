#!/usr/bin/python

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
        familyID                     = 'crv01'

        self.fProjectName            = project
        self.fDsid                   = familyID
        self.fStage                  = {}
#------------------------------------------------------------------------------
# s2:ds_resampler
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s2');
        job                          = s.new_job('ds_resampler');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_ds_resampler_'+familyID+'.fcl'

        job.fInputDsID               = "crv01s11b0" ;
        defname                      = 'sim.mu2e.crv01s11b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'crv01s11b0','local')
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment = 1
        job.fNEventsPerSegment       = 500000   # per Yuri's suggestion
        job.fResample                = 'yes'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'tgtStopOutput'               , 'detectorOutput'               ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'       , familyID+s.name()+'2b0'        ]  # crv01s21b0
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0], 'sim.mu2e.'+job.fOutputDsID[1] ]
        job.fOutputFormat            = [ 'art'                         , 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s2:ps_resampler
#------------------------------------------------------------------------------        
        job                          = s.new_job('ps_resampler');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_ps_resampler_'+familyID+'.fcl'

        job.fInputDsID               = "crv01s12b0" ;
        job.fInputDataset            = Dataset('sim.mu2e.crv01s12b0.Mu2eII.art','crv01s12b0','local')
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment = 1
        job.fNEventsPerSegment       = 500000   # to be desided 
        job.fResample                = 'yes'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'detectorOutput'               ]
        job.fOutputDsID              = [ familyID+s.name()+'3b0'        ]  # crv01s21b0
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s2:concat_ds_resampler : concatenating the CRV stream 
#------------------------------------------------------------------------------
        job                          = s.new_job('concat_ds_resampler');

        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_ds_resampler_'+familyID+'.fcl'

        job.fInputDsID               = 'crv01s22b0'                       # concatenation
        dsn                          = 'sim.mu2e.'+job.fInputDsID+'.Mu2eII.art'  # dataset: mu2e.crv01s22b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'crv01s22b0','local')  # dataset: mu2e.crv01s22b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  100                               # 
        job.fNEventsPerSegment       =  20000                             # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '2h'
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
# s2:concat_ps_resampler - concatenating the CRV stream
#------------------------------------------------------------------------------
        job                          = s.new_job('concat_ps_resampler');

        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_ps_resampler_'+familyID+'.fcl'

        job.fInputDsID               = 'crv01s23b0'                       # concatenation
        dsn                          = 'sim.mu2e.'+job.fInputDsID+'.Mu2eII.art'  # dataset: mu2e.crv01s22b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'crv01s23b0','local')  # dataset: mu2e.crv01s22b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  1000                              # 
        job.fNEventsPerSegment       =  20000                             # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '2h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'3b0'        ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s3:add_proton_time_map_ds (also perform filtering - require at least one flash hit T>200 ns)
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s3');
        job                          = s.new_job('add_proton_time_map_ds');

        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_add_proton_time_map_ds_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('sim.mu2e.crv01s22b0.Mu2eII.art','crv01s22b0','local');
        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment =  500
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath              = [ 'out'                          ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ 'crv01s32b0'                   ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s3:add_proton_time_map_ps (also perform filtering - require at least one flash hit T>200 ns)
#------------------------------------------------------------------------------        
        job                          = s.new_job('add_proton_time_map_ps');

        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_add_proton_time_map_ps_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('sim.mu2e.crv01s23b0.Mu2eII.art','crv01s23b0','local');
        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment =  500
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath              = [ 'out'                          ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ 'crv01s33b0'                   ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim_b1 : 1 batch mode pileup 1M events = (1000 segments) x (1000 events/segment) 
# limited by the time/segment
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_'+job.name()+'_'+familyID+'.fcl'
        job.fInputDataset            = Dataset('generator','crv01s00b1','local');
        job.fNInputFiles             = 1000                              # 1000
        job.fNEventsPerSegment       = 1000
        job.fMaxInputFilesPerSegment = 1

        define_mixing_inputs_crv(job);

        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '16h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath              = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim_b2 : 2-batch mode pileup 500K events = (1000 segments) x (500 events/segment)
# limited by the output file size / time per segment
#------------------------------------------------------------------------------        
        job                          = s.new_job('sim_b2');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_'+job.name()+'_'+familyID+'.fcl'
        job.fInputDataset            = Dataset('generator','crv01s00b2','local');
        job.fNInputFiles             = 1000                              # 1000
        job.fNEventsPerSegment       = 500
        job.fMaxInputFilesPerSegment = 1

        define_mixing_inputs_crv(job);

        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '16h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath              = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b2'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

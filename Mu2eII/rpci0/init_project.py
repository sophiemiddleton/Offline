#!/usr/bin/python

from local_classes import *
from mixing_inputs import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name] = Stage(name);
        return self.fStage[name]

    def __init__(self):

        project                      = 'Mu2eII'
        familyID                     = 'rpci0'

        self.fProjectName            = project
        self.fDsid                   = familyID;
        self.fStage                  = {}
#------------------------------------------------------------------------------
# s4:sim ; 
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_'+familyID+'.fcl'
        job.fInputStage              = 's0'

        job.fInputDataset            = Dataset('generator','rpci0s00b0','local');
        job.fNInputFiles             = 100                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  10000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim_b1 : 1-batch mode CE dataset, no concatenation needed, 1M = 500 (segments)x2000(events per segment)
#             limit on the number of events/segment comes from the output file size
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_b1_'+familyID+'.fcl'

        dsid                         = familyID+'s00b1';
        job.fInputDataset            = Dataset('generator',dsid,'local');
        job.fNInputFiles             = 500                      # number of the job segments for 1M events total
        job.fNEventsPerSegment       = 2000                     # for 1-batch mode signal+minbias events generate 2000 events/segment
        define_mixing_inputs(job);

        job.fMaxInputFilesPerSegment =  1                       # =1 for MC generator
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
# s4:sim_b2 : 2-batch mode CE dataset, no concatenation needed , 1M = 1000(segments)x1000(events per segment)
#------------------------------------------------------------------------------
        job                          = s.new_job('sim_b2');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_b2_'+familyID+'.fcl'

        idsid                        = familyID+'s00b2';
        job.fInputDataset            = Dataset('generator',idsid,'local');
        job.fNInputFiles             = 100                                # number of the job segments for 1M events total
        job.fNEventsPerSegment       = 10000                              # for 1-batch mode signal+minbias events generate 2000 events/segment
        define_mixing_inputs(job);

        job.fMaxInputFilesPerSegment =  1                                 # MC generator
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b2'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:concat 
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_'+familyID+'.fcl'


        dsid                         = familyID+'s41b0';
        defname                      = 'dig.mu2e.'+dsid+'.Mu2eII.art'     # dataset: mu2e.rpci0s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(defname,'','local')        # dataset: mu2e.rpci0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  10                                # 
        job.fNEventsPerSegment       =  20000                             # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ] # # the same as the input FamilyID
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'
#------------------------------------------------------------------------------
# s5:reco_stn
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'

        dsid                         = familyID+'s41b0';
        defname                      = 'dig.mu2e.'+dsid+'.Mu2eII.art'    # input dataset SAM definition name
        job.fInputDataset            = Dataset(defname,'','local')       # dataset: 
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                          # ifdh/xrootd

        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
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

        dsid                         = familyID+'s41b1';
        defname                      = 'dig.mu2e.'+dsid+'.Mu2eII.art'    # input dataset SAM definition name
        job.fInputDataset            = Dataset(defname,'','local')       # dataset: mu2e.cele0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  5
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
# s5:reco_stn_b2  : reconstruct and stntuple 1-batch mode dataset, don't write output, only stntuple's
#------------------------------------------------------------------------------
        job                          = s.new_job('reco_stn_b2');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_b2_'+familyID+'.fcl'

        idsid                        = familyID+'s41b2';                 # input dsid
        defname                      = 'dig.mu2e.'+idsid+'.Mu2eII.art'   # input dataset SAM definition name
        job.fInputDataset            = Dataset(defname,'','local')       # dataset: mu2e.cele0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  5
        job.fNEventsPerSegment       =  -1                               # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                          # ifdh/xrootd

        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b2'        ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

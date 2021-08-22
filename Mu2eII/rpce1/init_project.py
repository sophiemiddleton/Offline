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
        dsid                         = 'rpce1'

        self.fProjectName            = project
        self.fDsid                   = dsid
        self.fStage                  = {}
#------------------------------------------------------------------------------
# s4:sim
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_gen_sim_digi_'+dsid+'.fcl'
#        job.fInputStage              = 's0'

        job.fInputDataset            = Dataset('generator','rpce1s00b0','local');
        job.fNInputFiles             = 100                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  1000000
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
# init s4:concat 
# simulation produces about 300 events per file, 9 MByte large files
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'
#        job.fInputStage              = 's4'

        job.fInputDsID               = 'rpce1s41b0'                       # concatenation
        dsn                          = project+'.'+job.fInputDsID+'.art'  # dataset: mu2e.rpce1s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'rpce1s41b0','local')  # dataset: mu2e.rpce1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  10                                # 
        job.fNEventsPerSegment       =  20000                             # not used for concatenation
        job.fResample                = 'no'                               # yes/no
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
        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'
#------------------------------------------------------------------------------
# s5:reco_stn
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
#        job.fInputStage              = 's4'

        dsn                          = 'dig.mu2e.rpce1s41b0.Mu2eII.art'  # 
        job.fInputDataset            = Dataset(dsn,'rpce1s41b0','local') # dataset: mu2e.rpce1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

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
        job.fOutputFormat            = [ 'art:stn'                      ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'
#------------------------------------------------------------------------------
# s6:reco_stn - reprocessing of s5:reco_stn
#------------------------------------------------------------------------------
        s                            = self.new_stage('s6');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
#        job.fInputStage              = 's6'

        dsn                          = 'mcs.mu2e.rpce1s51b0.Mu2eII.art'  # 
        job.fInputDataset            = Dataset(dsn,'rpce1s51b0','sam') # dataset: mu2e.rpce1s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                            # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                      ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

#!/usr/bin/python

import os;
from local_classes import *

#------------------------------------------------------------------------------
# not used currently
def get_number_of_files(dataset):
    
    if (dataset == 'Mu2eII.1000_1998.g4s4_digi.art') :
        nfiles = 10;

class Project:
    #------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'bpim0'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# init first stage. a Stage can have one or several jobs associated with it
#------------------------------------------------------------------------------        
        s                            = Stage('s1');

        job                          = Job('sim');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/s1_pion_beam_'+dsid+'.fcl'

        job.fInputDataset            = None;           # generator
        job.fInputDsID               = 'bpim0s00b0';   # 's0' indicates the generator
        job.fNInputFiles             = 400

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputStream            = [ 'mubeamout'  ]
        job.fOutputFnPattern         = [ 'PS-mubeam'  ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0' ]
        job.fOutputPath              = [ 'trigmubeam' ]

        # grid output dir
        desc = project+'.bpim0s00b0.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# init second stage
#------------------------------------------------------------------------------        
        s                            = Stage('s2');

        job                          = Job('sim');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_pion_beam_'+dsid+'.fcl'

        job.fInputDsID               = 'bpim0s11b0'
        dsn                          = 'Mu2eII.bpim0s11b0.art'
        job.fInputDataset            = Dataset(dsn,'bpim0s11b0','local')
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment =  20
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputStream            = ['mubeamout'  ]
        job.fOutputFnPattern         = ['TS-mubeam'  ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0' ]
        job.fOutputPath              = ['trigmubeam' ]
        
        # grid output dir
        job.fDescription             = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# bpim3 3rd stage uses input from the 2nd stage of bpim0
#------------------------------------------------------------------------------        
        s                            = Stage('s3');

        job                          = Job('sim');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_pion_beam_'+dsid+'.fcl'

        job.fInputDsID               = 'bpim0s21b0'
        dsn                          = 'Mu2eII.bpim0s21b0.art'
        job.fInputDataset            = Dataset(dsn,'bpim0s21b0','local')
        job.fNInputFiles             = -1      # to be figured from the input dataset

        job.fMaxInputFilesPerSegment =  1      # stage 3 job takes more time than stage2, no point in grouping the files
        job.fResample                = 'no'    # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = ['tgtstops'     , 'ootstops'     ]
        job.fOutputFnPattern         = ['DS-TGTstops'  , 'DS-OOTstops'  ]
        job.fOutputPath              = ['tgtStopOutput', 'ootStopOutput']
        job.fOutputDsID              = [ dsid+s.name()+'1b0'   , dsid+s.name()+'2b0'   ]

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# bpim3 3rd stage uses input from the 2nd stage of bpim0 'tnt': target stop ntuple
#------------------------------------------------------------------------------        
        s                            = Stage('s3');

        job                          = Job('tgt_nt');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_pion_tgtstop_ntuple_'+dsid+'.fcl'

        job.fInputDsID               = 'bpim0s31b0'
        dsn                          = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset            = Dataset(dsn,'bpim0s31b0','local')
        job.fNInputFiles             = -1                                 # to be figured from the input dataset

        job.fMaxInputFilesPerSegment =  500                               # 
        job.fResample                = 'no'                               # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

#         job.fOutputStream            = ['tgtstops'     , 'ootstops'     ]
#         job.fOutputFnPattern         = ['DS-TGTstops'  , 'DS-OOTstops'  ]
#         job.fOutputPath              = ['tgtStopOutput', 'ootStopOutput']
#         job.fOutputDsID              = [ dsid+s.name()+'1b0'   , dsid+s.name()+'2b0'   ]

        # grid output dir
        desc = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

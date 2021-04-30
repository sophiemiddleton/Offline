#!/usr/bin/python

from local_classes import *

#------------------------------------------------------------------------------
def get_number_of_files(dataset):
    
    if (dataset == 'Mu2eII.1000_1998.g4s4_digi.art') :
        nfiles = 10;

class Project:
    #------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'bmum3'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# init first stage. a Stage can have one or several jobs associated with it
#------------------------------------------------------------------------------        
        s                          = Stage('s1');

        job                        = Job('sim');
        job.fRunNumber             = 1000;
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'
        job.fInputStage            = 's0'
        job.fInputStream           = 'gen'

        job.fInputDsID             = 'bmum3s00b0'              # generator
        job.fInputDataset          = None;                     # generator
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment =  1
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '12h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd
        job.fOutputStreams         = [ 'mubeamout'         ]
        job.fOutputFnPattern       = [ 'PS-mubeam'         ]
        job.fOutputDsID            = [ dsid+s.name()+'1b0' ]
        job.fOutputPath            = [ 'trigmubeam'        ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# init second stage
#------------------------------------------------------------------------------        
        s                          = Stage('s2');

        job                        = Job('sim');
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'
        job.fInputStage            = 's1'
        job.fInputStream           = 'mubeam'

        job.fInputDsID             = 'bmum3s11b0'
        dsn                        = project+'.'+job.fInputDsID+'art';
        job.fInputDataset          = Dataset(dsn,'bmum3s11b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 1
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '5h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd
        job.fOutputStream          = ['mubeamout'          ]
        job.fOutputFnPattern       = ['TS-mubeam'          ]
        job.fOutputDsID            = [ dsid+s.name()+'1b0' ]
        job.fOutputPath            = ['trigmubeam'         ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# bmum3 3rd stage uses input from the 2nd stage of bmum0
#------------------------------------------------------------------------------        
        s                            = Stage('s3');

        job                          = Job('sim');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'
        job.fInputStage              = 's1'
        job.fInputStream             = 'mubeam'

        job.fInputDsID               = 'bmum3s21b0'
        dsn                          = project+'.'+job.fInputDsID+'art';
        job.fInputDataset            = Dataset(dsn,'bmum3s21b0','local')
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment =  2
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = ['tgtstops'          , 'ootstops'          ]
        job.fOutputFnPattern         = ['DS-TGTstops'       , 'DS-OOTstops'       ]
        job.fOutputPath              = ['tgtStopOutput'     , 'ootStopOutput'     ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0', dsid+s.name()+'2b0' ]

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

#!/usr/bin/python

from local_classes import *

class Project:
    #------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'bmum1'

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

        job.fInputDsID             = 'bmum1s00b0';
        job.fInputDataset          = None;               # generator, input dataset is not defined, DsID - is
        job.fNInputFiles           = -1

        job.fNInputFilesPerSegment =  1
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '12h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd
        job.fOutputStream          = [ 'mubeamout'  ]
        job.fOutputFnPattern       = [ 'PS-mubeam'  ]
        job.fOutputDsID            = [ dsid+s.name()+'1b0' ]
        job.fOutputPath            = [ 'trigmubeam' ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# bmum1 reuses s2 from bmum0, so the 's2' part is not really needed, 
#       it is just easier to keep it, then to delete
#------------------------------------------------------------------------------        
        s                          = Stage('s2');

        job                        = Job('sim');
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'
        job.fInputStage            = 's1'
        job.fInputStream           = 'mubeam'

        job.fInputDsID             = job.fInputDataset.split('.')[1];
        job.fInputDataset          = Dataset('Mu2eII.bmum1s11b0.art','bmum1s11b0','local')
        job.fNInputFiles           = -1

        job.fNInputFilesPerSegment =  1
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '5h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd
        job.fOutputStream          = ['mubeamout'  ]
        job.fOutputFnPattern       = ['TS-mubeam'  ]
        job.fOutputDsID            = [ dsid+s.name()+'1b0' ]
        job.fOutputPath            = ['trigmubeam' ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# s3:sim 
#              3rd stage uses input from the 2nd stage of bmum0
#------------------------------------------------------------------------------        
        s                          = Stage('s3');

        job                        = Job('sim');
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'
        job.fInputStage            = 's1'
        job.fInputStream           = 'mubeam'

        job.fInputDataset          = Dataset('Mu2eII.bmum1s21b0.art','bmum1s21b0','local')
        job.fInputDsID               = job.fInputDataset.split('.')[1];
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment =  1
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '5h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream          = ['tgtstops'     , 'ootstops'     ]
        job.fOutputFnPattern       = ['DS-TGTstops'  , 'DS-OOTstops'  ]
        job.fOutputPath            = ['tgtStopOutput', 'ootStopOutput']
        job.fOutputDsID            = [ dsid+s.name()+'1b0'   , dsid+s.name()+'2b0'   ]

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# s4:tgt: make muon target stop ntuple
#------------------------------------------------------------------------------        
        s                          = Stage('s4');

        job                          = Job('tgt');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_tgtstop_ntuple_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('Mu2eII.bmum1s31b0.art','bmum1s31b0','local')
        job.fInputDsID               = job.fInputDataset.split('.')[1];
        job.fNInputFiles             = -1      # to be figured from the input dataset

        job.fMaxInputFilesPerSegment =  500    # stage 3 job takes more time than stage2, no point in grouping the files
        job.fResample                = 'no'    # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# s4:oot: make muon out-of-target stop ntuple
#------------------------------------------------------------------------------        
        job                          = Job('oot');
        job.fTarball                 = None;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_ootstop_ntuple_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('Mu2eII.bmum1s32b0.art','bmum1s32b0','local')
        job.fInputDsID               = job.fInputDataset.split('.')[1];
        job.fNInputFiles             = -1      # to be figured from the input dataset

        job.fMaxInputFilesPerSegment =  500    # stage 3 job takes more time than stage2, no point in grouping the files
        job.fResample                = 'no'    # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]           = job
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

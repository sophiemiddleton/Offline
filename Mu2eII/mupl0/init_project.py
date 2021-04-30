#!/usr/bin/python

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'mupl0'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
        #------------------------------------------------------------------------------
        # init s4; stage can have one or several jobs associated with it
        s                            = Stage('s4');

        #------------------------------------------------------------------------------
        # init s4:sim ; 

        job                          = Job('sim');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_gen_sim_digi_'+dsid+'.fcl'
        job.fInputStage              = 's0'

        job.fInputDataset            = None;
        job.fInputDsID               = 'mupl0s00b0'        # muon stops
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
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        s.fJob[job.name()]         = job
        #------------------------------------------------------------------------------
        # init s4:concat ; 

        job                          = Job('concat');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'
        job.fInputStage              = 's4'

        job.fInputDsID               = 'mupl0s41b0'                       # concatenation
        dsn                          = project+'.'+job.fInputDsID+'.art'  # dataset: mu2e.mupl0s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'mupl0s41b0','local')  # dataset: mu2e.mupl0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  5                                 # 
        job.fNEventsPerSegment       =  20000                             # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job

        self.fStage[s.name()]      = s;
        #------------------------------------------------------------------------------
        # init first stage. a Stage can have one or several jobs associated with it
        s                            = Stage('s5');

        job                          = Job('reco_stn');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
        job.fInputStage              = 's4'

        job.fInputDsID               = 'mupl0s41b0'                      # dsID
        dsn                          = project+'.'+job.fInputDsID+'.art' # dataset: mu2e.mupl0s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'mupl0s41b0','local') # dataset: mu2e.mupl0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                            # 1 input file, IFDH is more robust ... ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]

        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

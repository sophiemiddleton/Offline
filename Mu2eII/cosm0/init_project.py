#!/usr/bin/python

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def __init__(self):
        project                          = 'Mu2eII'
        dsid                             = 'cosm0'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:
#------------------------------------------------------------------------------
        s                            = Stage('s4');
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# s4:helix_filter ; 
#------------------------------------------------------------------------------
        job                          = Job('helix_filter');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/s4_helix_filter_'+dsid+'.fcl'

        job.fInputDsID               = 'cosm0s00b0'        # 
        job.fInputDataset            = Dataset('dig.mu2e.CRY-cosmic-general.2025lo.art','cosm0s00b0','sam');
        job.fNInputFiles             = -1                  # defined by the dataset

        job.fMaxInputFilesPerSegment = 50                  # 100 files per segment
        job.fNEventsPerSegment       = -1                  # 
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '5h'                # normally, it should be fast
        job.fIfdh                    = 'xrootd'            # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        s.add_job(job)
#------------------------------------------------------------------------------
# s4:concat ; 
#------------------------------------------------------------------------------
        job                          = Job('concat');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'

        job.fInputDsID               = 'cosm0s41b0'                       # concatenation
        dsn                          = project+'.'+job.fInputDsID+'.art'  # dataset: mu2e.cosm0s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'cosm0s41b0','local'); # dataset: 
        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  5                                 # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
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
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# s5: reco_stn
#------------------------------------------------------------------------------
        s                            = Stage('s5');
        self.fStage[s.name()]        = s;

        job                          = Job('reco_stn');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
        job.fInputStage              = 's4'

        job.fInputDataset            = Dataset('dig.mu2e.cosm0s41b0.Mu2eII.art','cosm0s41b0','sam') # dataset: mu2e.cosm0s41b0.Mu2eII.art
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
        job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                      ]

        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# s5: reco_stn_00, reco_stn_01, reco_stn_02
#     2000 events per segment... total ~5000 events/file
#     3 jobs trun on the same input file, process different parts 
#------------------------------------------------------------------------------
        segments = ['00','01','02']
        for seg in segments:
            job_name                     = 'reco_stn_'+seg;
            job                          = Job(job_name);
            job.fRunNumber               = 1000;
            job.fBaseFcl                 = project+'/'+dsid+'/s5_'+job_name+'_'+dsid+'.fcl'
            
            job.fInputDsID               = 'cosm0s41b0'                      # dsID
            job.fInputDataset            = Dataset('dig.mu2e.cosm0s41b0.Mu2eII.art','cosm0s41b0','sam') # dataset: mu2e.cosm0s41b0.Mu2eII.art
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
            job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
            job.fOutputFormat            = [ 'art:stn'                      ]

            # grid output dir
            desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
            job.fDescription             = desc;
            
            # directory where output is saved from scratch dcache
            job.fOutputTopDir            = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'
            
            s.add_job(job)
#------------------------------------------------------------------------------
# s6: re-reconstruction of s5
#------------------------------------------------------------------------------
        s                            = Stage('s6');
        self.fStage[s.name()]        = s;

        job                          = Job('reco_stn');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDsID               = 'cosm0s51b0'                      # dsID
        dsn                          = 'mcs.mu2e.cosm0s51b0.Mu2eII.art'  # 
        job.fInputDataset            = Dataset(dsn,'cosm0s51b0','sam')   # dataset: mu2e.cosm0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                          # ifdh/xrootd
        job.fOutputPath              = [ 'out'                          ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                      ]

        # grid output dir
        desc                         = project+'.'+job.fInputDsID+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# s7: re-reconstruction of s6
#------------------------------------------------------------------------------
        s                            = Stage('s7');
        self.fStage[s.name()]        = s;

        job                          = Job('reco_stn');
        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('mcs.mu2e.cosm0s61b0.Mu2eII.art','cosm0s61b0','local') 
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
        job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                      ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# s8: re-reconstruction of s7 with the '2025 light yield (no 30 "safety factor")
#------------------------------------------------------------------------------
        s                            = Stage('s8');
        self.fStage[s.name()]        = s;

        job                          = Job('reco_stn');
        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('mcs.mu2e.cosm0s71b0.Mu2eII.art','cosm0s71b0','sam') 
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'                              # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                          # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                      ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# s9:remake_digis : use just one, cosmic, time map
#------------------------------------------------------------------------------
        s                            = Stage('s9');
        self.fStage[s.name()]        = s;

        job                          = Job('remake_digis');
        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_remake_digis_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('mcs.mu2e.cosm0s81b0.Mu2eII.art','cosm0s81b0','sam') 
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'                              # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                          # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# s9:reco_stn : re-reconstruct with correct timing
#------------------------------------------------------------------------------
        job                          = Job('reco_stn');
        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('dig.mu2e.cosm0s91b0.Mu2eII.art','cosm0s91b0','local') 
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'                              # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                            # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'msc.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir            = '/mu2e/data/users/sophie/datasets'

        s.add_job(job)
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

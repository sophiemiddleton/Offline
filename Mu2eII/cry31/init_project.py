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
        dsid                             = 'cry31'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:helix_filter : strip events with at least one helix
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('helix_filter');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/s4_helix_filter_'+dsid+'.fcl'

#        job.fInputDsID               = 'cry31s00b0'        # 
        job.fInputDataset            = Dataset('dig.mu2e.CRY-cosmic-general.2025hi.art','cry31s00b0','sam');

        job.fInputDataset.add_fileset('cry31s01b0','dig.mu2e.cry31s01b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s02b0','dig.mu2e.cry31s02b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s03b0','dig.mu2e.cry31s03b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s04b0','dig.mu2e.cry31s04b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s05b0','dig.mu2e.cry31s05b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s06b0','dig.mu2e.cry31s06b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s07b0','dig.mu2e.cry31s07b0.Mu2eII.art')
        job.fInputDataset.add_fileset('cry31s08b0','dig.mu2e.cry31s08b0.Mu2eII.art')

        job.fNInputFiles             = -1                  # defined by the dataset
        job.fMaxInputFilesPerSegment = 50                  # 50 files per segment
        job.fNEventsPerSegment       = -1                  # 
        job.fResample                = 'no'                # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '5h'                # normally, it should be fast
        job.fIfdh                    = 'xrootd'            # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:concatenation
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_concat_'+dsid+'.fcl'

#        job.fInputDsID               = 'cry31s41b0'                       # concatenation
        dsn                          = project+'.cry31s41b0.art';
        job.fInputDataset            = Dataset(dsn,'cry31s41b0','local'); # dataset: mu2e.cry31s41b0.Mu2eII.art

        job.fNInputFiles             = 50                                 # the number of input files to be defined dynamically
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
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'
#------------------------------------------------------------------------------
# s5: reconstruction
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
        dsn                          = project+'.cry31s41b0.art'         # 
        job.fInputDataset            = Dataset(dsn,'cry31s41b0','local') # dataset: mu2e.cry31s41b0.Mu2eII.art
        job.fNInputFiles             = 400                               # placeholder, the real number is defined by the input dataset
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
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/'+os.getenv('USER')+'datasets'
#------------------------------------------------------------------------------
# s5: reco_stn_00, reco_stn_01, reco_stn_02, reco_stn_03
#     2000 events per segment... total 6400-7600 events/file
#     4 jobs trun on the same input file, process different parts 
#------------------------------------------------------------------------------
        segments = ['00','01','02','03']
        for seg in segments:
            job_name                     = 'reco_stn_'+seg;
            job                          = s.new_job(job_name); # 
            job.fRunNumber               = 1000;
            job.fBaseFcl                 = project+'/'+dsid+'/s5_'+job_name+'_'+dsid+'.fcl'
            
            job.fInputDataset            = Dataset('Mu2eII.cry31s41b0.art','cry31s41b0','local') # 
            job.fNInputFiles             = -1                                # placeholder, the real number defined by the input dataset
        
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
            job.fOutputTopDir            = '/mu2e/data/users/'+os.getenv('USER')+'/datasets'
#------------------------------------------------------------------------------
# s6: re-reconstruction of s5
#------------------------------------------------------------------------------
        s                            = self.new_stage('s6');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('Mu2eII.cry31s51b0.art','cry31s51b0','local') 
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

        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        # directory where output is saved from scratch dcache
        job.fOutputTopDir          = '/mu2e/data/users/sophie/datasets'
#------------------------------------------------------------------------------
# s7: re-reconstruction of s6
#------------------------------------------------------------------------------
        s                            = self.new_stage('s7');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('mcs.mu2e.cry31s61b0.Mu2eII.art','cry31s61b0','local') 
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
#------------------------------------------------------------------------------
# s8: re-reconstruction of s7 with the '2025 light yield (no 30 "safety factor")
#------------------------------------------------------------------------------
        s                            = self.new_stage('s8');

        job                          = s.new_job('reco_stn');
        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('mcs.mu2e.cry31s71b0.Mu2eII.art','cry31s71b0','sam') 
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
        s                            = self.new_stage('s9');

        job                          = s.new_job('remake_digis');
        job.fRunNumber               = 1000;                             # not used
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_remake_digis_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('mcs.mu2e.cry31s81b0.Mu2eII.art','cry31s81b0','sam') 
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

        job.fInputDataset            = Dataset('dig.mu2e.cry31s91b0.Mu2eII.art','cry31s91b0','local') 
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

        s.add_job(job)
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

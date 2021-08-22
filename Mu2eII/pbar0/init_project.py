#!/usr/bin/python3

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name] = Stage(name);
        return self.fStage[name]

    def __init__(self):

        project                      = 'Mu2eII'
        familyID                     = 'pbar0'

        self.fProjectName            = project
        self.fDsid                   = familyID
        self.fStage                  = {}
#------------------------------------------------------------------------------
# init s0 - stage0 job
#------------------------------------------------------------------------------
        s                            = self.new_stage('s0');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s0_vertex_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','pbar0s00b0','local');
        job.fNInputFiles             = 1000                # number of the job segments
        job.fNEventsPerSegment       = 2000000

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

#         job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'filteredOutput'               ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s1 - stage1 job
#------------------------------------------------------------------------------
        s                            = self.new_stage('s1');
        job                          = s.new_job('gen_sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s1_gen_sim_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s01b0.Mu2eII.art'        # muon stops
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  1000000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd
#        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'filteredOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s1 - stage1 job
#------------------------------------------------------------------------------
        job                          = s.new_job('flat_sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s1_flat_sim_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s01b0.Mu2eII.art'        # muon stops
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                 # number of the job segments

        job.fMaxInputFilesPerSegment =  1                  # MC generator
        job.fNEventsPerSegment       =  2000000
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd
#        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'filteredOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'0b0'         ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0]  ]
        job.fOutputFormat            = [ 'root'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s1:concat - do it by 100
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');

#        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s11b0.Mu2eII.art'        # muon stops
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  100                               # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s1:reprocess ReadAntiProtonSteps tree 
#------------------------------------------------------------------------------
        job                          = s.new_job('newtree');

#        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_newtree_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s11b0.Mu2eII.art'        # previous stage 1
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  1000                                 # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'2b0'        ] # one more than the original 
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'root'                         ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s2 - stage2 job
#------------------------------------------------------------------------------
        s                            = self.new_stage('s2');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s2_sim_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s11b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             =  -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  1                # resampling 
        job.fNEventsPerSegment       =  3000000
        job.fResample                = 'yes'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'filteredOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s2 - stage2 job (low resampling)
#------------------------------------------------------------------------------
        job                          = s.new_job('low');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s2_low_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s11b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             =  -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  1                # resampling 
        job.fNEventsPerSegment       =  300000
        job.fResample                = 'yes'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'filteredOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s2:concat - do it by 100
#------------------------------------------------------------------------------
        job                          = s.new_job('concat');

#        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s21b0.Mu2eII.art'        
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  100                               # 
        job.fNEventsPerSegment       =  -1                                # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'filteredOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s3 - stage3 job
#------------------------------------------------------------------------------
        s                            = self.new_stage('s3');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s3_sim_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s21b0.Mu2eII.art'      
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  100                  # MC generator
        job.fNEventsPerSegment       =  -1
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'filteredOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s3 - stage3 job (low resampling)
#------------------------------------------------------------------------------
        job                          = s.new_job('low');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s3_low_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s21b0.Mu2eII.art'      
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  100                  # MC generator
        job.fNEventsPerSegment       =  -1
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'filteredOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s4 - stage4 job
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s4_sim_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s31b0.Mu2eII.art'      
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  -1
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'tgtStopOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'root'                         ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# init s4 - stage4 job (low resampling)
#------------------------------------------------------------------------------
        job                          = s.new_job('low');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/s4_low_'+familyID+'.fcl'

        defname                      = 'sim.mu2e.pbar0s31b0.Mu2eII.art'      
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1               # number of the job segments

        job.fMaxInputFilesPerSegment =  10
        job.fNEventsPerSegment       =  -1
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

        job.fOutputStream            = [ 'tgtStopOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'sim.mu2e.'+job.fOutputDsID[0] ]
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
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_concat_'+familyID+'.fcl'
        job.fInputStage              = 's4'

        defname                      = 'sim.mu2e.pbar0s41b0.Mu2eII.art'        # muon stops
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                                 # the number of input files to be defined dynamically

        job.fMaxInputFilesPerSegment =  10                                # 
        job.fNEventsPerSegment       =  20000                             # not used for concatenation
        job.fResample                = 'no'                               # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '1h'
        job.fIfdh                    = 'xrootd'                           # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ] 
        job.fOutputDsID              = [ familyID+s.name()+'1b0'            ] # # the same as the input DsID
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------q|
# s5:sim : make digi datasets (for pbars the stage numbers are shifted by one)
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_gen_sim_digi_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator','pbar0s41b0','local');
        job.fNInputFiles             = 100                       # number of the job segments

        job.fMaxInputFilesPerSegment =  1                       # MC generator
        job.fNEventsPerSegment       =  50000
        job.fResample                = 'no'                     # yes/no
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
# s6:reco_stn (for pbars there are more intermediate stages than for other datasets)
#------------------------------------------------------------------------------
        s                            = self.new_stage('s6');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.pbar0s51b0.Mu2eII.art'        # muon stops
        job.fInputDataset            = Dataset(defname,'','local')

        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]

        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

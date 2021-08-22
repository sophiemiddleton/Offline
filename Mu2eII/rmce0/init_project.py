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
        project                          = 'Mu2eII'
        familyID                         = 'rmce0'

        self.fProjectName                = project
        self.fDsid                       = familyID
        self.fStage                      = {}
#----------------------------------------------------------------------------------------
# s4:sim : generate weighted photons from muons stop to create a conversion map
#----------------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_weighted_photons_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator', 'rmce0s00b0', 'local');
        job.fNInputFiles             =  60 # need ~5e9 generated photons, 60 jobs = 10%

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  10000000 
        job.fResample                = 'no'                         # yes/no
        job.fMaxMemory               = '2500MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                     # ifdh/xrootd
        job.fOutputPath              = [ ]

        job.fOutputStream            = [ ]
        job.fOutputDsID              = [ ]
        job.fOutputFnPattern         = [ ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s5:sim : zero luminosity generate conversions from conversion map
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s5');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_conversions_gen_sim_dig_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator', 'rmce0s40b0', 'local')
        job.fNInputFiles             = 100 #need ~7e6-3e7 total events, so 10% ~50-100 jobs

        job.fMaxInputFilesPerSegment = 1
        job.fNEventsPerSegment       = 15000 # ~1.5 minutes / 1000 on mu2egpmv02, but next stage is very slow 
        job.fMaxMemory               = '2500MB' # default conversion map is large currently
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '3h'
        job.fIfdh                    = 'xrootd' # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDataset.id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s5:sim_b1 : 1-batch mode, generate conversions from conversion map
#------------------------------------------------------------------------------        
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_conversions_gen_sim_dig_b1_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator', 'rmce0s40b1', 'local');
        job.fNInputFiles             = 1000 #1.5M total events
        define_mixing_inputs(job);

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  1500 #
        job.fResample                = 'no'                         # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '15h'
        job.fIfdh                    = 'xrootd'                     # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.fInputDataset.id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s6:reco_stn: reconstruct output from fifth stage
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s6');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.rmce0s51b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'','local')
        job.fNInputFiles             = -1

        job.fMaxInputFilesPerSegment = 1
        job.fNEventsPerSegment       = 250000              # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out'                          ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s6:reco_stn_b1  : reconstruct and stntuple 1-batch mode dataset, don't write output, only stntuple's
#------------------------------------------------------------------------------
        job                          = s.new_job('reco_stn_b1');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_b1_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.rmce0s51b1.Mu2eII.art'  # files in PNFS, not on tape yet
        job.fInputDataset            = Dataset(defname,'','local')       # 
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  2
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
# end
#------------------------------------------------------------------------------

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
        familyID                         = 'rmce2'

        self.fProjectName                = project
        self.fDsid                       = familyID
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:sim : zero luminosity, generate flat photons from muon stops
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_flat_photons_gen_sim_dig_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator', 'rmce2s00b0', 'local');
        job.fNInputFiles             =  100    

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  300000 #~0.03 s/event in interactive, ~0.1% survive filtering
        job.fResample                = 'no'                         # yes/no
        job.fMaxMemory               = '2500MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                     # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.fInputDataset.id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s4:sim_b1 : 1-batch mode, generate flat photons from muon stops
#------------------------------------------------------------------------------        
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_flat_photons_gen_sim_dig_b1_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('generator', 'rmce2s00b1', 'local');
        job.fNInputFiles             = 100
        define_mixing_inputs(job);

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  500000 #~0.01 s/event on profile build, ~0.1% survive filtering
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
# s5:reco_stn: zero luminosity, reconstruct output from fourth stage
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s5');
        job                          = s.new_job('reco_stn');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.rmce2s41b0.Mu2eII.art'
        job.fInputDataset            = Dataset(defname,'','local')
        job.fNInputFiles             = -1                 

        job.fMaxInputFilesPerSegment = 10                  # takes ~4 s / event, stage 4 is ~0.1% efficient
        job.fNEventsPerSegment       = 250000              # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.fInputDataset.id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# s5:reco_stn_b1  : reconstruct and stntuple 1-batch mode dataset, don't write output, only stntuple's
#------------------------------------------------------------------------------
        job                          = s.new_job('reco_stn_b1');
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_b1_'+familyID+'.fcl'

        defname                      = 'dig.mu2e.rmce2s41b1.Mu2eII.art'  # files in PNFS, not on tape yet
        job.fInputDataset            = Dataset(defname,'','local')       # dataset: mu2e.cele0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset

        job.fMaxInputFilesPerSegment =  2 # ~2s / event interactive, stage 4 is 0.2% efficient
        job.fNEventsPerSegment       =  -1                               # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'        ]
        job.fOutputFnPattern         = [ 'nts.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]

        # grid output dir
        desc                         = project+'.'+job.fInputDataset.id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

#!/usr/bin/python

import os;
from local_classes import *

class Project:
    def new_stage(self,name):
        self.fStage[name] = Stage(name);
        return self.fStage[name]
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'rmce5'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# init fourth stage, generate flat photons from muon stops
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s4');

        job                          = s.new_job('sim');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_flat_photons_gen_sim_dig_'+dsid+'.fcl'
        job.fRunNumber               = 1000
        job.fInputStage              = 's0'

        job.fInputDataset            = Dataset('generator', 'rmce5s00b0', 'local'); # generator
        job.fNInputFiles             =  1000    #need ~1e9 total, 1e3 is 10% 

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  300000 #~0.03 s/event in interactive, ~0.1% survive filtering
        job.fResample                = 'no'                         # yes/no
        job.fMaxMemory               = '2500MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                     # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput' ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0']
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]

        # grid output dir
        desc                         = project+'.'+dsid+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# init fifth stage: reconstruct output from fourth stage
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s5');

        job                          = s.new_job('reco_stn');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
        job.fRunNumber               = 1000;
        job.fInputStage              = 's4'

        job.fInputDsID               = 'rmce5s41b0'        # conversions
        dsn                          = 'dig.mu2e.'+job.fInputDsID+'.Mu2eII.art' # dataset: dig.mu2e.rmce5s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'rmce5s41b0','local')      
        job.fNInputFiles             = -1                                       # placeholder, defined by input dataset

        job.fMaxInputFilesPerSegment = 10                  # takes ~3 s / event, stage 4 is ~0.1% efficient
        job.fNEventsPerSegment       = 250000              # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput' ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0']
        job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+dsid+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

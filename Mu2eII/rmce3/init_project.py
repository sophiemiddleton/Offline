#!/usr/bin/python

import os;
from local_classes import *

class Project:
    #------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
    def __init__(self):

        project                          = 'Mu2eII'
        dsid                             = 'rmce3'

        self.fProjectName                = project
        self.fDsid                       = dsid
        self.fStage                      = {}
#------------------------------------------------------------------------------
# init fourth stage, generate flat photons from muon stops
#------------------------------------------------------------------------------        
        s                            = Stage('s4');

        job                          = Job('sim');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_flat_photons_gen_sim_dig_'+dsid+'.fcl'
        job.fRunNumber               = 1000
        job.fInputStage              = 's0'

        job.fInputDataset            = None;               # generator
        job.fInputDsID               = 'rmce3s00b0'        # muon stops
        job.fNInputFiles             =  100    

        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  50000 #~0.2 s/event in interactive, ~10% survive filtering
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

        s.fJob[job.name()]           = job
        self.fStage[s.name()]        = s;
#------------------------------------------------------------------------------
# init fifth stage: reconstruct output from fourth stage
#------------------------------------------------------------------------------        
        s                            = Stage('s5');

        job                          = Job('reco_stn');
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'
        job.fRunNumber               = 1000;

        job.fInputStage              = 's4'

        job.fInputDsID               = 'rmce3s41b0'        # conversions
        dsn                          = project+'.mu2e.'+job.fInputDsID+'.art' # dataset: mu2e.fele0s41b0.Mu2eII.art
        job.fInputDataset            = Dataset(dsn,'rmce3s41b0','local')      # dataset: mu2e.fele0s41b0.Mu2eII.art
        job.fNInputFiles             = -1                 

        job.fMaxInputFilesPerSegment = 5                  # takes ~.2 s / event, stage 4 is ~10% efficient
        job.fNEventsPerSegment       = 250000              # placeholder
        job.fResample                = 'no'   # yes/no
        job.fMaxMemory               = '3000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput' ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0']
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        
        # grid output dir
        desc                         = project+'.'+dsid+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        s.fJob[job.name()]         = job
        self.fStage[s.name()]      = s;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

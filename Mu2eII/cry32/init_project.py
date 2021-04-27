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
        familyID                             = 'cry32'

        self.fProjectName                = project
        self.fDsid                       = familyID
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:helix_filter : strip events with at least one helix
#------------------------------------------------------------------------------
        s                            = self.new_stage('s4');
        job                          = s.new_job('helix_filter');

        job.fRunNumber               = -1;
        job.fBaseFcl                 = project+'/'+familyID+'/s4_helix_filter_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('dig.oksuzian.CRY-cosmic-general.cry3-digi-lo.art','cry32s00b0','local');

        job.fNInputFiles             = -1                  # defined by the dataset
        job.fMaxInputFilesPerSegment = 20                  # 20 files per segment
        job.fNEventsPerSegment       = -1                  # 
        job.fResample                = 'no'                # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '5h'                # normally, it should be fast
        job.fIfdh                    = 'xrootd'            # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s5:reco_stn : re-reconstruct with correct timing
#------------------------------------------------------------------------------
        s                            = self.new_stage('s5');

        job                          = s.new_job('reco_stn');
        job.fRunNumber               = -1;                               # not used
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_reco_stn_'+familyID+'.fcl'

        job.fInputDataset            = Dataset('dig.mu2e.cry32s41b0.Mu2eII.art','cry32s41b0','sam') 
        job.fNInputFiles             = -1                                # placeholder, the real number is defined by the input dataset
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       =  250000                           # placeholder
        job.fResample                = 'no'                              # yes/no
        job.fMaxMemory               = '2000MB'
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                          # ifdh/xrootd
        job.fOutputPath              = [ 'out' ]
        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b0'        ]
        job.fOutputFnPattern         = [ 'msc.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'stn'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dsid()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

        s.add_job(job)
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

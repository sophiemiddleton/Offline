#!/usr/bin/python

from local_classes import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly

    def new_stage(self,name):
        self.fStage[name]            = Stage(name,self);
        return self.fStage[name]

    def __init__(self):

        project                      = 'Mu2eII'
        dsid                         = 'bmup1'

        self.fProjectName            = project
        self.fDsid                   = dsid
        self.fStage                  = {}
#------------------------------------------------------------------------------
# init first stage. a Stage can have one or several jobs associated with it
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s1');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('generator','bmup1s00b0','local');
        job.fNInputFiles             = 4000                     # N(job segments)
        job.fNEventsPerSegment       = 250000

        job.fMaxInputFilesPerSegment =  1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'mubeamout'                    ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s2: trace up to TS5
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s2');

        job                        = s.new_job('sim');
        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDsID             = 'bmup1s11b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup1s11b0','local')

        job.fNInputFiles           = -1
        job.fMaxInputFilesPerSegment =  20
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '12h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'mubeamout'                    ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s3: make 'mothers'
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s3');
        job                        = s.new_job('sim');

        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDsID             = 'bmup1s21b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup1s21b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 5
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '12h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream          = ['tgtstops'          , 'ootstops'         , 'mothers'           ]
        job.fOutputDsID
        for i in range(0,3):
            job.fOutputDsID.append('bmup1'+s.name()+str(i+3)+'b0');
            job.fOutputFnPattern.append('dig.mu2e.'+job.fOutputDsID[i]);
            job.fOutputFormat.append('art')

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s4:sim_digi
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s4');
        job                        = s.new_job('sim_digi');

        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_sim_digi_'+dsid+'.fcl'

        job.fInputDataset          = Dataset('sim.mu2e.bmup1s33b0.Mu2eII.art','','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# s5:reco_stn
#------------------------------------------------------------------------------        
        s                          = self.new_stage('s5');
        job                        = s.new_job('reco_stn');

        job.fBaseFcl               = project+'/'+dsid+'/'+s.name()+'_reco_stn_'+dsid+'.fcl'

        job.fInputDsID             = 'bmup1s41b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup1s41b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'mcs.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art:stn'                      ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

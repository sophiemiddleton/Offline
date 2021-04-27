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
        dsid                         = 'bmup3'

        self.fProjectName            = project
        self.fDsid                   = dsid
        self.fStage                  = {}
#------------------------------------------------------------------------------
# s1:sim 
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s1');
        job                          = s.new_job('sim');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDataset            = Dataset('generator','bmup3s00b0','local');
        job.fNInputFiles             = 500                              # N(job segments)
        job.fNEventsPerSegment       = 1000000

        job.fMaxInputFilesPerSegment =  1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '24h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream            = [ 'mubeamout'                    ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

#------------------------------------------------------------------------------
# s2: trace up to TS5
        s                            = self.new_stage('s2')
        job                          = s.new_job('sim')
        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'

        job.fInputDsID               = 'bmup3s11b0'
        dsn                          = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset            = Dataset(dsn,'bmup3s11b0','local')

        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment = 5
        job.fResample                = 'no'  # yes/no
        job.fRequestedTime           = '12h'
        job.fIfdh                    = 'xrootd'          # ifdh/xrootd

        job.fOutputStream            = [ 'mubeamout'                    ]
        job.fOutputDsID              = [ dsid+s.name()+'1b0'            ]
        job.FOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# s2_resample: resample s1 output and trace into DS
#------------------------------------------------------------------------------        
        #s                            = self.new_stage('s2');

        job                          = s.new_job('sim_resample');
        #job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_muon_beam_'+dsid+'.fcl'
        job.fBaseFcl                 = project+'/'+dsid+'/'+s.name()+'_resample_'+dsid+'.fcl'

        job.fInputDsID               = 'bmup3s11b0'
        dsn                          = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset            = Dataset(dsn,'bmup3s11b0','local')

        job.fRunNumber               = 1000
        job.fNInputFiles             = -1
        job.fMaxInputFilesPerSegment =  1
        job.fResample                = 'yes'   # yes/no
        job.fNEventsPerSegment       =  1000000
        job.fRequestedTime           = '48h'
        job.fIfdh                    = 'ifdh'                 # ifdh/xrootd

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

        job.fInputDsID             = 'bmup3s21b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup3s21b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 10
        job.fResample              = 'no'   # yes/no
        job.fRequestedTime         = '12h'
        job.fIfdh                  = 'xrootd'                 # ifdh/xrootd

        job.fOutputStream          = ['tgtstops'          , 'ootstops'         , 'mothers'           ]
        job.fOutputDsID
        for i in range(0,3):
            job.fOutputDsID.append('bmup3'+s.name()+str(i+3)+'b0');
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

        job.fInputDataset          = Dataset('Mu2eII.bmup3s11b0.s2_sim.art','','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 10
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

        job.fInputDsID             = 'bmup3s41b0'
        dsn                        = project+'.'+job.fInputDsID+'.art'
        job.fInputDataset          = Dataset(dsn,'bmup3s41b0','local')
        job.fNInputFiles           = -1

        job.fMaxInputFilesPerSegment = 1
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '24h'
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

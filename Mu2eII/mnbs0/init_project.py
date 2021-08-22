#!/usr/bin/python

from local_classes import *
from mixing_inputs import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name] = Stage(name);
        return self.fStage[name]

    def __init__(self):

        project                          = 'Mu2eII'
        familyID                         = 'mnbs0'

        self.fProjectName                = project
        self.fDsid                       = familyID
        self.fStage                      = {}
#------------------------------------------------------------------------------
# s4:sim_b1 : 1M events = (1000 segments) x (1000 events/segment) ; limited by the output file size
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s4');
        job                          = s.new_job('sim_b1');

        job.fRunNumber               = 1000;
        job.fBaseFcl                 = project+'/'+familyID+'/'+s.name()+'_'+job.name()+'_'+familyID+'.fcl'
        job.fInputDataset            = Dataset('generator','mnbs0s00b1','local');
        job.fNInputFiles             = 1000                              # 1000
        job.fNEventsPerSegment       = 1000
        job.fMaxInputFilesPerSegment = 1

        define_mixing_inputs(job);

        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '5h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd

        job.fOutputPath              = [ 'defaultOutput'                ]
        job.fOutputDsID              = [ familyID+s.name()+'1b1'        ]
        job.fOutputFnPattern         = [ 'dig.mu2e.'+job.fOutputDsID[0] ]
        job.fOutputFormat            = [ 'art'                          ]

        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------

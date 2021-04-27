#
from local_classes import *
#------------------------------------------------------------------------------
# define common inputs for Mu2eII mixing
#------------------------------------------------------------------------------
def define_mixing_inputs(job):
    job.fAuxInputs               = {}

    # name, dataset, number_of_files per segment

    job.fAuxInputs['deuteronMixerTrkCal'] = ('physics.filters.deuteronMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.deut0s41b0.Mu2eII.art','deut0s41b0','local'),
                                             1);

    job.fAuxInputs['dioMixerTrkCal'     ] = ('physics.filters.dioMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.dio00s41b0.Mu2eII.art','dio00s41b0','local'),
                                             1);
    
    job.fAuxInputs['flashMixerTrkCal'   ] = ('physics.filters.flashMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.flsh1s51b0.Mu2eII.art','flsh1s51b0','local'),
                                             1);
    
    job.fAuxInputs['ootMixerTrkCal'     ] = ('physics.filters.ootMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.ootm0s41b0.Mu2eII.art','ootm0s41b0','local'),
                                             1);

    job.fAuxInputs['neutronMixerTrkCal' ] = ('physics.filters.neutronMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.neut0s41b0.Mu2eII.art','neut0s41b0','local'),
                                             1);

    job.fAuxInputs['photonMixerTrkCal'  ] = ('physics.filters.photonMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.epho0s41b0.Mu2eII.art','epho0s41b0','local'),
                                             1);

    job.fAuxInputs['protonMixerTrkCal'  ] = ('physics.filters.protonMixerTrkCal.fileNames', 
                                             Dataset('sim.mu2e.prot0s41b0.Mu2eII.art','prot0s41b0','local'),
                                             1);

#------------------------------------------------------------------------------
# inputs for mixing in CRV
#------------------------------------------------------------------------------
def define_mixing_inputs_crv(job):
    job.fAuxInputs               = {}

    # name, dataset, number_of_files per segment

    job.fAuxInputs['dioMixerCRV'        ] = ('physics.filters.dioMixerCRV.fileNames', 
                                             Dataset('sim.mu2e.dio00s42b0.Mu2eII.art','dio00s42b0','local'),
                                             1);
    
    job.fAuxInputs['neutronMixerCRV'    ] = ('physics.filters.neutronMixerCRV.fileNames', 
                                             Dataset('sim.mu2e.neut0s42b0.Mu2eII.art','neut0s42b0','local'),
                                             1);

    job.fAuxInputs['DSMixerCRV'         ] = ('physics.filters.DSMixerCRV.fileNames', 
                                             Dataset('sim.mu2e.crv01s32b0.Mu2eII.art','crv01s32b0','local'),
                                             1);
    
    job.fAuxInputs['PSMixerCRV'         ] = ('physics.filters.PSMixerCRV.fileNames', 
                                             Dataset('sim.mu2e.crv01s33b0.Mu2eII.art','crv01s33b0','local'),
                                             1);


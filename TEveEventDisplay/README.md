# Mu2e TEve Event Display Development Branch

This branch contains the prototype for the TEve based Mu2e Event Display.

Here are the details.

## The Module

TEveEventDisplay/src/TEveEventDisplay_module.cc is the Analyzer mdoule which currently controls the TEveEventDisplay. It is here where the navigation panel is drawn, events are accessed and plotting code is called. This is your "main" function.

## The fcl file

The prolog.fcl file resides in TEveEventDisplay/fcl and contains whats known as module instances for the TEveEventDisplay. Currently these are for : helix tracks, cosmic tracks or calo only events. We can add more.

## Running the code

You should access some of the Art examples in TEveEventDisplay/ArtExamples. You should access fcl files in TEveEventDisplay/CallerFcls

to run: ```$ mu2e -c PATH_TO_CALLER_FCL/File.fcl PATH_TO_ART/art.art --nevts 100 (for 100 events)```

The TEve Browser will appear. It is currently very slow.

## The TEve Event Display Infrastructure

This branch is very much a "work-in-progress" way it is currently structured is not optimal. Currently there are a number of Interfaces defined in ``dict classes`` in the ``src`` directory:

### gdml

The GDML file used here can be regenerated using: ```mu2e -c mu2eG4/fcl/gdmldump.fcl```. It contains the entire Mu2e World.

### Geom Interface

Contains callers for access to Tracker and Calo geometry. This class also contains functions to set visability of different elements absed on their names within the gdml.

### Drawing

Contains drawing options for both hits and clusters. These are cosmetic features.

### Collection Interface

Not yet used but will be as we move towards a more complex design. bsed on EventDisplay infrastructure.

### Particle Interface

The physics or PID world.

### NavState and EvtUtils

These contian the event navigation. It currently does not work as planned. This needs to be adapted to the ART world.



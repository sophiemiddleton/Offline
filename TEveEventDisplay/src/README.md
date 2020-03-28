This branch contains the prototype for the TEve based Mu2e Event Display.

Here are the details.

1) The Module

TEveEventDisplay/src/TEveEventDisplay_module.cc is the Analyzer mdoule which currently controls the TEveEventDisplay. It is here where the navigation panel is drawn, events are accessed and plotting code is called. This is your "main" function.

2) The fcl file

The prolog.fcl file resides in TEveEventDisplay/fcl and contains whats known as module instances for the TEveEventDisplay. Currently these are for : helix tracks, cosmic tracks or calo only events. We can add more.

3) Running the code

You should access some of the Art examples in TEveEventDisplay/ArtExamples. You should access fcl files in TEveEventDisplay/CallerFcls

to run: $ mu2e -c PATH_TO_CALLER_FCL/File.fcl PATH_TO_ART/art.art --nevts 100 (for 100 events)

The TEve Browser will appear. It is currently slow.



#### genomicsProcessing

Updated script runs from genomics.Rmd


Originally put together by Harshal Inamdar.

The script `dem.pl` will identify the sequening run to be processed, edit samplesheet 
to create new samplesheets which will then be used to demultiplex the run. The jobs would automatically be fired onto cluster.
The script is located on athens at `/ifs/data/Hiseq/Runs/` ; execute `perl dem.pl`  to start demultiplexing.


four.qsub,five.qsub,and six.qsub all are run on the SRA....404, SRA....405, SRA....406 respectively to get the barcodes.

freq.qsub was called multiple times to get the frequency of barcodes.

index.qsub and mapper.qsub make the index and map for Alevin.

The .R file makes the graphs.

The command_line.txt contains commands performed on the results of the qsub output:
The first command (tr) gets the counts for unique barcodes and writes it to another file (used in the freq.qsub).
The second block of awk gets the mean-standard deviation, and filters barcodes by counts.
The third two awk and grep commands extract files to map items for the Alevin Map variable.


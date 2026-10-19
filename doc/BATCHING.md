# Running in Batch Jobs
This page contains some tips and suggestions for how to process large projects in a batch queue environment with multiple processing machines available.

There are many different systems for batch queues and spreading work across multiple machines. In addition, there is little standardization of their APIs. So, rather than try to write code to interface with all of them, instead we provide information on the key areas in which our code can benefit from such systems, and leave it to the user to decide how that best fits their favourite approach.

## Concurrency
The main speedup available is by running the processing of individual tiles in some concurrent manner. Most of the scripts which process the tiled data have a simple loop over all tiles. However, in most cases, the processing for each tile in the loop is independant of the others. This means that the simplest way to benefit from concurrency is to split the set of tiles into groups and run the same tile loop for each group, in a separate batch job.

## Dividing The Tile List
Each of the scripts which would benefit from this approach takes a command line option `--groupMofN`. This will tell that invocation of the script that although it is processing tiles from the given input directory, it is expected that it will be only the M-th of N groups of tiles. This would mean N batch jobs, all running with the same `--indir` and `--outdir`, but each with a different value for `M`. For example, giving `--groupMofN 2 5` on the command line would mean that this invocation is for the 2nd group, dividing the whole set of tiles into 5 groups in total, one in each batch job. This would give approximately a 5-fold speedup, overall.

Using this does imply that all batch jobs have direct read/write access to the relevant directories, so there needs to be a file system shared between the machines.

It should be noted that the `alsv_tile_products` and `alsv_canopyheightmodel` scripts are independent of each other, and work from the same input files (indexed LAZ files with Height Above Ground added). If both of these are desired, they can be run concurrently as separate jobs, or sets of jobs.

## Mosaics
The mosaic products must, of course, process all tiles in one job. However, each product can be run independantly of the others. The `alsv_product_mosaic` command will, by default, loop over all products, creating each mosaic in turn. However, it also takes a `--stagecode` option, to limit to just a single product, thus allowing one batch job for each stage code, running concurrently.

The `alsv_flowacc` script runs on the complete DEM mosaic. It must wait until that has completed, and so cannot really be run concurrently with any other mosaic job.
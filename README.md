# IMPORTANT: THIS PIPELINE IS DEPRECATED

We recommend using the [VGP pipeline](https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish) based on FreeBayes instead of Pilon for polishing.

# PilonGrid

The distribution is a parallel wrapper around the [Pilon](https://github.com/broadinstitute/pilon) framework The pipeline is composed of bash scripts, an example mapping.fofn which shows how to input your fastq files (you give paths to the R1 file), and how to launch the pipeline. 

The current pipeline has been designed to run on the SGE scheduling system and has hard-coded grid resource request parameters. You must edit pilon.sh to match your grid options. It is, in principle, possible to run on other grid engines but will require editing all shell scripts to not use SGE_TASK_ID but the appropriate variable for your grid environment and editing the qsub commands in pilon.sh to the appropriate commands for your grid environment.

To run the pipeline you need to:

1. You must have a working java version (1.8+ recommended). Run ```javac *.java``` in the working folder
2. You must have a working installation of Pilon and symlink the jar file from wherever you have placed PilonGrid. You must also have a working bwa mem and samtools in your path.

3. Create a mapping.fofn input file specifying which fastq files you want to use as input. They can be gzipped/bzipped. The pipeline assumes that paired-end reads are used with one file named \*\_R1\_\* and the other name \*\_R2\_\*. It will fail if this is not the case.

4. Run the pipeline specifying the input fasta assembly.


```
sh pilon.sh asm.fasta
```

The pipeline is very rough and has undergone limited testing so user beware.

### CITE
If you find this pipeline useful, please cite the original Pilon paper:<br>
Walker BJ et. al. [Pilon: an integrated tool for comprehensive microbial variant detection and genome assembly improvement.](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963) PLoS One, 2014.

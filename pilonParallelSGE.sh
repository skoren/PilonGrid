#!/usr/bin/env bash

######################################################################
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################

ASM=`cat asm`
PREFIX=`cat prefix`
LEN=`cat lens`
SCRIPT_PATH=`cat scripts`

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
jobid=$1
fi

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

#make sure we have something to process
NUM_SAM=`ls $PREFIX*.sorted.bam 2>/dev/null |wc -l |awk '{print $1}'`
if [ $NUM_SAM -le 0 ]; then
   echo "Error: couldn't find any files named $PREFIX*.sorted.bam"
   echo "Please align some libraries to your referenc using bowtie2 or bwa and re-run pilon.sh"
   exit
fi

# now figure out which contig we are
NUM_JOBS=`wc -l $LEN |awk '{print $1}'`
if [ $jobid -le 0 ]; then
   echo "Invalid job id, must be 1 or greater"
   exit
fi

if [ $jobid -gt $NUM_JOBS ]; then
   echo "Invalid job id, max is $NUM_JOBS"
   exit
fi

line=`cat $LEN |head -n $jobid |tail -n 1`
utg=`echo $line |awk '{print $1}'`
len=`echo $line |awk '{print $NF}'`

if [ -e "$utg.pilon.vcf" ]; then
   echo "Already done"
else
   BAMLIST=""
   # first we need to subset bam
   for file in `ls $PREFIX*.sorted.bam`; do
      echo "Subsetting $file to be $utg $len"
      prefix=`echo $file |sed s/.sorted.bam//g`
      outputName="$jobid.$prefix.$jobid.sorted"
      samtools view -b -h $file $utg:0-$len | samtools sort -T $jobid.$prefix - -o  $outputName.bam
      samtools index $outputName.bam

      BAMLIST="$BAMLIST --frags $outputName.bam"
   done
   echo "$utg 0 $len" > $utg.$jobid.cut
   java -Xmx60g SubFasta $utg.$jobid.cut $ASM > $utg.$jobid.fasta

   java -Xmx60g -jar $SCRIPT_PATH/pilon.jar --fix bases --genome $utg.$jobid.fasta $BAMLIST --output $utg.pilon --changes --vcf --diploid

   for file in `ls *.$jobid.sorted.bam`; do
      echo "Removing file $file"
      rm $file
      rm $file.bai
   done
   rm $utg.$jobid.cut
   rm $utg.$jobid.fasta
fi

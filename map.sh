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

ASM=`cat asm |sed s/.fasta//g`

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
jobid=$1
fi

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

NUM=`cat mapping.fofn |wc -l`

if [ $jobid -gt $NUM ]; then
   echo "Invalid $jobid, max is $NUM"
   exit
fi

input=`cat mapping.fofn |head -n $jobid |tail -n 1`
file1=$input
file2=`echo $input |sed s/_R1_/_R2_/g`

if [ ! -e $file1 ]; then
   echo "Error file $file1 not found"
   exit
fi
if [ ! -e $file2 ]; then
   echo "Error file $file2 not found"
   exit
fi

prefix=`basename $file1`
prefix=`echo $prefix | awk -F "_R1_" '{print $1"_"$2}'|sed s/.fastq.gz//g`

echo "bwa mem -t 16 $ASM.fasta  $file1 $file2 2> $jobid.err  |samtools view -S -b - > $ASM.$prefix.unsorted.bam"
echo "samtools sort $ASM.$prefix.unsorted.bam $ASM.$prefix.sorted"
echo "samtools index $ASM.$prefix.sorted.bam"
bwa mem -t 16 $ASM.fasta  $file1 $file2 2> $jobid.err  |samtools view -S -b - > $ASM.$prefix.unsorted.bam
samtools sort $ASM.$prefix.unsorted.bam $ASM.$prefix.sorted
samtools index $ASM.$prefix.sorted.bam


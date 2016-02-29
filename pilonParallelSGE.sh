#!/usr/bin/env bash

######################################################################
#Copyright (C) 2015, Battelle National Biodefense Institute (BNBI);
#all rights reserved. Authored by: Sergey Koren
#
#This Software was prepared for the Department of Homeland Security
#(DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
#part of contract HSHQDC-07-C-00020 to manage and operate the National
#Biodefense Analysis and Countermeasures Center (NBACC), a Federally
#Funded Research and Development Center.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#* Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
#* Neither the name of the Battelle National Biodefense Institute nor
#  the names of its contributors may be used to endorse or promote
#  products derived from this software without specific prior written
#  permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      samtools view -b -h $file $utg:0-$len | samtools sort -  $outputName
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

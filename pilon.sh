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
echo "Usage: pilon.sh <assembly fasta>"
echo "This scrip will run pilon in parallel on the grid."
echo "It will index and sort alignments"

MACHINE=`uname`
PROC=`uname -p`
SCRIPT_PATH=$BASH_SOURCE
SCRIPT_PATH=`dirname $SCRIPT_PATH`
JAVA_PATH=$SCRIPT_PATH:.

ASM=$1
PREFIX=`echo $1 |sed s/.fasta//g |sed s/.fna//g |sed s/.fa//g`
LEN="$PREFIX.lens"

syst=`uname -s`
arch=`uname -m`
name=`uname -n`

if [ "$arch" = "x86_64" ] ; then
  arch="amd64"
fi
if [ "$arch" = "Power Macintosh" ] ; then
  arch="ppc"
fi

if [ x$ASM == x"" ]; then
   echo "Error: you must specify an assembly fasta file"
   exit
fi

if [ ! -e $ASM ]; then
   echo "Error: couldn't find $ASM, please try again"
   exit
fi
if [ ! -e $LEN ]; then
   java SizeFasta $ASM > $LEN
fi
NUM_CTG=`wc -l $LEN |awk '{print $1}'`
if [ $NUM_CTG -le 0 ]; then
   echo "Error: there are no contigs in the provided fasta file $ASM. Please try again"
   exit
fi
NUM_MAP=`wc -l mapping.fofn |awk '{print $1}'`
if [ $NUM_MAP -le 0 ]; then
   echo "Error: there are no input files to map provided in the mapping.fofn file. Please try again"
   exit
fi

echo "$ASM" > asm
echo "$PREFIX" > prefix
echo "$LEN" > lens
echo "$SCRIPT_PATH" > scripts

echo "Running with $PREFIX $ASM $NUM_MAP mappings on $NUM_CTG contigs"
qsub -V -q low.q -pe make-dedicated 1 -l mem_free=10g -cwd -N "${PREFIX}index" -j y -o `pwd`/index.out $SCRIPT_PATH/index.sh
qsub -V -q low.q -pe make-dedicated 16 -l mem_free=1g -t 1-$NUM_MAP -cwd -N ${PREFIX}map" -hold_jid "${PREFIX}index" -j y -o `pwd`/\$TASK_ID.map.out $SCRIPT_PATH/map.sh 
qsub -V -q low.q -pe make-dedicated 1 -tc 400  -l mem_free=60G -t 1-$NUM_CTG -cwd -N "${PREFIX}pilon" -hold_jid "${PREFIX}map" -j y -o `pwd`/\$TASK_ID.out $SCRIPT_PATH/pilonParallelSGE.sh

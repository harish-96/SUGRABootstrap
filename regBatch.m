(* ::Package:: *)

ClearAll[RunningFile]
RunningFile[fname_,nmax_,mreg_]:=Block[{},
"#!/bin/bash -l
#SBATCH --verbose
#SBATCH --partition=defq
#SBATCH --nodes 2
#SBATCH --job-name=N"<>ToString[nmax]<>"_M"<>ToString[Log10[mreg]]<>"
#SBATCH --output=/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/L10_2.out
#SBATCH --error=/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/L10_2.err
#SBATCH --mem 160192
#SBATCH --array=10-10
#SBATCH --time 24:00:00
MAT_NUM=$(( $SLURM_ARRAY_TASK_ID * 1 ))

module load sdpb/ubuntu-18.04/2.4.0

math -script /home/hmurali/string_bootstrap/SUGRABootstrap3/truncObj_reg2.m "<>ToString[fname]<>" $MAT_NUM "<>ToString[nmax]<>" 100 "<>ToString[mreg] ];

Nmin=20; Nmax=45;
Ns = Range[Nmin, Nmax];
Ms = Table[10^i, {i, 25, 25}];
(*Ns = {22, 25, 29, 33};*)
(*Ns = {10};*)

Do[Export["sdpb_run_"<>ToString[nn]<>"_"<>ToString[mreg]<>".txt",RunningFile[truncRegnp250ImTp,nn,mreg]],{nn,Ns},{mreg,Ms}];


Do[Run["sbatch sdpb_run_"<>ToString[nn]<>"_"<>ToString[mreg]<>".txt"],{nn,Ns},{mreg,Ms}];


Do[DeleteFile["sdpb_run_"<>ToString[nn]<>"_"<>ToString[mreg]<>".txt"],{nn,Ns},{mreg,Ms}];

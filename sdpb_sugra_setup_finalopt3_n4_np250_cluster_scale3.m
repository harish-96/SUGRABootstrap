(* ::Package:: *)

ClearAll[RunningFile]
RunningFile[fname_,nmax_]:=Block[{},
"#!/bin/bash -l
#SBATCH --verbose
#SBATCH --partition=amdq
#SBATCH --nodes 2
#SBATCH --job-name=L10_N"<>ToString[nmax]<>"
#SBATCH --output=/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/L10_2.out
#SBATCH --error=/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/L10_2.err
#SBATCH --mem 160192
#SBATCH --array=10-10
#SBATCH --time 24:00:00
MAT_NUM=$(( $SLURM_ARRAY_TASK_ID * 1 ))

module load sdpb/ubuntu-18.04/2.4.0

math -script /home/hmurali/string_bootstrap/SUGRABootstrap3/truncObj_reg.m "<>ToString[fname]<>" $MAT_NUM "<>ToString[nmax]<>" 100" ];

Nmin=36; Nmax=45;
Ns = Range[Nmin, Nmax];
(*Ns = {22, 25, 29, 33};*)
(*Ns = {10};*)

Do[Export["sdpb_run_"<>ToString[nn]<>"_"<>ToString[c]<>".txt",RunningFile[primalTruncObjReg1e6np250ImTp2,nn]],{nn,Ns},{c,0,0}];


Do[Run["sbatch sdpb_run_"<>ToString[nn]<>"_"<>ToString[c]<>".txt"],{nn,Ns},{c,0,0}];


Do[DeleteFile["sdpb_run_"<>ToString[nn]<>"_"<>ToString[c]<>".txt"],{nn,Ns},{c,0,0}];

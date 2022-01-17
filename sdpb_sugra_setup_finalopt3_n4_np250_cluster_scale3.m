(* ::Package:: *)

ClearAll[RunningFile]
RunningFile[fname_,nmin_,nmax_,lmax_]:=Block[{},
"#!/bin/bash -l
#SBATCH --verbose
#SBATCH --partition=defq
#SBATCH --nodes 3
#SBATCH --job-name=L"<>ToString[lmax]<>"
#SBATCH --output=/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/L"<>ToString[lmax]<>"_2.out
#SBATCH --error=/gpfs/hmurali/string_bootstrap/SUGRABootstrap3/L"<>ToString[lmax]<>"_2.err
#SBATCH --mem 160192
#SBATCH --array="<>ToString[nmin]<>"-"<>ToString[nmax]<>"
#SBATCH --time 24:00:00
MAT_NUM=$(( $SLURM_ARRAY_TASK_ID * 1 ))

module load sdpb/ubuntu-18.04/2.4.0

math -script /home/hmurali/string_bootstrap/SUGRABootstrap3/fullObj_cluster.m "<>ToString[fname]<>" "<>ToString[lmax]<>" $MAT_NUM 100" ];

Lmax=20;
Nmin=20; Nmax=35;
Ns = Range[Nmin, Nmax];
(*Ns = {22, 25, 29, 33};*)
(*Ns = {10};*)

Export["sdpb_run_"<>ToString[Nmin]<>"_"<>ToString[Nmax]<>".txt",RunningFile[primalFullObjnp250ImTp2,Nmin,Nmax,Lmax]];

Run["sbatch sdpb_run_"<>ToString[Nmin]<>"_"<>ToString[Nmax]<>".txt"]

DeleteFile["sdpb_run_"<>ToString[Nmin]<>"_"<>ToString[Nmax]<>".txt"]

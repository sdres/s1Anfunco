#!/bin/bash

sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_Nifti/ -f convertall -s 01 02 03 04 05 06 07 08 10 -ss 001 -c none --overwrite



sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_Nifties/ -f /base/heudiconv_heuristic.py -s 06 07 08 09 10 -ss 001 -c dcm2niix -b --overwrite


sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_Nifties/ -f /base/heudiconv_heuristic.py -s 02 -ss 002 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_Nifties/ -f /base/heudiconv_heuristic.py -s 11 12 -ss 001 -c dcm2niix -b --overwrite


sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_BIDS/ -f /base/heudiconv_heuristic.py -s 06 07 09 10 12 15 16 17 18 -ss 001 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_BIDS/ -f /base/heudiconv_heuristic.py -s 05 -ss 002 -c dcm2niix -b --overwrite

sudo docker run --rm -it -v /media/sebastian/Data/S1_anfunco/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_BIDS/ -f /base/heudiconv_heuristic.py -s 02 -ss 003 -c dcm2niix -b --overwrite


sudo docker run --rm -it -v /media/sebastian/Data/S1ANFUNCO/raw_data:/base nipy/heudiconv:latest -d /base/S1ANFUNCO_S{subject}/ses-{session}/*.IMA -o /base/S1ANFUNCO_BIDS/ -f /base/heudiconv_heuristic.py -s 02 -ss 001 -c dcm2niix -b --overwrite

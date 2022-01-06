#!/bin/bash

copyvectors=BrahuiYorubaSimulation.copyvectors.txt
idfile=BrahuiYorubaSimulation.idfile.txt
paramfile=BrahuiYorubaSimulation.SourcefindParamfile.txt

SOURCEFIND=../sourcefindv2.R

Rscript ${SOURCEFIND} \
	--chunklengths ${copyvectors} \
	--parameters ${paramfile} \
	--target BrahuiYorubaSimulation \
	--output BrahuiYorubaSimulation \
	--idfile BrahuiYorubaSimulation.idfile.txt

# HuaSmallTools
This repository contains some home-made bioinformic scripts. 
## Why I write this repository?
I am a student working on non-model plants. Unlike model organism, reams of annotation need to be done on non-model plants mannually. So
I put my home-made scripts in this repository.
## What's in this repository
For now, I assort the scripts into three functions:convert, parse and spider.
### parse scripts
These scripts mainly used for parse results for multiple Databases. eg.KEGG, GO, NCBI
1. parse_eggNOG.py
2. parse_KEGGkolist.py
This scripts was used for parse KOALA web server results to a more user-friendly format. And could be directly used for clusterProfile.
### convert scripts
Actually, just some simple id-mapping scripts
1. convert_SRX.py
### spider scripts
Though I named this module Spider, most of them just the wrapper of mainstream database API.

# HuaSmallTools
This repository contains some home-made bioinformic scripts. 
## Why I write this repository?
I am a student working on non-model plants. Unlike model organism, reams of annotation need to be done on non-model plants mannually. So
I put my home-made scripts in this repository.
## What's in this repository
For now, I assort the scripts into four functions:convert, parse, spider and some tools for sequence.
### convert
1. [convert_aspera.py](convert/convert_aspera.py)  
Convert SRP/SRX/SRS accession to SRR url which could be used for aspera directly.   
*Notes:* I sincerely recommend use the fastq mode. 
### parse
These scripts mainly used for parse results for multiple Databases. eg.KEGG, GO, NCBI
1. [parse_eggNOG.py](parse/parse_eggNOG.py)  
Parse KEGG and GO result from eggNOG 5.0 web server results to a table format.
2. [parse_KEGGkolist.py](parse/parse_KEGGkolist.py)  
Parse KOALA web server results with pathway information. The result file could be used in clusterProfiler.
3. [parse_IPSGO.py](parse/parse_IPSGO.py)   
Parse GO terms in InterPro Scan results.
4. [parse_BioProject.py](parse/parse_BioProject.py)  
Parse Biosample information from the directly downloaded BioProject txt file.
5. [parse_go_obofile.py](parse/parse_go_obofile.py)  
Parse GO database OBO file to a more user-friendly format. Some scripts in this repository may use its result.
### spider
1. [spider_chinese_name.py](spider/spider_chinese_name.py)  
Convert plant latin name to chinese name based on [植物智](http://www.iplant.cn/)  
*Notes:* I used the pyquery module in this scripts, which is not a default module in anaconda.
### sequence
These scripts was used for simple sequence processing (including MSA, GFF and other related files)
1. [calculate_psy_chemi.py](sequence/calculate_psy_chemi.py)  
[ExPASy - ProtParam tool](https://web.expasy.org/protparam/) is recommended for calculating physicochemical properties
of DNA sequence. However, it could not be applied to batch processing.
2. [trim_align.py](sequence/trim_align.py)  
A simple scirpt to trim alignment before construct the tree.
3. [extract_cds_prokka](sequence/extract_cds_prokka.py)  
Prokka results exclude the cds sequence file. This script could do this job.  
*Notes:* This script could only be used for prokaryote. And gff file from different software,in my experience
were different more or less. So, if you use it in places other than Prokka results, Please be careful. 
## Contact me
If you have any question in using these scripts or any advices, please feel free to contact me by njbxhzy@126.com.
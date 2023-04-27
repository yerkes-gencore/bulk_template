# Using the template #

Click the green button in the upper-right corner to make a new repo using this as the template, or use GitHub CLI to make a new repo from this template.

*The QC template relies on functions from the [gencore-bulk](https://github.com/yerkes-gencore/gencore-bulk) repo. See the link for install instructions*

Assuming the config is set up correctly (see next section), you can execute the script by running the `scripts/bulk_rnaseq_qc_template.runfile.Rmd`. This should load all data, set up the analysis objects, and render the initial QC report in the reports folder.  

 ### Experimental Design Table ###
 
The template relies on metadata provided in a table format. The table should include the following fields:

* **FileID:** sample-specific prefix to the STAR counts file. The redundant part of the file name (e.g. " _ReadsPerGene.out.tab") can be specified in the config as STARreadSuffix

* **SampleID** unique identifier for each sample, ideally short and easy to read as it is used for figure labels

Other fields can be supplied as relevant. You can refer to these columns by their name in the `QC_config.yml` to add layers to plots.

An example design table

|FileID  	| SampleID |  	Group |	Intgroup 	| Label |
| --- | --- | --- | --- | --- |
|p22074-s007_LiproxposFAC1_S147	| p22074-s007	| LiproxposFAC  |	LiproxposFAC |	LiproxposFAC_S147 |
|p22074-s008_LiproxposFAC2_S148 |	p22074-s008	| LiproxposFAC	| LiproxposFAC	| LiproxposFAC_S148 |
|p22074-s012_MTposFAC3_S152     |	p22074-s012	| MTposFAC      |	MTposFAC |	MTposFAC_S152 |
| … | … | … | … | … |
| #Required, unique and specific based on file names | #Required, unique, anything	| Not required	| Not required	| Not required |

### Config ###

A yaml-style config specifies variables to be used in the script. Here you specify the STAR output directory, STAR reference directory (used to rename genes from IDs to symbols), PCA mapping aesthetics, etc. 



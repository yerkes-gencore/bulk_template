Set up the sample metadata table in R for record of the transformation
Table should be called exp.design and have several columns specified

**FileID** should have the sample-specific prefix to the STAR counts file. The redundant part of the file name (e.g. " _ReadsPerGene.out.tab") can be specified in the config as STARreadSuffix

**SampleID** should be a unique identifier for each sample, Ideally short and easy to read as it is used for figure labels

**Group** should have the primary experimental condition(s) to be considered. May not apply in more complex studies

**Label** A short name for samples to be used in figures as labels for samples. E.g. S1 may not be informative enough for your intentions, but Day1_Sample1 might be better.

**Intgroup** A grouping of samples that will be evaluated together. In a more complex study with interaction terms, this might represent the interaction terms of interest. E.g., a time series study with 2 treatments might have values such as `D1_control` in this field. In a simpler study, this may be redundant with the `Group` term, but the template will still look for `Intgroup` when using the sample filter function, baseline normalization for heatmaps, etc. 

Additional metadata can be included in additional columns. Those columns can be named anything and can be referenced for figure mapping, etc. in the config

An example design table

|FileID  	| SampleID |  	Group |	Intgroup 	| Label |
| --- | --- | --- | --- | --- |
|p22074-s007_LiproxposFAC1_S147	| p22074-s007	| LiproxposFAC  |	LiproxposFAC |	LiproxposFAC_S147 |
|p22074-s008_LiproxposFAC2_S148 |	p22074-s008	| LiproxposFAC	| LiproxposFAC	| LiproxposFAC_S148 |
|p22074-s012_MTposFAC3_S152     |	p22074-s012	| MTposFAC      |	MTposFAC |	MTposFAC_S152 |
| … | … | … | … | … |
| #Required, unique and specific based on file names | #Required, unique, anything	| #Required, non unique	| #Required, non unique	| #Required, anything |
				
				

Some data may need to be specifically cast to ensure proper ggplot mapping. E.g. values in Batch may be categorical and should be factors, while continuous values (e.g. "Dosage") may want to be cast as numeric. 


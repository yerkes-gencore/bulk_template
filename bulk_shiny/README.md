# Shiny app
Open either the `ui.r` or `server.r` file in Rstudio, then click 'run app' in the middle right of rstudio (right next to a green play button).
Upload an analysis.RData object in the upload tab under analysis data. Yours is `p22136_Stacey_analysis.RData`
Give the app a few seconds to load the data
Check the QC tab and wait for the Read Counts figure to render

# Upload
Mostly for initialization, shouldn't have to come back after loading the Rdata object. You can upload genelists here with a bit more control over the input list delimiters. These lists can be used for heatmaps.

# DGE Analysis
Select your comparison groups. Log fold changes will be reported such that postive LFC indicate increased expression in the 'numerator' comparison group. 
The volcano plot is interactive. You can zoom it, hover points for more details, and lasso points for more information in the table below. See the tools wheel in the top right corner of the volcano plot to reset zoom, use the lasso, or export to svg (currently doesn't seem to be working).
Use the checkbox below the volcano plot to make the subsequent table show only lassoed genes, otherwise all gene data will be loaded.
You can apply filters in the table or search gene names. Click on a gene row (shift click for multiple) to highlight a gene and add it to your gene list using the button below the table.

# Heatmaps
Under `Genes of interest` tab
Under `heatmap` tab
Heatmaps can be generated with your genelist by checking the box. See the `My genelist` tab to curate your list and add genes manually. You can also import a genelist as a text file with one gene name per line to render a heatmap, just be sure to unselect the 'use my genelist' box. 

# GSEA
Select a database 
Choose your comparison groups
A positive enrichment score will correspond to more genes with positive LFC, or more expression in the numerator group.
The collection tab will let you look at individual genesets from your selected databases. You can see the GSEA enrichment plot and DGE data for the leading edge genes.
X-compare allows you to plot GSEA enrichment for a single geneset for multiple comparisons.
Heatmaps will generate a heatmap of relative log expression for all samples, focusing on genes from the selected geneset. Filtering the heatmap for N genes will filter by variance of RLE, not by a leading edge union, but should still be useful for distinguishing the most affected genes.

# Troubleshooting
If your app unexpectedly crashes, try to rerun the app from Rstudio by pressing run app again. If the green play button was replaced by 'reload app' It is possible that it crashed but R is still executing it. You can either try clicking in the console on the bottom of Rstudio and pressing `esc` or `cntrl+c`, or terminate R and restart (R studio might prompt you if R is not responding when you try to relaunch). 
If you have issues starting the app, you may have lost connection to Sblab. Make sure you have sblab mounted. If you had to reconnect to the server, you may have to close Rstudio and reopen the app.
If the app crashes while in full window mode, you might not be able to exit and will be stuck with a blank screen. Use `CMD+tab` to navigate outside the app's shell and try restarting using the above methods. 

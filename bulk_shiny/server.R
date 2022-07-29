options(shiny.maxRequestSize = 99*1024^2)
server <- function(input, output, session){
  ##### Load data
  #analysis <- list()
  #load('p22136_Stacey_Analysis.RData')
  analysis <- reactive({
    get(load(input$analysis_obj$datapath)) 
    #get(load(filename))
  })%>% bindCache(input$analysis_obj$datapath, cache='session') 
  sample_order <- reactive({fct_inorder(as.character(analysis()$ddsDrop@colData$SampleID))})
  pairs <- reactive({
    a <- utils::combn(unique(analysis()$sampleTable$Group),2)
    return(matrix(a, ncol=ncol(a)))
  })
  
  pairs.formatted <- reactive({apply(pairs(), 2, function(x){paste0(x[1],"_vs_", x[2])})})
  pairs.formatted.inverted <- reactive({apply(pairs(), 2, function(x){paste0(x[2],"_vs_", x[1])})})
  
  # eventReactive(
  #   analysis,
  #   {
  #     #load(file=input$analysis_obj$datapath, .GlobalEnv)
  #     
  #     
  #     #analysis()$rldDrop <- rlog(analysis()$ddsDrop, blind = FALSE, fitType = "parametric")
  #   }
  # 
  # )
  # load_data <- function(failed=FALSE){
  #   modalDialog(
  #     fileInput("analysis_obj", "Upload Analysis.RData object",
  #               multiple=FALSE, accept=c(".RData", ".rdata", ".Rdata")), 
  #     span('Upload a valid RData analysis object from the core'),
  #     if (failed){
  #       div(tags$b("Invalid object", style = "color: red;"))
  #     },
  #     footer=actionButton('load_analysis', "Load")
  #   )
  # }
  # showModal(load_data())
  # observeEvent(input$load_analysis,
  #              {
  #                showModal(modalDialog(print(input$analysis_obj)))
  #                if (is.null(input$analysis_obj)){
  #                  load_data(failed=TRUE)
  #                }
  #                e <- new.env()
  #                x <- load(input$analysis_obj$datapath, envir = e)
  #                analysis <- e[[x]]
  #                if (exists(analysis)){
  #                  removeModal()
  #                } else(
  #                  load_data(failed=TRUE)
  #                )
  #              })
  
  #showModal(modalDialog("Loading data", footer=NULL))
  #load(file="p22136_Stacey_Analysis.RData")
  #observeEvent(analysis()$config,
  #             {removeModal()})
  
  ##### import files
  
  # exp.design <- reactive({
  #   req(input$file2)
  #   df <- read.csv(input$file1$datapath,
  #                  header = input$header,
  #                  sep = input$sep,
  #                  quote = input$quote,
  #                  stringsAsFactors=FALSE)
  #   #df <- df %>% dplyr::arrange(across(all_of(analysis()$config$sampleGrouping)))
  #   df <- as.data.frame(sapply(df, as.factor))
  #   return(df)
  # })
  # output$exp.design <- renderDataTable({
  #   req(input$file1)
  #   
  #   if(input$disp == "head") {
  #     
  #     tmp <- head(exp.design())
  #   }
  #   else {
  #     tmp <- exp.design()
  #   }
  #   return(tmp)
  # })
  # config <- reactive({
  #   read_yaml(input$file2$datapath, fileEncoding = "UTF-8")
  # })
  
  #####
  ## sample mapping
  # output$mapping <- renderPlot({
  #   return(
  #     rownames_to_column(as_tibble(analysis()$mapBins, rownames = NA), var = "map_result") %>%
  #       pivot_longer(!map_result, names_to = "SampleID", values_to = "count") %>%
  #       #  mutate(SampleID = as.factor(str_sub(SampleID,end = -18))) %>%
  #       ggplot(aes(x = SampleID, y = count, fill = factor(map_result, levels = c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous","N_identified")))) +
  #       geom_bar(stat = "identity") +
  #       theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #       guides(fill=guide_legend(title="Map Result")) +
  #       ggtitle(paste0(analysis()$config$analysis," Mapping"))# + aes(fct_inorder(SampleID)) 
  #   )
  # }  )
  #####
  ## PCA
  #pca <- reactive({analysis()$pca})
  output$mapping <- renderPlot({analysis()$mapping_plot})
 
  output$RLE_raw <- renderPlot({analysis()$RLE_raw})
  output$RLE_norm <- renderPlot(analysis()$RLE_norm)
  output$PCA.plot <- renderPlot({analysis()$pca})
  
  #####
  ## sample comparisons
  output$SelectGroup1 <- renderUI({
    selectInput("comparison1", "Comparison numerator", unique(analysis()$sampleTable$Group))
  })
  output$SelectGroup2 <- renderUI({
    selectInput("comparison2", "Comparison denominator", unique(analysis()$sampleTable$Group))
  })
  #voltable <- comptable <- data.frame()
  pairwise.comparison <- eventReactive(
    input$pairwisecomp, {
      if (input$comparison1==input$comparison2){
        showModal(modalDialog("Select two different groups for comparison", easyClose=TRUE))
        return(NULL)
      } else {
        #req(input$comparison1!=input$comparison2)
        result <- as.data.frame(results(analysis()$ddsDrop, contrast=c("Group", input$comparison1, input$comparison2)))#, filter=maxMinFilter(analysis()$ddsDrop, "Group", c(numer_contrast,denom_contrast)))
        return(result)
      }
    })
  observe(output$summary_out <- renderPrint({req(is.data.frame(pairwise.comparison()))
    summary(pairwise.comparison())[,c(1,2,3,5,6)]}))
  volData <- eventReactive(
    pairwise.comparison(),{
      req(is.data.frame(pairwise.comparison()))
      result <- pairwise.comparison()
      volData <- result[!is.na(result$padj),]
      volData <- volData[order(volData$padj),]
      return(volData)
    })
  
  output$volcano_plot <- renderScatterD3({
    df <- volData() %>% mutate(threshold = padj < analysis()$config$alpha) %>% mutate(log10padj = -log10(padj+.Machine$double.xmin))
    #df_final <- sample_n(volData, 100)
    # display subsample of points that did not meet threshold (reduces lag)
    df_false_subsample <- sample_n(subset(df,threshold==FALSE),1000)
    # display all points that did meet threshold
    df_true_all <- subset(df,threshold==TRUE)
    # combine two df's above
    df_final <- rbind(df_true_all,df_false_subsample)
    # create tooltip for display
    tooltips <- paste("<b>Name: </b>", rownames(df_final),"<br /> <b>logFC: </b> ", round(df_final$log2FoldChange,3), "<br /> <b>padj: </b>",round(df_final$padj,3))#,"<br /> <b>Product: </b>",df_final$product)
    #urls <- paste0("http://ncbi.nlm.nih.gov/protein/",substr(df_final$Parent,6,nchar(df_final$Parent[1])))
    scatterD3(x=df_final$log2FoldChange,
              y=df_final$log10padj,
              xlab="logFC",
              ylab="-log10(padj)",
              axes_font_size="160%",
              col_var=df_final$threshold,
              col_lab="Significant?",
              #colors = c("TRUE" = "#C02942", "FALSE"="#53777A"),
              lines=data.frame(slope=c(Inf, Inf),
                               intercept=c(-2,2),
                               stroke="#000",
                               stroke_width=1,
                               stroke_dasharray=5),
              # #symbol_var=df_final$isHypothetical,
              # #symbol_lab="Hypothetical protein?",
              ylim=c(0,max(df_final$log10padj, na.rm=TRUE)+2),
              point_size=50,
              point_opacity=0.7,#,
              tooltip_text=tooltips,
              click_callback = "function(id, d) {
                      if(id && typeof(Shiny) != 'undefined') {
                          Shiny.onInputChange('selected_point', d.key_var);
                      }
                    }",
              lasso=TRUE,
              lasso_callback = "function(sel) {
                      Shiny.onInputChange('lassod_points', sel.data().map(function(d) {return d.key_var}));
                    }")
  })
  
  comptable <- eventReactive(
    volData(),{
      volData() %>% select(baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
        mutate_if(is.numeric,round,3) %>% rownames_to_column("Gene")
    })
  # observeEvent(input$reset_comptable,{
  #   input$lassod_points <- NULL
  # })
  # observeEvent(input$view_lassod_genes,
  #              {renderDataTable({comptable()[input$lassod_points,]}, rownames=FALSE)})
  output$comptable.rendered <- renderDataTable({
    if (input$view_lassod_genes){
      if (!is.null(input$lassod_points)){
        comptable()[input$lassod_points,]
      } else{
        comptable()
      }
    } else {
      comptable()
    }
    
  }, rownames=FALSE, filter="top")
  observe({
    req(is.data.frame(pairwise.comparison()))
    output$add.selected.button <- renderUI({actionButton("add.selected.rows", "Add selected genes to genelist")})
    #output$add.lassod.button <- renderUI({actionButton("add.lassod.points", "Add selected genes to genelist")})
  })
  # voltable <<- df_final
  # comptable <<- volData
  # result <<- result
  # return(volData)
  
  # })
  #       
  #       #url_var=urls)
  #       
  #       
  #       
  #      # return(TRUE)
  #      #}
  # 
  # 
  # selected.gene.volplot <- reactive({input$selected_point})
  # observeEvent(input$selected_point,
  #              {output$click_selected <- renderTable({comptable()[selected.gene.volplot(),c(1,3,6)]})}
  # )
  # lassod.genes <- reactive({comptable()[input$lassod_points,1]})
  # observeEvent(input$lassod_points,
  #              {output$lassod.points <- DT::renderDataTable(comptable()[input$lassod_points,c(1,3,6)], rownames=FALSE)}
  # )
  selected.genes <- reactive({comptable()[input$comptable.rendered_rows_selected,1]})
  output$selected.genes <- renderPrint({cat(selected.genes())})
  
  #####
  ## gene boxplots
  gene.comparison <- eventReactive(
    input$renderGenePlot, {
      validate(need(input$genename %in% rownames(analysis()$ddsDrop), "Gene name not in data"))
      validate(need(mean(plotCountsGKT(analysis()$ddsDrop, input$genename, intgroup = colnames(analysis()$ddsDrop@colData), returnData =TRUE)$count, na.rm=TRUE)>0, "Gene counts are all 0"))
      plotCountsGKT(analysis()$ddsDrop, input$genename, intgroup = colnames(analysis()$ddsDrop@colData), returnData =TRUE) %>% ggplot(aes(x=Group, y=count)) + geom_violin() + geom_boxplot(width=0.1, outlier.shape = NA) + xlab("Group") + geom_sina() + geom_text_repel(aes(label=Label)) + theme_bw() + labs(title = input$genename) + scale_y_log10()
      
    }
  )
  output$gene.comparison <- renderPlot({gene.comparison() })
  
  #####
  ## heatmap
  ##REMOVE?
  heatmap.input <- eventReactive(
    input$heatmap.render,
    {
      genes <- read.table(input$heatmap.genes$datapath)
      height <- unit(100,"mm") + nrow(genes)*unit(5, "mm")
      return(list(genes=genes, height=height))
    }
  )
  heatmap.plot <- eventReactive(
    input$heatmap.render,
    {
      if (input$use.genelist){
        genes.oi <- geneList()
      } else {
        genes.oi <- read.table(input$heatmap.genes$datapath)$V1
      }
      heatmap.height <- 350 + length(genes.oi)*13
      #analysis()$assayRlogForGSEA <- assay(analysis()$rldDrop)
      # filter low/no-expression genes
      #analysis()$assayRlogForGSEA <- analysis()$assayRlogForGSEA[rowMeans(analysis()$assayRlogForGSEA)>0,]
      #write_tsv(data.frame(Name = str_remove(rownames(analysis()$assayRlogForGSEA), "[.].*"), Description = "na", analysis()$assayRlogForGSEA), path = paste0("rlog_forGSEA_", analysis()$config$analysis, ".txt"))
      #analysis()$clsLinesGroup <- c(paste0(c(length(analysis()$rldDrop$Group),length(unique(analysis()$rldDrop$Group)),1), collapse = " "), paste0(c("#",unique(as.vector(analysis()$rldDrop$Group))), collapse = " "), paste0(analysis()$rldDrop$Group, collapse = " "))
      #write_lines(analysis()$clsLinesGroup, path = paste0("Group_",analysis()$config$analysis,".cls"))
      # read in the DESeq2 regularized log expression for Seq1
      hmap <- analysis()$rldDrop[rownames(analysis()$rldDrop) %in% genes.oi,]
      validate(need(nrow(hmap)>0, "No data for the provided genes"))
      baseline_means <- rowMeans2(assay(hmap[,hmap@colData$Group %in% analysis()$config$designBaseline]))
      hmap <- assay(hmap) - baseline_means
      plot <- Heatmap(hmap, show_row_names = TRUE, heatmap_legend_param = list(title="Relative log\nexpression"),
                      width = ncol(hmap)*unit(5, "mm"), 
                      height = nrow(hmap)*unit(5, "mm"),
                      column_order = as.character(sample_order()),
                      border='black',
                      #column_labels=rep("", ncol(hmap)),
                      top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(col=NA),#fill = 2:5),
                                                                          labels = unique(analysis()$ddsDrop@colData$Group), 
                                                                          labels_gp = gpar(col = "black", fontsize = 10), labels_rot=60, height=unit(2,"cm"))),
                      bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(col=NA),#fill = 2:5),
                                                                             labels = unique(analysis()$ddsDrop@colData$Group), 
                                                                             labels_gp = gpar(col = "black", fontsize = 10), labels_rot=60, height=unit(2,"cm"))),
                      # bottom_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=1:length(levels(as.factor(analysis()$ddsDrop@colData$Group)))),
                      #                                                      #labels=levels(as.factor(analysis()$ddsDrop@colData$Group)),
                      #                                                      labels_gp = gpar(col='white', fontsize=10))), 
                      column_split=analysis()$ddsDrop@colData$Group,
                      column_title = " ",
                      show_column_names = FALSE,
                      #column_title_rot = 60,
                      col=colorRamp2(c(min(hmap), 0, max(hmap)), c("blue", "white", "red")))
      
      
      #, clustering_distance_rows="spearman")#, height=heatmap.height)#, column_order=as.character(sample_order()))
      return(list(height=heatmap.height, plot=plot))
    }
  )
  observe(output$heatmap.plot <- renderPlot({heatmap.plot()$plot}, height={height=heatmap.plot()$height}))#[,!(analysis()$rldDrop@colData$Group %in% analysis()$config$designBaseline)])
  #shinyFileSave(input, "save.heatmap", roots=c(wd='.'))
  output$download.heatmap <- downloadHandler(filename='heatmap.png',
                                             content = function(file){
                                               png(file)
                                               heatmap.plot()$plot
                                               dev.off()
                                             }, contentType = "image/png")
  
  get_leading_ratio <- function(x){
    a <- unlist(str_split(x, ","))[3]
    b <- unlist(str_split(a, "="))[2]
    c <- as.numeric(sub("%", "", b))
    c/100
  }
  ##### 
  ## GSEA
  msigdb_labels <- read.table('msigdb_symbols.txt', header=TRUE)
  gmt.file <- list()
  for (a in msigdb_labels[,2]){
    gmt.file <- append(gmt.file, gmtPathways(file.path("./GSEA_collections/", a)))
  }
  #updateSelectizeInput(session, 'gsea.heatmap.collection', choices = names(gmt.file), server = TRUE)
  output$gsea.heatmap.collection.UI <- renderUI({selectizeInput("gsea.heatmap.collection", "Select collection",choices=names(gmt.file))})
  
  
  my_GSEA_dotplot <- function(res){
    
    
    res$leading_ratio <- round(unlist(lapply(res$leadingEdge, length)) / res$size, 2)
    ggplot(res) + geom_point(aes(x=leading_ratio, y=pathway, size=-log10(padj), color=NES)) + labs(y="", x="Fraction of genes in leading edge") + scale_color_gradient2() +
      theme_classic(11) +
      theme(panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_line(colour = "grey92"),
            panel.grid.major.y = element_line(colour = "grey92"),#element_blank(),
            panel.grid.minor.y = element_line(colour = "grey92")) +  
      scale_radius(name="NOM p-val", range=c(1,8), breaks=-log10(c(0.5,0.1,0.01,0.001)), limits=c(0,3), labels=c(0.5,0.1,0.01,0.001)) +
      scale_color_gradient2(low="blue", mid="white",high="red", breaks=c(-2,-1,0,1,2))#, limits=c(-2,2))
  }
  output$gsea.group1.UI <- renderUI({
    selectInput("gsea.group1", "GSEA comparison numerator", unique(analysis()$sampleTable$Group))
  })
  output$gsea.group2.UI <- renderUI({
    selectInput("gsea.group2", "GSEA comparison denominator", unique(analysis()$sampleTable$Group))
  })
  output$gsea.groups <- renderText({paste0('Comparing ', input$gsea.group1, " over ", input$gsea.group2)})
  #output$gsea.inputs.UI <- renderUI({selectizeInput("gsea.inputs", "Select comparisons",pairs.formatted())})
  output$gsea.collections.UI <- renderUI({selectizeInput("gsea.database", "Select Database",
                                                         choices=c(msigdb_labels$Symbol), multiple=TRUE)})
  gsea_input_data <- eventReactive(
    input$gsea.render,
    {
      #req(input$gsea.group1)
      #req(input$gsea.group2)
      validate(need(!is.null(input$gsea.group1) && !is.null(input$gsea.group2), "Select a comparison for GSEA"))
      validate(need(input$gsea.group1!=input$gsea.group2, "Select different gsea comparison groups"))
      #comp1 <- unlist(strsplit(input$gsea.inputs, "_"))
      comp2 <- input$gsea.group1#comp1[1]
      comp3 <- input$gsea.group2#comp1[length(comp1)]
      a <- as.data.frame(results(analysis()$ddsDrop, contrast=c("Group", comp2, comp3)))
      #a <- results(analysis()$ddsDrop, name=input$gsea.inputs)
      test_results <- as.data.frame(a) %>% na.omit() %>% mutate(ranking_stat = sign(log2FoldChange)*-log10(pvalue+.Machine$double.xmin)) %>% rownames_to_column("Gene")
      return(test_results)  
    }
  )
  
  observeEvent(
    input$gsea.render,
    {
      test_results <- gsea_input_data()
      ranks <- test_results$ranking_stat
      names(ranks) <- test_results$Gene
      ranks <- sort(ranks, decreasing=TRUE)
      ranks <<- ranks
      gmt.file <- list()
      
      
      for (db in input$gsea.database){
        gmt.file <- append(gmt.file, gmtPathways(file.path("GSEA_collections/", msigdb_labels[msigdb_labels$Symbol==db,2])))
      }
      # gmt.file <- dplyr::case_when(
      #   input$gsea.database=="H" ~ "mouse_hallmark.gmt"
      #   input$gsea.database=="C2" ~ "m5.all.v0.3.symbols.gmt"
      # )
      # gmt.file <- gmtPathways("mouse_hallmark.gmt")
      #gmt.file <<- gmt.file
      fgseaRes <- fgseaSimple(pathways=gmt.file,
                              stats=ranks,
                              minSize=15,
                              maxSize=500,
                              nperm=10000)
      
      topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
      topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
      output$GSEA_textsum <- renderText({paste0(length(gmt.file), " genesets tested, ",nrow(fgseaRes %>% filter(padj<0.05)), 
                                                " had an adjusted pvalue below 0.05. ",nrow(fgseaRes %>% filter(pval<0.05)),
                                                " genesets had a nominal pvalue below 0.05.")})
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
      output$gsea.toppaths <- renderPlot({
        validate(need(nrow(fgseaRes)>=1, "GSEA yielded no results"))
        plotGseaTable(gmt.file[topPathways], ranks, fgseaRes, 
                      gseaParam=0.5, render = TRUE)})
      output$select.gsea.geneset <- renderUI({
        selectizeInput("selected.gsea.geneset", "Geneset", names(gmt.file))
      })
      pathways_oi <- topPathways
      output$GSEA_overview_dotplot <- renderPlot({
        my_GSEA_dotplot(fgseaRes %>% filter(pathway %in% pathways_oi))
      })
      fgseaRes <<- fgseaRes
      output$gsea.genes.button <- renderUI({actionButton("add.gsea.genes", "Add selected genes to genelist")})
    }
  )
  observeEvent(
    input$selected.gsea.geneset,
    {
      output$gsea.geneset.plot <- renderPlot({plotEnrichment(gmt.file[[input$selected.gsea.geneset]], ranks)})
      #output$tmp <- renderText({unname(unlist(fgseaRes[fgseaRes$pathway==input$selected.gsea.geneset,8]))})
      gsea.LE_table <<- gsea_input_data() %>% filter(Gene %in% unname(unlist(fgseaRes[fgseaRes$pathway==input$selected.gsea.geneset,8]))) %>% mutate_if(is.numeric,round,3)
      output$gsea.LE_table <- renderDataTable({
        gsea.LE_table }, rownames=FALSE)
      
    }
  )
  
  #output$selected.genes <- renderPrint({cat(selected.genes())})
  
  
  selected.gsea.genes <- reactive({gsea.LE_table[input$gsea.LE_table_rows_selected,1]})
  
  output$gsea.xcomp.inputs.UI <- renderUI({selectizeInput("gsea.xcomp.inputs", "Select comparisons",c(pairs.formatted(), pairs.formatted.inverted()), multiple=TRUE)})
  
  output$gsea.xcomp.collection.UI <- renderUI({selectInput("gsea.xcomp.collection", "Select collection",unique(fgseaRes$pathway))})
  
  my_plotEnrichment <- function(pathway, stats,
                                gseaParam=1,
                                ticksSize=0.2) {
    
    #rnk <- rank(-stats)
    #ord <- order(-stats)#order(rnk)
    
    statsAdj <- stats / max(abs(stats))
    #statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    #statsAdj <- statsAdj / max(abs(statsAdj))
    
    pathway1 <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    genes <- names(statsAdj)[pathway1]
    pathway1 <- sort(pathway1)
    
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway1,
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway1 - 1, pathway1))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x=c(xs), y=c(ys), genes = as.vector(rbind(names(tops),names(tops))), tops=as.vector(rbind(rep(FALSE, length(tops)), rep(TRUE,length(tops)))))
    return(toPlot)
  }
  
  ### This should probably already be in the analysis object 
  # analysis()$rldDrop <- rlog(analysis()$ddsDrop, blind = FALSE, fitType = "parametric")
  # analysis()$assayRlogForGSEA <- assay(analysis()$rldDrop)
  # # filter low/no-expression genes
  # analysis()$assayRlogForGSEA <- analysis()$assayRlogForGSEA[rowMeans(analysis()$assayRlogForGSEA)>0,]
  observeEvent(
    input$gsea.xcomp.render,
    {
      gsea_res <- list()
      for (comp in input$gsea.xcomp.inputs){
        comp1 <- unlist(strsplit(comp, "_"))
        comp2 <- comp1[1]
        comp3 <- comp1[length(comp1)]
        a <- as.data.frame(results(analysis()$ddsDrop, contrast=c("Group", comp2, comp3)))
        b <- as.data.frame(a) %>% rownames_to_column('NAME') %>% filter(NAME %in% rownames(analysis()$assayRlogForGSEA)) %>% mutate(ranking_stat=sign(log2FoldChange)*-log10(pvalue+.Machine$double.xmin)) %>% column_to_rownames('NAME')
        #b <- as.data.frame(a) %>% na.omit() %>% mutate(ranking_stat = sign(log2FoldChange)*-log10(pvalue+.Machine$double.xmin))
        ranks <- b$ranking_stat
        names(ranks) <- rownames(b)
        ranks <- sort(ranks, decreasing=TRUE)
        # d <- list()
        # 
        # msigdb_labels <- read.table('msigdb_symbols.txt', header=TRUE)
        # for (db in ""){
        #   d <- append(d, gmtPathways(msigdb_labels[msigdb_labels$Symbol==db,2]))
        # }
        pathway <- gmt.file[input$gsea.xcomp.collection]
        fgseaRes <- fgseaSimple(pathways=pathway,
                                stats=ranks,
                                minSize=0,
                                maxSize=5000,
                                nperm=10000)
        RES <- my_plotEnrichment(gmt.file[[input$gsea.xcomp.collection]], ranks)
        RES$LE <- RES$genes %in% unlist(fgseaRes$leadingEdge)
        RES$Source <- comp
        gsea_res[[comp]] <- RES
      }
      gsea_xcomp <- do.call(rbind, gsea_res)
      
      colnames(gsea_xcomp) <- c("RANK.IN.GENE.LIST", "RUNNING.ES", "GENE", "actual","CORE.ENRICHMENT",  "source")
      p <- ggplot(gsea_xcomp, aes(x=RANK.IN.GENE.LIST)) +
        theme_classic(11) +
        theme(panel.grid.major = element_line(colour = "grey92"),
              panel.grid.minor = element_line(colour = "grey92"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) +
        scale_x_continuous(expand=c(0.01,0.01)) +
        geom_line(aes(y=RUNNING.ES, color=source),size=1, alpha=0.8)+
        geom_point(aes(y=RUNNING.ES, color=source, shape=CORE.ENRICHMENT), size=3, alpha=0.8) +
        scale_shape_manual(values=c(1, 16)) +
        #scale_color_discrete(labels=c()) +
        theme(legend.position = c(.8, .8),# legend.title = element_blank(),
              legend.background = element_rect(fill = "transparent"))+ ylab("Running Enrichment Score") + xlab("") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x=element_blank(),
              plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm")) +
        labs(color="Dataset") + guides(shape='none')
      
      gsea_xcomp.single <- gsea_xcomp %>% filter(actual==TRUE)
      gsea_xcomp.single$height <- abs(gsea_xcomp.single$RUNNING.ES - (gsea_xcomp[gsea_xcomp$actual==FALSE,'RUNNING.ES']))
      gsea_xcomp.single$height <- gsea_xcomp.single$height/max(gsea_xcomp.single$height)
      gsea_xcomp.single$y <- (2*as.integer(as.factor(gsea_xcomp.single$source)))-1
      #non scaled heights
      gsea_xcomp.single$height <- 1#gsea_xcomp.single$height/max(gsea_xcomp.single$height)
      gsea_xcomp.single$y <- (1*as.integer(as.factor(gsea_xcomp.single$source)))-1
      
      
      
      
      p2 <- ggplot(gsea_xcomp.single, aes(x = RANK.IN.GENE.LIST)) +
        geom_linerange(aes(ymin=y, ymax=y+height, color=source)) +
        xlab('Rank in Ordered Dataset') + ylab(NULL) + theme_classic(11) +
        theme(legend.position = "none",
              plot.margin = margin(t=-.5, b=0,unit="cm"),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) + #,
        #    axis.line.x = element_blank()) +
        scale_x_continuous(expand=c(0.01,0.01)) +
        scale_y_continuous(expand=c(0,0)) +
        #scaled
        #geom_hline(yintercept=2*(unique(as.integer(as.factor(gsea_xcomp$source)))-1))
        #non scaled
        geom_hline(yintercept=1*(unique(as.integer(as.factor(gsea_xcomp$source)))-1))
      output$gsea.xcomp.plot <- renderPlot({plot_grid(plotlist=list(p,p2), ncol=1, align="v", rel_heights = c(1+.1*length(input$gsea.xcomp.inputs), .1*length(input$gsea.xcomp.inputs)))})
    })
  #     
  #     
  #     output$gsea.xcomp.plot <-
  #   }
  
  ## GSEA heatmap
  
  
  gsea.heatmap.plot <- eventReactive(
    input$gsea.heatmap.render,
    {
      genes.oi <- unlist(gmt.file[[input$gsea.heatmap.collection]])
      
      
      #analysis()$assayRlogForGSEA <- assay(analysis()$rldDrop)
      # filter low/no-expression genes
      #analysis()$assayRlogForGSEA <- analysis()$assayRlogForGSEA[rowMeans(analysis()$assayRlogForGSEA)>0,]
      #write_tsv(data.frame(Name = str_remove(rownames(analysis()$assayRlogForGSEA), "[.].*"), Description = "na", analysis()$assayRlogForGSEA), path = paste0("rlog_forGSEA_", analysis()$config$analysis, ".txt"))
      #analysis()$clsLinesGroup <- c(paste0(c(length(analysis()$rldDrop$Group),length(unique(analysis()$rldDrop$Group)),1), collapse = " "), paste0(c("#",unique(as.vector(analysis()$rldDrop$Group))), collapse = " "), paste0(analysis()$rldDrop$Group, collapse = " "))
      #write_lines(analysis()$clsLinesGroup, path = paste0("Group_",analysis()$config$analysis,".cls"))
      # read in the DESeq2 regularized log expression for Seq1
      hmap <- analysis()$rldDrop[rownames(analysis()$rldDrop) %in% genes.oi,]
      
      baseline_means <- rowMeans2(assay(hmap[,hmap@colData$Group %in% analysis()$config$designBaseline]))
      hmap <- assay(hmap) - baseline_means
      if (nrow(hmap)>input$gsea.heatmap.genes.N){
        hmap <- head(hmap[order(rowVars(hmap), decreasing = TRUE),],input$gsea.heatmap.genes.N)
      } 
      #   
      # }
      validate(need(nrow(hmap)>0, "No data for genes in this pathway"))
      heatmap.height <- 350 + nrow(hmap)*13
      plot <- Heatmap(hmap, show_row_names = TRUE, heatmap_legend_param = list(title="Relative log\nexpression"),
                      width = ncol(hmap)*unit(5, "mm"), 
                      height = nrow(hmap)*unit(5, "mm"),
                      column_order = as.character(sample_order()),
                      border='black',
                      #column_labels=rep("", ncol(hmap)),
                      top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(col=NA),#fill = 2:5),
                                                                          labels = unique(analysis()$ddsDrop@colData$Group), 
                                                                          labels_gp = gpar(col = "black", fontsize = 10), labels_rot=60, height=unit(2,"cm"))),
                      bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(col=NA),#fill = 2:5),
                                                                             labels = unique(analysis()$ddsDrop@colData$Group), 
                                                                             labels_gp = gpar(col = "black", fontsize = 10), labels_rot=60, height=unit(2,"cm"))),
                      # bottom_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=1:length(levels(as.factor(analysis()$ddsDrop@colData$Group)))),
                      #                                                      #labels=levels(as.factor(analysis()$ddsDrop@colData$Group)),
                      #                                                      labels_gp = gpar(col='white', fontsize=10))), 
                      column_split=analysis()$rldDrop@colData$Group,
                      column_title = " ",
                      show_column_names = FALSE,
                      #column_title_rot = 60,
                      col=colorRamp2(c(min(hmap), 0, max(hmap)), c("blue", "white", "red")))
      
      
      #, clustering_distance_rows="spearman")#, height=heatmap.height)#, column_order=as.character(sample_order()))
      return(list(height=heatmap.height, plot=plot))
    }
  )
  observe(output$gsea.heatmap.plot <- renderPlot({gsea.heatmap.plot()$plot}, height={height=gsea.heatmap.plot()$height}))
  #shinyFileSave(input, "save.heatmap", roots=c(wd='.'))
  output$gsea.download.heatmap <- downloadHandler(filename='gsea.heatmap.png',
                                                  content = function(file){
                                                    png(file)
                                                    heatmap.plot()$plot
                                                    dev.off()
                                                  }, contentType = "image/png")
  
  
  
  ##### 
  ## Update genelist
  # output$manual_genelist_addition_button <- renderUI({
  #   selectizeInput("manual_genelist_addition", "Add a gene",rownames(analysis()$ddsDrop))
  # })
  searched_genes <- eventReactive(input$search_genelist,
                                  {
                                    data.frame(genes=rownames(analysis()$ddsDrop), a='a') %>% filter(grepl(input$manual_genelist_addition, genes, ignore.case = TRUE))
                                  })
  output$searched_genes_table <- renderDataTable({data.frame(gene=searched_genes()[,1])}, rownames=FALSE, options=list(searching=FALSE))
  selected.genes.genelist <- reactive({searched_genes()[input$searched_genes_table_rows_selected,1]})
  geneList <- reactiveVal(list())
  observeEvent(
    input$add_selected_searched_genes,
    {
      geneList(unique(append(geneList(), selected.genes.genelist())))
    }
  )
  observeEvent(
    input$add.lassod.points,
    {
      geneList(unique(append(geneList(), lassod.genes())))
      # geneList(unique(append(geneList(), selected.gene.volplot())))
      #output$genelist.table <- renderPrint({geneList()})
    }
  )
  observeEvent(
    input$add.selected.rows,
    {
      geneList(unique(append(geneList(), selected.genes())))
      #output$genelist.table <- renderPrint({geneList()})
      
    }
  )
  # observeEvent(geneList,
  #              {output$genelist.table <- DT::renderDataTable(as.data.table(geneList))}
  observeEvent(
    geneList(),
    {
      validate(need(length(geneList())>0, "Select genes in the analysis tab"))
      output$remove.selected.button <- renderUI({actionButton("remove.selected.rows", "Remove selected genes")})
      output$genelisttable <- DT::renderDataTable({data.frame(matrix(unlist(geneList()), byrow = TRUE), stringsAsFactors = FALSE)}, rownames=FALSE,colnames="")
    }
  )
  observeEvent(
    input$remove.selected.rows,
    {
      if (length(input$genelisttable_rows_selected)>0){
        geneList(geneList()[-c(input$genelisttable_rows_selected)])
      }
      
    }
  )
  observeEvent(
    input$add.gsea.genes,
    {
      geneList(unique(append(geneList(), selected.gsea.genes())))
      # geneList(unique(append(geneList(), selected.gene.volplot())))
      #output$genelist.table <- renderPrint({geneList()})
    }
  )
  ## import genelist
  observeEvent(
    input$load_genelist_file,
    {
      req(input$genelist_file)
      #if (input$genelist_sep == "\n"){
      df <- scan(input$genelist_file$datapath, what="character", sep=input$genelist_sep)
      #} else {
      #  df <- read.csv(input$genelist_file$datapath,
      #                 header = input$FALSE,
      #                 sep = input$genelist_sep,
      #                 stringsAsFactors=FALSE)
      #}
      
      #df <- df %>% dplyr::arrange(across(all_of(analysis()$config$sampleGrouping)))
      #df <- as.data.frame(sapply(df, as.factor))
      modalDialog(paste("Loaded", length(df), "genes", sep=" "))
      geneList(unique(append(geneList(), df)))
      
    })
  output$download.genelist <- downloadHandler(filename='genelist.txt',
                                              content = function(file){
                                                write.table(geneList(), file, quote=FALSE, row.names = FALSE, col.names=FALSE)
                                              })
  # output$raw.table <- DT::renderDT({
  #   analysis()$
  # })
  
  # analysis()$samplefileIDs = dir(config$alignmentDir)
  # analysis()$sampleSTARReads <- sapply(analysis()$samplefileIDs, function(sid) {read_tsv(paste0(
  #   dir(config$alignmentDir, pattern = sid, full.names = TRUE),
  #   "/",sid,
  #   config$STARreadSuffix
  # ), col_names = c("gene_id","unstranded_count","sense_count","antisense_count"))},
  # simplify = FALSE,
  # USE.NAMES = TRUE)
  # analysis()$sampleTable <- exp.design
  # output$config <- renderUI({
  #   req(input$file2)
  #   return(knitr::kable(t(as.data.frame(config())), format='html'))
  # })
  
}



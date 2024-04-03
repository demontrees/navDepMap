require(shiny)
require(bslib)
require(ComplexHeatmap)
require(matrixStats)
require(pathfindR)
require(shinyjs)

load(file = '/srv/shiny-server/navDepMap/data/cells_DM.RData')
load(file = '/srv/shiny-server/navDepMap/data/cellsCN_DM.RData')
load(file = '/srv/shiny-server/navDepMap/data/muts_DM.RData')
load(file = '/srv/shiny-server/navDepMap/data/celldatDM.RData')
load(file = '/srv/shiny-server/navDepMap/data/dnames_DM.RData')
load(file = '/srv/shiny-server/navDepMap/data/dr_DM.RData')
load(file = '/srv/shiny-server/navDepMap/data/kd_DM.RData')
load(file = '/srv/shiny-server/navDepMap/data/ctypes.RData')
load(file = '/srv/shiny-server/navDepMap/data/prot.RData')

dnames<-dnames[rownames(dr),]

#-----------UI
ui <- fluidPage(
  useShinyjs(),
  theme = bs_theme(bootswatch = 'darkly'),
  
  
  sidebarLayout(
    sidebarPanel( radioButtons(
                       'SEP', 'Separate By', c('RNA','Gene Copy Number','Protein','Mutation'),selected = 'Gene Copy Number'
                     ),
                  selectInput('CAN','Cancer Type',choices = unique(celldat$tcga_code),multiple = T),
                  selectInput('GEN','Genes',NULL,multiple = T),
                  tabsetPanel(
                       id = "params",
                       type = "hidden",
                       tabPanel("expression",
                                sliderInput('CUT','Cut',2,25,6),
                                checkboxGroupInput(
                                  'SEL', 'Select Clusters', seq(6)
                                )
                       ),
                       tabPanel("Mutation",
                                selectInput('MUT','Mutation type',choices = unique(muts$VariantInfo),multiple = T),
                       )),
                  textInput('DIR','Save Directory:',value = '~/navDepMap/results/'),
                  actionButton('SAV','Save Results')),
    mainPanel(navset_card_underline(
          title = "Graphs",
          nav_panel("Heatmap", plotOutput("HMP"),plotOutput("CTY")),
          nav_panel("Violin Plots", plotOutput("VIO"),plotOutput("VI2")),
          nav_panel("GSEA", plotOutput("GSE"))
        ))
                 
  
  
))
#--------------------------------Server

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

prot <- prot[which(colSums(is.na(prot))<200)]
cellsCN_ <- as.data.frame(cellsCN)
cells_ <- as.data.frame(cells)

server <- function(input, output, session) {

  observeEvent(input$SEP, {
    if (input$SEP %in% c('RNA',"Gene Copy Number",'Protein')){sep<-'expression'}
    else if (!input$SEP %in% c('RNA',"Gene Copy Number",'Protein')){sep<-'Mutation'}
    updateTabsetPanel(inputId = "params", selected = sep)
    updateTabsetPanel(inputId = "plots", selected = input$SEP)
  })
  sep<-reactive({input$SEP})
  cells_<-reactive({
    if (sep() == "RNA"){updateSelectizeInput(session, 'GEN',label = 'Gene', choices = colnames(cells), server = TRUE)
                        cells_.<-cells
    } else if (sep() == "Gene Copy Number"){updateSelectizeInput(session, 'GEN',label = 'Gene', choices = colnames(cellsCN), server = TRUE)
                        cells_.<-cellsCN
    } else if (sep() == "Mutation"){updateSelectizeInput(session, 'GEN',label = 'Protein', choices = unique(muts$HugoSymbol), server = TRUE)
                                    
    } else if (sep() == "Protein"){updateSelectizeInput(session, 'GEN',label = 'Protein', choices = colnames(prot), server = TRUE)
                        cells_.<-prot}
    cells_.
      })
  clls <- reactive({cl<-celldat$depMapID[which(celldat$tcga_code %in% input$CAN)]
                    cl<-cl[cl %in% rownames(cells_())]})
#-----------------------HEATMAP
  d <- reactive({dist(cells_()[clls(),input$GEN],  method = "euclidean", diag = FALSE)})
  hc <- reactive({hclust(d(), method="complete", members=NULL)}) # Clusters rows (i.e. CL).
  ct <- reactive({cutree(tree=hc(), k=input$CUT)})
  grp <- reactive({split(x = names(ct()), f = ct())})
  cols<-reactive({
    cls <- palette("Alphabet")[seq(input$CUT)]
    names(cls) <- unique(ct())
    cls
  })
  hmp<-reactive({if (input$SEP %in% c('RNA',"Gene Copy Number",'Protein')){Heatmap(t(cells_()[clls(),input$GEN]), top_annotation =HeatmapAnnotation(foo=ct(),col=list(foo=cols())), clustering_method_rows = "complete", cluster_rows = FALSE,show_column_names = FALSE, column_dend_height =  unit(25, "mm") , column_title = 'Cell Lines Clustered by Gene Count', show_row_names = TRUE,)} else{'No Heatmap for Mutation'}})
  output$HMP <- renderPlot(hmp())
  observeEvent(input$CUT, {
    updateCheckboxGroupInput(session,'SEL', 'Select Clusters', choices=seq(input$CUT))
  }) 
  grp_sel <- reactive({if (input$SEP %in% c('RNA',"Gene Copy Number",'Protein')){l=lapply(grp(), function(z){for (i in input$SEL){return(grp()[[i]])}})
                        unlist(l)} else {unique(muts$ModelID[which(muts$VariantInfo %in% input$MUT & muts$HugoSymbol %in% input$GEN & muts$ModelID %in% celldat$depMapID[which(celldat$tcga_code %in% input$CAN)])])}})
  grp_uns <- reactive({if (input$SEP %in% c('RNA',"Gene Copy Number",'Protein')){l=lapply(grp(), function(z){for (i in which(!seq(input$CUT) %in% input$SEL)){return(grp()[[i]])}})
                        unlist(l)} else {unique(muts$ModelID[which(!muts$VariantInfo %in% input$MUT & !muts$HugoSymbol %in% input$GEN & !muts$ModelID %in% grp_sel() & muts$ModelID %in% celldat$depMapID[which(celldat$tcga_code %in% input$CAN)])])}})
#------------------------------Ctype %
  cty<-reactive({
  grpsel<-grp_sel()
  grpsel<-grpsel[which(grpsel %in% unlist(ctypes))]
  grpuns<-grp_uns()
  grpuns<-grpuns[which(grpuns %in% unlist(ctypes))]
  m=list()
  ctyps<-ctypes[which(names(ctypes) %in% c(input$CAN))]
  for (s in list(grpsel,grpuns)){
    c=list()
    for (l in ctyps){m <- append(m, (length(which(s %in% l & !s%in%c))/length(which(s %in% unique(unlist(ctyps))))))
    c <- append(c, l)}}
  vals<-matrix(m, ncol = 2)
  colnames(vals)<-c('Selected','Unselected')
  rownames(vals)<-names(ctyps)
  
  par(mar=c(5, 4, 4, 8), xpd=TRUE)
  barplot(vals, col = palette('Alphabet')[1:length(ctyps)], main = 'Percent Cancer Type In\nHigh and Low Enriched Clusters' , cex.main = 1.2, beside = FALSE , cex.lab = 1,horiz = T)
  legend(1.05,2, legend=names(ctyps), title="C-Type",palette('Alphabet')[1:length(ctyps)],cex = 1)
  })
  output$CTY <- renderPlot(cty())
#-------------------------------DRUG VIO
    require(doMC)
    require(foreach)
    
  lisx_drugs = reactive({
    registerDoMC(cores = 48)
      grpsel <- grp_sel()
      grpsel <- grpsel[which(grpsel %in% colnames(dr))]
      grpuns <- grp_uns()
      grpuns <- grpuns[which(grpuns %in% colnames(dr))]
      rs <- (rowSums(is.na(dr[,grpsel]))<(dim(dr[,grpsel])[2]-2) | rowSums(is.na(dr[,grpuns]))<(dim(dr[,grpuns])[2]-2))
      dr.<-dr[which(rs),]
      lisx=foreach(tt = seq(dim(dr.)[1]), .inorder=T, .combine=c) %dopar%{
        aa = median(as.numeric(dr.[tt,grpuns]),na.rm = T)-median(as.numeric(dr.[tt,grpsel]),na.rm = T)
      }
      names(lisx)<-rownames(dr.)
      unregister_dopar()
      
      lisx
    })
    
  lisy_drugs = reactive({
    registerDoMC(cores = 48)
      grpsel <- grp_sel()
      grpsel <- grpsel[which(grpsel %in% colnames(dr))]
      grpuns <- grp_uns()
      grpuns <- grpuns[which(grpuns %in% colnames(dr))]
      rs <- (rowSums(is.na(dr[,grpsel]))<(dim(dr[,grpsel])[2]-2) | rowSums(is.na(dr[,grpuns]))<(dim(dr[,grpuns])[2]-2))
      dr.<-dr[which(rs),]
      lisy=foreach(tt = seq(dim(dr.)[1]), .inorder=T, .combine=c) %dopar%{
        aa = wilcox.test(x = as.numeric(dr.[tt,grpsel]), y = as.numeric(dr.[tt,grpuns]),
                  alternative = "two.sided")
        aa$p.value
      }
      names(lisy)<-rownames(dr.)
      unregister_dopar()
      
      lisy
    })
  eff_drugs<-reactive({l=lisx_drugs()[which(lisy_drugs()<0.05 & lisx_drugs() > 0)]
                        names(sort(l, decreasing = T))})
  res_drugs<-reactive({l=lisx_drugs()[which(lisy_drugs()<0.05 & lisx_drugs() < 0)]
                        names(sort(l, decreasing = F))})
  vio <- reactive({
    grpsel <- grp_sel()
    grpsel <- grpsel[which(grpsel %in% colnames(dr))]
    grpuns <- grp_uns()
    grpuns <- grpuns[which(grpuns %in% colnames(dr))]
    par(mar = c(18,4,4,4))
    vioplot(t(dr[eff_drugs()[1:20],grpsel]),col = 'red', plotCentre = "line", side = "left", na.rm=T, names=paste(as.vector(dnames[c(eff_drugs()[1:20]),c(1)]), as.vector(dnames[c(eff_drugs()[1:20]),c(2)]), sep ='\n' ),ylim=c(max(dr[eff_drugs()[1:20],c(grpsel,grpuns)],na.rm = T),(min(dr[eff_drugs()[1:20],c(grpsel,grpuns)], na.rm = T))), las = 3, cex.names = 0.5,cex.main = 1,cex = 0.7, cex.axis = 0.7, main = 'Significant Drug Performance Differences', ylab = 'Drug Response')
    vioplot(t(dr[eff_drugs()[1:20],grpuns]), col = 'blue', plotCentre = "line", side = "right", add = T, na.rm=T)
  })
  output$VIO <- renderPlot(vio())
#-----------------------------ESS VIO
  
  lisx_ess = reactive({
    registerDoMC(cores = 48)
    grpsel <- grp_sel()
    grpsel <- grpsel[which(grpsel %in% colnames(kd))]
    grpuns <- grp_uns()
    grpuns <- grpuns[which(grpuns %in% colnames(kd))]
    rs <- (rowSums(is.na(kd[,grpsel]))<(dim(kd[,grpsel])[2]-2) | rowSums(is.na(kd[,grpuns]))<(dim(kd[,grpuns])[2]-2))
    kd.<-kd[which(rs),]
    kd.<-kd.[which(rowVars(as.matrix(kd.[,grpsel]))>0 & rowVars(as.matrix(kd.[,grpuns]))>0),]
    lisx=foreach(tt = seq(dim(kd.)[1]), .inorder=T, .combine=c) %dopar%{
      aa = median(as.numeric(kd.[tt,grpuns]),na.rm = T)-median(as.numeric(kd.[tt,grpsel]),na.rm = T)
    }
    names(lisx)<-rownames(kd.)
    unregister_dopar()
    print(lisx)
    lisx
  })
  
  lisy_ess = reactive({
    registerDoMC(cores = 48)
    grpsel <- grp_sel()
    grpsel <- grpsel[which(grpsel %in% colnames(kd))]
    grpuns <- grp_uns()
    grpuns <- grpuns[which(grpuns %in% colnames(kd))]
    rs <- (rowSums(is.na(kd[,grpsel]))<(dim(kd[,grpsel])[2]-2) | rowSums(is.na(kd[,grpuns]))<(dim(kd[,grpuns])[2]-2))
    kd.<-kd[which(rs),]
    kd.<-kd.[which(rowVars(as.matrix(kd.[,grpsel]))>0 & rowVars(as.matrix(kd.[,grpuns]))>0),]
    lisy=foreach(tt = seq(dim(kd.)[1]), .inorder=T, .combine=c) %dopar%{
      aa = wilcox.test(x = as.numeric(kd.[tt,grpsel]), y = as.numeric(kd.[tt,grpuns]),
                  alternative = "two.sided")
      aa$p.value
    }
    names(lisy)<-rownames(kd.)
    unregister_dopar()
    print(lisy)
    lisy
  })
  eff_ess<-reactive({l=lisx_ess()[which(lisy_ess()<0.05 & lisx_ess() > 0 & !is.na(lisx_ess()) & !is.na(lisy_ess()))]
                      names(sort(l, decreasing = T))})
  res_ess<-reactive({l=lisx_ess()[which(lisy_ess()<0.05 & lisx_ess() < 0 & !is.na(lisx_ess()) & !is.na(lisy_ess()))]
                      names(sort(l, decreasing = F))})
  vi2 <- reactive({
    grpsel <- grp_sel()
    grpsel <- grpsel[which(grpsel %in% colnames(kd))]
    grpuns <- grp_uns()
    grpuns <- grpuns[which(grpuns %in% colnames(kd))]
    par(mar = c(6,4,4,4))
    vioplot(t(kd[eff_ess()[1:20],grpsel]),col = 'red', plotCentre = "line", side = "left", na.rm=T, names=eff_ess()[1:20],ylim=c(max(kd[eff_ess()[1:20],c(grpsel,grpuns)],na.rm = T),(min(kd[eff_ess()[1:20],c(grpsel,grpuns)], na.rm = T))), las = 3, cex.names = 0.5,cex.main = 1,cex = 0.7, cex.axis = 0.7, main = 'Significant Essentiality Differences', ylab = 'Essentiality')
    vioplot(t(kd[eff_ess()[1:20],grpuns]), col = 'blue', plotCentre = "line", side = "right", add = T, na.rm=T)
  })
  output$VI2 <- renderPlot(vi2())
#--------------------------------------GSEA
  gse<-reactive({
    both <-c(eff_ess()[1:200],res_ess()[1:200])
    both<-both[which(both %in% names(lisy_ess()))]
    p=as.numeric(unlist(lisy_ess()[both]))
    n=which(!is.na(p))
    both[n]
    p=p[n]
    fc=as.numeric(unlist(lisx_ess()[both]))
    fc=fc[n]
    df<-as.data.frame(list(both,fc,p), col.names = c('nm','fc','p'))
    run_pathfindR(df,gene_sets = "GO-All")
  })
  
  output$GSE <- renderPlot(gse())
  #---------------------Saving
  observeEvent(input$SAV,{
    
    shinyjs::html("SAV", "Saving...")
    
    both <- c(eff_drugs(),res_drugs())
    drugz<-unname(dnames[both,1])
    p=unname(unlist(lisy_drugs()[both]))
    fc=unname(unlist(lisx_drugs()[both]))
    df<-as.data.frame(list(drugz,fc,p), col.names = c('nm','fc','p'))
    write.csv(df,paste(input$DIR,'drugs.csv',sep=''))
    
    both <-c(eff_ess(),res_ess())
    p=unname(unlist(lisy_ess()[both]))
    fc=unname(unlist(lisx_ess()[both]))
    df<-as.data.frame(list(both,fc,p), col.names = c('nm','fc','p'))
    write.csv(df,paste(input$DIR,'genes.csv',sep=''))
    
    sepr <- paste('Analysis:\n',paste((input$SEP)),collapse = ' ')
    ctyp <- paste('\n\nCancer Types:\n',paste(c(input$CAN),collapse = ', '),collapse = '')
    gens <- paste('\n\nGenes:\n',paste(c(input$GEN),collapse = ', '),collapse = '')
    effd <- paste('\n\nTop Drugs:\n',paste(dnames$Drug.Name[which(rownames(dnames) %in% eff_drugs()[1:20])],collapse = ', '),collapse = '')
    resd <- paste('\n\nResistant Drugs:\n',paste(dnames$Drug.Name[which(rownames(dnames) %in% res_drugs()[1:20])],collapse = ', '),collapse = '')
    effg <- paste('\n\nTop Genes:\n',paste(eff_ess()[1:20],collapse = ', '),collapse = '')
    resg <- paste('\n\nResistant Genes:\n',paste(eff_ess()[1:20],collapse = ', '),collapse = '')
    write(paste(c(sepr,ctyp,gens,effd,resd,effg,resg),collapse = ''),file =paste(input$DIR,"results.txt", sep = ''))
    
    tiff(filename = paste(input$DIR,'VioDrugs.tif',sep=''),height = 10000,width = 10000, res = 1200, pointsize = 12, compression = "lzw+p")

    print(vio())
    dev.off()
    
    tiff(filename = paste(input$DIR,'VioGenes.tif',sep=''),height = 5000,width = 10000, res = 1200, pointsize = 12, compression = "lzw+p")

    print(vi2())
    dev.off()
    
    tiff(filename = paste(input$DIR,'Heatmap.tif',sep=''),height = 5000,width = 10000, res = 1200, pointsize = 12, compression = "lzw+p")
 
    print(hmp())
    dev.off()
    
    tiff(filename = paste(input$DIR,'GSEA.tif',sep=''),height = 10000,width = 10000, res = 1200, pointsize = 12, compression = "lzw+p")

    gse()
    dev.off()
    
 
    shinyjs::html("SAV", "Save Results")
    })
}
#----------END

shinyApp(ui, server)

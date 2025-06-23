## qiimera functions

fmtBlast <- function(bo_in = blasted) {
    
    blast_outs <- bo_in
    blast_outs <- split(blast_outs,seq(nrow(blast_outs)))
    blast_outs <- lapply(blast_outs, function(X){
    orig <- X
    n.len  <- unlist(strsplit(X[,7],"\n"))  

    fl <-  gsub("\n.*", "",X)
    subs <- paste0('^', paste0(fl[7],'\n'))
    X <- X[,7]
    X <- gsub(subs,'',X)


    if(length(n.len)==1){
        return(orig)
    } else {
        X <- unlist(strsplit(X,"\n"))
        X <- strsplit(X,",")
        #X <- lapply(X, strsplit,",")
        X <- rbind(fl,do.call(rbind,X))
        colnames(X) <- colnames(orig)
        return(data.frame(X, row.names=NULL))
    }
})
blast_outs <- data.frame(do.call(rbind,blast_outs), row.names=NULL)
return(blast_outs)
}


### gets the taxa with the smallest evalue for each query. returns multiple hits when there are ties
bhit_sum <- function(id, input ){
        indx <- which(input$qseqid == id)
        hits <- input[indx,][which( input[indx,'evalue'] == min(as.numeric(input[indx,'evalue']) ) ), ]
        
        species_sci <- unique(hits$sscinames)
        species_com <- unique(hits$scomnames)

        is.fish <- ifelse(species_sci %in% fish, 'YES', '')
        is.fish <- paste(is.fish, collapse = ';')

        species_sci <- paste(species_sci, collapse = ';')
        species_com <- paste(species_com, collapse = ';')

        species_king <- paste(unique(hits$sskingdoms), collapse = ';')
        species_com <- paste(unique(hits$scomnames), collapse = ';')
        
        evalue <- paste(unique(hits$evalue), collapse = ';')
        pident <- paste(unique(hits$pident), collapse = ';')

        return(data.frame(species_king,species_com,species_sci,evalue,pident,is.fish))

}

### gets the taxa with the smallest evalue for each query. returns multiple hits when there are ties
bhit_sums <- function(id, input ){
        indx <- which(input$qseqid == id)
        hits <- input[indx,][which( input[indx,'evalue'] == min(as.numeric(input[indx,'evalue']) ) ), ]
        
        species_sci <- unique(hits$sscinames)
        species_com <- unique(hits$scomnames)

        species_sci <- paste(species_sci, collapse = ';')
        species_com <- paste(species_com, collapse = ';')

        species_king <- paste(unique(hits$sskingdoms), collapse = ';')
        species_com <- paste(unique(hits$scomnames), collapse = ';')
        
        evalue <- paste(unique(hits$evalue), collapse = ';')
        pident <- paste(unique(hits$pident), collapse = ';')

        return(data.frame(species_king,species_com,species_sci,evalue,pident))

}

### count up reads for each OTU
collapseWallASVs <- function(tax_sub,feat_table){
  ## OTU_out
  OTU_out <- lapply(unique(tax_sub$Taxon), function(x) {
      #rownames(tax_sub)[tax_sub$Taxon == x] 
      asvs <- rownames(tax_sub)[tax_sub$Taxon == x] 
      colSums(data.frame(feat_table[asvs,,drop=F], drop = FALSE, check.names = FALSE))
      })
  OTU_out <- do.call(rbind, OTU_out)

  ## OTU_Tax_out
  OTU_Tax_out <- data.frame(Feature.ID = unique(tax_sub$Taxon) ,Taxon = unique(tax_sub$Taxon)) %>% parse_taxonomy()
  OTU_Tax_out$Common <- com[OTU_Tax_out$Species]
  rownames(OTU_Tax_out) <- NULL


  return(data.frame(OTU_Tax_out, OTU_out, check.names = FALSE))
}


add_tax <- function(Updated_ASV = all_ASVs, samps = all){
  ### Collapse ASVs to species OTUs and write to file 
  tax_paste <- apply(Updated_ASV[,colnames(update_ptax)[1:7]], 1, paste, collapse = '; ')

  ## Table of ASVs to OTUs
  OTU_out <- lapply(unique(tax_paste), function(x) {
    #rownames(tax_sub)[tax_sub$Taxon == x] 
    asvs <- names(tax_paste)[tax_paste == x]
    colSums(data.frame(Updated_ASV[asvs,samps,drop=F], check.names = FALSE))
  })
  OTU_out <- do.call(rbind, OTU_out)
  rownames(OTU_out) <- unique(tax_paste)

  OTU_Tax_out <- data.frame(Feature.ID = unique(tax_paste) ,Taxon = unique(tax_paste), row.names = unique(tax_paste) ) %>% parse_taxonomy()
  Updated_OTUs <- data.frame(OTU_Tax_out[unique(tax_paste),],Common = com[OTU_Tax_out[unique(tax_paste),'Species']],OTU_out[unique(tax_paste),], check.names = FALSE, row.names = NULL)
  Updated_OTUs <- Updated_OTUs[rowSums(Updated_OTUs[,samps]) > 0,]
  return(Updated_OTUs)
}

add_otu_taxa <- function(Updated_ASV = target_ASVs, tax = update_ptax, samps = colnames(feat_table), com = blast_names){
  ### Collapse ASVs to species OTUs and write to file 

  txlevels <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species','Common')

  tax_paste <- apply(tax[,txlevels], 1, paste, collapse = '; ')

  ## Table of ASVs to OTUs
  OTU_out <- lapply(unique(tax_paste), function(x) {
    #rownames(tax_sub)[tax_sub$Taxon == x] 
    asvs <- names(tax_paste)[tax_paste == x]
    colSums(data.frame(feat_table[asvs,samps,drop=F], check.names = FALSE))
  })
  OTU_out <- do.call(rbind, OTU_out)
  rownames(OTU_out) <- unique(tax_paste)



  OTU_Tax_out <- data.frame(Feature.ID = unique(tax_paste) ,Taxon = unique(tax_paste), row.names = unique(tax_paste) ) %>% parse_taxonomy()
  OTU_Tax_out$Common <- com[OTU_Tax_out$Species,'com']
  Updated_OTUs <- data.frame(OTU_Tax_out[unique(tax_paste),], OTU_out[unique(tax_paste),], check.names = FALSE, row.names = NULL)
  Updated_OTUs <- Updated_OTUs[rowSums(Updated_OTUs[,samps]) > 0,]
  return(Updated_OTUs)
}
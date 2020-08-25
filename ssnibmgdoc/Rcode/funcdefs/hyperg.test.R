# Define the function for hypergeometric test
require(Category)
hyperg <- Category:::.doHyperGInternal
hyperg.test <-
  function(pathway.genes, genes.of.interest, all.geneIDs, over=TRUE)
  {
    white.balls.drawn <- length(intersect(genes.of.interest, pathway.genes))
    white.balls.in.urn <- length(pathway.genes)
    total.balls.in.urn <- length(all.geneIDs)
    black.balls.in.urn <- total.balls.in.urn - white.balls.in.urn
    balls.pulled.from.urn <- length(genes.of.interest)
    hyperg(white.balls.in.urn, black.balls.in.urn,
         balls.pulled.from.urn, white.balls.drawn, over)
  }


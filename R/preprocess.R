
#' Pre-processing function for sex classification
#'
#' The purpose of this function is to process a single cell counts matrix into
#' the appropriate format for the \code{classifySex} function.
#'
#' This function will filter out cells that are unable to be classified due to
#' zero counts on *XIST/Xist* and all of the Y chromosome genes. If
#' \code{qc=TRUE} additional cells are removed as identified by the
#' \code{perCellQCMetrics} and \code{quickPerCellQC} functions from the
#' \code{scuttle} package. The resulting counts matrix is then log-normalised
#' and scaled.
#'
#' @param x the counts matrix, rows are genes and columns are cells. Row names
#' must be gene symbols.
#' @param genome the genome the data arises from. Current options are
#' human: genome = "Hs" or mouse: genome = "Mm".
#' @param qc logical, indicates whether to perform additional quality control
#' on the cells. qc = TRUE will predict cells that pass quality control only
#' and the filtered cells will not be classified. qc = FALSE will predict
#' every cell except the cells with zero counts on *XIST/Xist* and the sum
#' of the Y genes. Default is TRUE.
#'
#' @return outputs a list object with the following components
#' \item{tcm.final }{A transposed count matrix where rows are cells and columns
#' are the features used for classification.}
#' \item{data.df }{The normalised and scaled \code{tcm.final} matrix.}
#' \item{discarded.cells }{Character vector of cell IDs for the cells that are
#' discarded when \code{qc=TRUE}.}
#' \item{zero.cells }{Character vector of cell IDs for the cells that can not
#' be classified as male/female due to zero counts on *Xist* and all the
#' Y chromosome genes.}
#'
#' @importFrom AnnotationDbi select
#' @importFrom stringr str_to_title
#' @importFrom scuttle perCellQCMetrics
#' @importFrom scuttle quickPerCellQC
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export preprocess
#'
#' @examples
#'
#' library(speckle)
#' library(SingleCellExperiment)
#' library(CellBench)
#' library(org.Hs.eg.db)
#'
#' # Get data from CellBench library
#' sc_data <- load_sc_data()
#' sc_10x <- sc_data$sc_10x
#'
#' # Get counts matrix in correct format with gene symbol as rownames
#' # rather than ENSEMBL ID.
#' counts <- counts(sc_10x)
#' ann <- select(org.Hs.eg.db, keys=rownames(sc_10x),
#'              columns=c("ENSEMBL","SYMBOL"), keytype="ENSEMBL")
#' m <- match(rownames(counts), ann$ENSEMBL)
#' rownames(counts) <- ann$SYMBOL[m]
#'
#' # Preprocess data
#' pro.data <- preprocess(counts, genome="Hs", qc = TRUE)
#'
#' # Look at counts on XIST and superY.all
#' plot(pro.data$tcm.final$XIST, pro.data$tcm.final$superY)
#'
#' # Cells that are identified as low quality
#' pro.data$discarded.cells
#'
#' # Cells with zero counts on XIST and all Y genes
#' pro.data$zero.cells
#'
preprocess<- function(x, genome=genome, qc=qc){

    x <- as.matrix(x)
  
  if (length(unique(colnames(x))) != ncol(x)){
    message("Cell names are missing/duplicated. Cells are renamed to cell1 - cell", ncol(x))
    colnames(x) = paste(rep("cell", ncol(x)), seq(1, ncol(x)), sep="")
  }

  # genes located in the X chromosome that have been reported to escape
  # X-inactivation
  Xgenes<- c("Xist","Eif2s3x","Timm17b","Bcor","Kdm6a",
             "Lancl3","5530601H04Rik","Cybb","Pbdc1","Kdm5c",
             "Jpx","Ddx3x","Ftx","Gm14719","Firre",
             "Mospd1","Utp14a","Atp11c","Pnma3","Slc16a2",
             "Plp1","Klhl15","Gprasp1","Gdi1","Syp",
             "Gpm6b","Gla","Tmem47","Vsig4","Amot",
             "Cfp","Fgd1","Tasl","Fam120c","5730416F02Rik",
             "Phf8","Bgn","Rnf128","Klf8","Car5b",
             "Tmsb4x","Ezhip","Mid1","Bmp15","Enox",
             "Lamp2","Rp2h","AU015836","Gm14820","Kif4",
             "Rlim","Sh3bgrl","Fam199x","Tmem164","Alg13",
             "Tmem29","Pdha1","Flna")

  # genes belonging to the  chromosome Y (unique genes, mm39)
  Ygenes<-c("1700020D14Rik","Gm29089","Zfy1","Uba1y","Gm28588",
            "Gm28587","Kdm5d","Eif2s3y","Gm29650","Tspy-ps",
            "Uty","Ddx3y","Usp9y","Zfy2","H2al2c",
            "Sry","H2al2b","Gm4064","Rbmy","Gm10256",
            "Gm10352","Gm29289","Gm21677","Gm21693",
            "Gm21704","Gm21708","Gm3376","Gm21064","Gm28242",
            "ENSMUSG00000121460","Gm28357","Gm29351","Gm29349",
            "Gm20918","Gm29353","Gm29158","Gm21854","Gm20914",
            "Gm28171","Gm28173","Gm29194","Gm21778","Gm29198",
            "Gm29193","Gm28430","Gm28593","Gm20873","Gm28919",
            "Gm21719","Gm28444","Gm20830","Gm28442","Gm28445",
            "Gm28999","Gm21244","Gm28998","Gm28398","Gm28395",
            "Gm28394","Gm28393","Gm20826","Gm28576","Gm28575",
            "Gm28571","Gm28570","Gm29049","Gm20825","ENSMUSG00000121461",
            "Gm29046","Gm29043","Gm29044","Gm20824","Gm28147",
            "Gm28331","Gm20815","Gm29364","Gm29363","Gm21292",
            "Gm21721","Gm29329","Gm21812","Gm21874","Gm20821",
            "Gm21310","Gm28509","Gm20834","Gm28510","Gm20737",
            "Gm29527","Gm29522","Gm20777","Gm28955","Gm28954",
            "Gm20828","Gm28656","Gm29554","Gm20812","Gm37739",
            "Gm20807","Gm38084","Gm37231","Gm21425","Gm37252",
            "Gm21440","Gm21454","Gm37147","Gm38159","Gm37870",
            "Gm37378","Gm20773","Gm37561","Gm37263","Gm37130",
            "Gm37467","Gm38363","Gm30174","Gm30353","Gm30686",
            "Gm30737","Gm37654","Gm37572","Gm37344","Gm37059",
            "Gm31942","Gm37440","Gm37577","Gm20822","Gm21809",
            "Gm20877","Gm38003","Rbm31y","Gm37865","Gm36941",
            "Gm37937","Gm20772","Gm20831","Gm37286","Gm38028",
            "Ssty1","Gm38072","Gm37473","Gm37656","Gm38127",
            "Gm38361","Gm38209","Gm37734","Gm38013","Gm37948",
            "Gm34550","Gm37544","Gm34716","Gm37635","Gm35070",
            "Gm35134","Gm37998","Gm38168","Gm37740","Gm36261",
            "Gm37898","Gm36467","Gm37927","Gm36929","Gm38174",
            "Gm29866","Gm36950","Gm20909","Gm20865","Gm37462",
            "Gm38370","Gm38054","Gm37690","Gm31571","Gm32114",
            "Gm37075","Gm38371","Gm37657","Gm37721","Gm21366",
            "Gm33954","Gm37538","Gm34217","Gm37952","Gm37454",
            "Gm37840","Gm37574","Gm35670","Gm37687","Gm38296",
            "Gm36345","Gm37157","Gm20809","Gm37434","Gm30045",
            "Gm37112","Gm37071","Gm37547","Gm30705","Gm37451",
            "Gm37875","Gm31422","Gm37798","Gm37346","Gm38136",
            "Gm37808","Gm32181","Gm28853","Gm28852","Gm28858",
            "Gm20890","Gm28245","Gm28246","Gm28249","Gm28250",
            "Gm28244","Gm29021","Gm28880","Gm21488","Gm28454",
            "Gm28878","Gm29399","Gm28799","Gm28798","Gm20894",
            "Gm28938","Gm29122","Gm28985","Gm28725","Gm28984",
            "Gm28987","Gm28986","Gm28129","Gm28127","Gm28128",
            "Gm28132","Gm28130","Gm28131","Gm21679","Gm29316",
            "Gm29315","Gm28216","Gm28217","Gm28092","Gm28091",
            "Gm29580","Gm28089","Gm28088","Gm29579","Gm29582",
            "Gm21529","Gm29581","Gm21539","Gm29578","Gm21780",
            "Gm28486","Gm28487","Gm28488","Gm28491","Gm28490",
            "Gm21562","Gm28298","Gm21572","Gm28297","Gm21723",
            "Gm28296","Gm21821","Gm28295","Gm21842","Gm28170",
            "Gm28764","Gm28762","Gm28763","Gm21094","Gm21588",
            "Gm28761","Gm21599","Gm28765","Gm21916","Gm29191",
            "Gm28210","Gm29074","Gm28518","Gm28519","Gm28520",
            "Gm20838","Gm21626","Gm28521","Gm21633","Gm28522",
            "Gm28517","Gm29532","Gm29531","Gm28464","Gm29656",
            "Gm29655","Gm29653","Gm29654","Gm29080","Gm29081",
            "Gm29077","Gm29078","Gm29225","Gm29226","Gm20855",
            "Gm29386","Gm29384","Gm28260","Gm28264","Gm28463",
            "Gm28585","Gm20896","Gm28470","Gm28469","Gm28332",
            "Gm28333","Gm28336","Gm28338","Gm28811","Gm28810",
            "Gm20835","Gm29027","Gm29028","Gm29023","Gm29025",
            "Gm28212","Gm20905","Gm28820","Gm29671","Gm20738",
            "Gm28545","Gm28547","Gm28771","Gm28772","Gm20897",
            "Gm28774","Gm28775","Gm28176","Gm28786","Gm29360",
            "Gm29497","Gm29498","Gm36782","Gm29839","Gm28595",
            "Gm21865","Gm29433","Gm29215","Gm21801","Gm28850",
            "Gm21914","Gm28617","Gm28619","Gm28993","Gm28994",
            "Gm20795","Gm28427","Gm28426","Gm28425","Gm28834",
            "Gm29537","Gm29274","Gm29275","Gm29276","Gm29645",
            "Gm28312","Gm29042","Gm28311","Gm28310","Gm28568",
            "Gm29286","Gm29285","Gm28280","Gm28278","Gm28279",
            "Gm21739","Gm28886","Gm29426","Gm28887","Gm29444",
            "Gm29445","Gm28457","Gm28458","Gm28461","Gm28460",
            "Gm28851","Gm28456","Gm29648","Gm29646","Gm28549",
            "Gm29210","Gm29213","Gm28550","Gm28604","Gm28565",
            "Gm28566","Gm28697","Gm28554","Gm21256","Gm29056",
            "Gm29098","Gm29373","Gm29370","Gm29166","Gm28194",
            "Gm20792","Gm28754","Gm28753","Gm29070","Gm28553",
            "Gm28726","Gm29270","Gm29271","Gm21209","Gm28681",
            "Gm28684","Gm28679","Gm28870","Gm29511","Gm29117",
            "Gm20883","Gm29116","Gm29321","Gm21117","Gm28317",
            "Gm28318","Gm28313","Gm28315","Gm28316","Gm28944",
            "Gm28945","Gm28947","Gm20920","Gm28948","Gm28950",
            "Gm29063","Gm28197","Gm28735","Gm21258","Gm28732",
            "Gm28733","Gm28233","Gm28235","Gm28236","Gm20747",
            "Gm29250","Gm29252","Gm29248","Gm28908","Gm20929",
            "Gm28276","Gm28274","Gm28291","Gm28293","Gm29305",
            "Gm29303","Gm29557","Gm29301","Gm20931","Gm29555",
            "Gm28355","Gm28352","Gm21858","Gm28201","Gm28202",
            "Gm28206","Gm28207","Gm28208","Gm29222","Gm29221",
            "Gm29220","Gm29405","Gm29406","Gm29404","Gm28964",
            "Sly","Gm28709","Gm28965","ENSMUSG00000121313","Gm21943",
            "Gm28532","ENSMUSG00000121314","Gm28966","Gm28961","Gm28962",
            "Gm28963","ENSMUSG00000121557","Gm29466","Gm29467","ENSMUSG00000121556",
            "Gm20963","Gm28326","ENSMUSG00000121543","Gm28325","Gm28284",
            "Gm29368","Gm21861","Gm21413","Gm28073","Gm21419",
            "Gm28074","Gm21428","Gm28072","Gm29672","Gm29131",
            "Gm29130","Gm20978","Gm29132","Gm21450","Gm28889",
            "Gm28890","Gm28891","Gm28691","Gm28692","Gm28690",
            "Gm28689","Gm28687","Gm20987","1700040F15Rik","Gm29207",
            "Gm29204","Gm29206","Gm29203","Gm21497","Gm29612",
            "Gm29343","Gm29342","Gm20879","Gm28600","Gm21518",
            "Gm20823","ENSMUSG00000121321","Gm29547","Gm29209",
            "Gm29549","Gm29421","Gm21627","Gm29060","Gm29061",
            "Gm29625","Gm29628","ENSMUSG00000121323","Gm29622","Gm29265",
            "Gm20908","Gm28252","ENSMUSG00000121324","Gm28254","Gm20924",
            "Gm29457","Gm20736","Gm28561","Gm29584","ENSMUSG00000121327",
            "Gm28134","Gm28135","ENSMUSG00000121328","Gm28138","Gm28137",
            "Gm28133","ENSMUSG00000121329","Gm29379","Gm29380","Gm29381",
            "Gm20852","Gm29003","ENSMUSG00000121330","Ssty2","Orly",
            "Gm28612","Gm28613","ENSMUSG00000121333","Gm29219","Gm29632",
            "Gm29082","Gm29217","Gm29636","ENSMUSG00000121334","Gm29662",
            "Gm29660","Gm20937","Gm29450","Gm29449","Gm20816",
            "Gm29447","Gm29446","Gm20817","Gm28633","ENSMUSG00000121282",
            "Gm28462","Gm29409","Gm28632","Gm28824","Gm28081",
            "Gm28082","Gm21118","Gm29424","Gm29425","Gm20888",
            "Gm28238","Gm28431","Gm28432","Gm28840","Gm28842",
            "Gm28839","Gm28668","Gm20869","Gm28482","Gm28485",
            "Gm20870","ENSMUSG00000121289","Gm20843","Gm29182","ENSMUSG00000121290",
            "ENSMUSG00000121291","ENSMUSG00000121292","Gm29255","Gm20903","Gm28259",
            "ENSMUSG00000121293","Gm29616","Gm29108","Gm29110","Gm28174",
            "Gm29416","Gm28789","Gm20814","Gm29278","Gm28606",
            "Gm20850","Gm21173","Gm29313","Gm29311","Gm29309",
            "Gm28702","Gm28701","Gm28704","Gm28540","Gm21650",
            "Gm28541","Gm28538","Gm29569","Gm20867","Gm29568",
            "Gm29566","Gm29565","Gm29564","Gm29146","Gm28664",
            "Gm28663","Gm28817","Gm28816","Gm20806","Gm28813",
            "Gm20917","Gm28421","Gm28422","Gm20916","Gm29606",
            "Gm29190","Gm21760","Gm28157","Gm28696","Gm28898",
            "Gm28897","Gm29324","Gm20911","Gm28741","Gm28743",
            "Gm28124","Gm28126","Gm21317","Gm29472","Gm29473",
            "Gm28226","Gm28225","Gm28827","Gm28832","Gm28472",
            "Gm21095","Gm28475","Gm21394","Gm29338","Gm21409",
            "Gm28758","Gm29071","Gm28348","Gm28972","Gm28970",
            "Gm28971","Gm20854","Gm28346","Gm28345","Gm29090",
            "Gm28152","Gm28823","Gm28220","Gm20820","Gm28218",
            "Gm29302","Gm28440","Gm29091","Gm29297","Gm28365",
            "Gm21477","Gm28367","Gm20906","Gm28109","Gm28108",
            "Gm28104","Gm28103","Gm28102","Gm28718","Gm29162",
            "Gm28866","Gm28405","Gm28406","Gm29436","Gm28977",
            "Gm28407","Gm29392","Gm29393","Gm21294","Gm28672",
            "Gm28670","Gm28673","Gm28674","Gm21996","Gm28930",
            "Gm29504","Gm20837","Gm28300","Gm28301","Gm21860",
            "Gm47283")
 
  # build artificial genes
  Xgene.set <-Xgenes[Xgenes %in% row.names(x)]
  Ygene.set <-Ygenes[Ygenes %in% row.names(x)]
  cm.new<-as.data.frame(matrix(rep(0, 3*ncol(x)), ncol = ncol(x),nrow = 3))
  row.names(cm.new) <- c("Xist","superX","superY")
  colnames(cm.new) <- colnames(x)

  if ("Xist" %in% row.names(x)) {
    cm.new["Xist", ]<- x["Xist", ]
  }else{

    cm.new["Xist", ]<- 0
  }

  if (length(Xgene.set)>0){
    cm.new["superX", ] <-colSums(x[Xgene.set,,drop = FALSE])
  }
  if (length(Ygene.set)>0){
    cm.new["superY", ] <-colSums(x[Ygene.set,,drop = FALSE])
  }


  ############################################################################
  # Pre-processing
  # perform simple QC
  # keep a copy of library size
  discarded.cells <- NA
  if (qc == TRUE){
    #data.sce <-SingleCellExperiment(assays = list(counts = x))
    qcstats <- scuttle::perCellQCMetrics(x)
    qcfilter <- scuttle::perCellQCFilters(qcstats)
    # save the discarded cells
    discarded.cells <- colnames(x[,qcfilter$discard])

    # cm.new only contains cells that pass the quality control
    cm.new <-cm.new[,!qcfilter$discard]
  }

  tcm.final <- t(cm.new)
  tcm.final <- as.data.frame(tcm.final)

  #Do Not Classify
  zero.cells <- NA
  dnc <- tcm.final$superY==0 & tcm.final$superX==0

  if(any(dnc)==TRUE){
    zero.cells <- row.names(tcm.final)[dnc]
    message(length(zero.cells), "cell/s are unable to be classified
              due to an abundance of zeroes on X and Y chromosome genes\n")
  }
  tcm.final <- tcm.final[!dnc, ]

  cm.new <- cm.new[,!dnc]

  cm.lib.size<- colSums(x[,colnames(cm.new)], na.rm=TRUE)

  # log-normalisation performed for each cell
  # scaling performed for each gene
  normsca.cm <- data.frame(lognormCounts(cm.new, log = TRUE,
                                         prior.count = 0.5,lib.size=cm.lib.size))
  data.df <- t(normsca.cm)
  data.df <- as.data.frame(data.df)
  row.names(data.df) = row.names(tcm.final)
  return(list(tcm.final=tcm.final, data.df=data.df, discarded.cells=discarded.cells,
       zero.cells=zero.cells))
}

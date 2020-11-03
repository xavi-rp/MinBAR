#' Presences (occurrences) of Linaria alpina
#'
#' A dataset containing the presences (1064) of Linaria alpina in Europe and North Africa.
#' Coord. ref. : +init=EPSG:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0.
#'
#'
#' @format A data frame with 1064 rows and 3 variables.
#' \describe{
#'   \item{decimalLongitude}{DecimalLongitude, in degrees}
#'   \item{decimalLatitude}{DecimalLatitude, in degrees}
#'   \item{species}{Name of the species}
#' }
#' @references
#' GBIF.org (07 March 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.phqgk3.
#'
#' @source \url{https://www.gbif.org/}
"sprecords"


#' Climate variables
#'
#' A raster brick containing 3 climate variables (resolution: 5 minutes) to be used as predictors for modelling species distributions
#'#' Coord. ref. : +init=EPSG:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0.
#'
#' @format A raster brick with 3 variables:
#' \describe{
#'   \item{bio1}{Annual Mean Temperature}
#'   \item{bio7}{Temperature Annual Range}
#'   \item{bio12}{Annual Precipitation}
#' }
#'
#' @examples
#' bioscrop <- raster::brick(paste0(system.file(package='MinBAR'), "/extdata/bioscrop.tif"))
#' names(bioscrop) <- c("bio1", "bio7", "bio12")
#' bioscrop
#'
#' @references
#' Fick, S.E. and R.J. Hijmans, 2017. Worldclim 2: New 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology.
#'
#' @source \url{https://worldclim.org}
"bioscrop"



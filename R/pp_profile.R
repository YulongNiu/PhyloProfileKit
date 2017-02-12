##' Plot profiles
##'
##' Plot phylogenetic profile
##'
##' @title pp plot tiles
##' @param x A character vector
##' @param ... Parameters passed to \code{geom_tile()} in the ggplot2 package.
##' @return A \code{gg} class object
##' @examples
##' pp_tile(rep(letters[1:3], 4))
##' pp_tile(rep(LETTERS[1:2], 3)) + scale_fill_manual(values = c('red','blue', 'green'))
##' pp_tile(rep(letters[1:3], 4)) + coord_flip()
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot geom_tile labs scale_x_continuous scale_y_continuous scale_y_reverse aes_string
##' @seealso \code{\link[ggplot2]{geom_tile}}
##' @export
##' 
pp_profile <- function(x, ...) {

}

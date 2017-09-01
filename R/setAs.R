##' @include PPTreeIdx.R
NULL

## transfer 'PPTreeIdx' to 'PPIdx'
setAs(from = 'PPTreeIdx',
      to = 'PPIdx',
      def = function(from) {
        p <- PPData(from)
        idx <- from@idx
        return(new('PPIdx', p, idx = idx))
      })

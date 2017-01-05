## validate names
setClass(Class = 'pp',
         contains = 'matrix')

tmp1 <- new('pp', matrix(sample(0:1, 10 * 20, replace = TRUE), ncol = 20))
tmp2 <- new('pp', matrix(rnorm(10 * 10), ncol = 20))

## validate nrow and ncol
## validate rownames
setClass(Class = 'idx',
         contains = 'matrix')

tmp3 <- new('idx', matrix(sample(1:10, 3 * 2), ncol = 2))

## validate boundary
## validate integer
setClass(Class = 'ppBin',
         slots = c(ppData = 'pp', selectIdx = 'idx'))

tmp4 <- new('ppBin', ppData = tmp1, selectIdx = tmp3)


## validate boundary
## validate numeric
setClass(Class = 'ppCont',
         slots = c(ppData = 'pp', selectIdx = 'idx'))
tmp5 <- new('ppCont', ppData = tmp2, selectIdx = tmp3)

setClassUnion(name = 'ppIdx',
              members = c('ppBin', 'ppCont'))

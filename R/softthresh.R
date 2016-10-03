####### Soft-threshhold
softsh <- function(u, eta)
{
    if( eta < abs(u) ){
        if( u > 0 ) {
            ss <- u - eta
        } else {
            ss <- u + eta
        }
    } else {
        ss <- 0
    }
    return(ss)
}

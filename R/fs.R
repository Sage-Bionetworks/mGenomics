fs <-
function (x) 
{
    u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
    u$d <- u$d^2/sum(u$d^2)
    u
}

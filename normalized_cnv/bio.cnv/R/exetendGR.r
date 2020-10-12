exetendGR <- function(gr, size=5000){
	ret <- gr
	start(ret) <- start(ret) - size
	end(ret) <- end(ret) + size
	ret
}
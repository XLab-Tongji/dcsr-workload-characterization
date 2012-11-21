##
## Copyright (C) 2008       Marco Guazzone
##                          [Distributed Computing System (DCS) Group,
##                           Computer Science Institute,
##                           Department of Science and Technological Innovation,
##                           University of Piemonte Orientale,
##                           Alessandria (Italy)]
##
## This file is part of dcsr-workload-characterization
##
## dcsr-workload-characterization is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## dcsr-workload-characterization is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with dcsr-workload-characterization  If not, see <http://www.gnu.org/licenses/>.
##

####### DEPRECATED ##########

##
## AUTHORS
##  Marco Guazzone (marco.guazzone@gmail.com)
##

## MG_PLOT_PP
mg_plot_pp <- function(x, qdist=qnorm, probs=NULL, line=TRUE, xlab=NULL, ylab="Probability", ...)
{
    DOTARGS <- as.list(substitute(list(...)))[-1]
    DOTARGS <- paste(names(DOTARGS), DOTARGS, sep="=", collapse=", ")

    xlab=deparse(substitute(x))

    x <- sort(x)
    QNAME <- deparse(substitute(qdist))
    DOTS <- list(...)
    qdist <- match.fun(qdist)
    QFUN <- function(p){
        args=DOTS
        args$p=p
        do.call("qdist", args)
    }
    
    y <- QFUN(ppoints(length(x)))

    if(is.null(probs)){
        probs <- c(.01, .05, seq(.1,.9, by=.1), .95, .99)
        if(length(x)>=1000)
            probs <- c(0.001, probs, .999)
    }

    qprobs <- QFUN(probs)
      
    plot(x, y, axes=FALSE, type="n", ylim=range(c(y,qprobs)), xlab=xlab, ylab=ylab );
    box()
    
    abline(h=qprobs, col="grey")
    axis(1)
    #axis(2, at=qprobs, labels=100*probs)
    axis(2, at=qprobs, labels=probs)

    points(x, y)

    QTEXT <- paste("Quantile: ", QNAME, sep="")
    if(nchar(DOTARGS))
        QTEXT <- paste(QTEXT, DOTARGS, sep=", ")
    mtext(QTEXT, side=1, line=3, adj=1)
    
    xl <- quantile(x, c(0.25, 0.75))
    yl <- qdist(c(0.25, 0.75), ...)
    slope <- diff(yl)/diff(xl)
    int <- yl[1] - slope * xl[1]
    
    if(line){
        abline(int, slope, col="red")
    }
    z <- list(qdist=QFUN, int=int, slope=slope)
    class(z) <- "probplot"
    invisible(z)
}

mg_plot_pp2 <- function(x,pdist,...)
{
	x <- sort(x);
	n <- length(x);
	p <- ((1 : n) - 0.5) / n;
#	p <- ecdf(x)(x);
#	p <- (1 : n) / (n+1);
	dots <- list(...);
	Fn <- function(q){
		args=dots;
		args=list(...);
		args$q=q;
		do.call(pdist, args);
	}
  	y = Fn(x);
	plot(p,y, pch="+", lwd=2,col="royalblue");
	Axis( x=c(0,1), side=1 );
	Axis( x=c(0,1), side=2 );

	l<-lm(y~p);
	abline(l$coefficients[1],l$coefficients[2], lty="dashed");
}

lines.probplot <- function(x, h=NULL, v=NULL, bend=FALSE, ...)
{
    if(is.null(h) & is.null(v)){
        abline(x$int, x$slope, ...)
    }

    pu <- par("usr")

    if(!is.null(h)){
        h <- x$qdist(h)
        if(!bend){
            abline(h=h, ...)
        } else{
            v <- c(v, (h-x$int)/x$slope)
        }
    }

    if(!is.null(v)){
        if(!bend){
            abline(v=v, ...)
        } else{
            h <- v*x$slope+x$int
            segments(v, pu[3], v, h, ...)
            segments(pu[1], h, v, h, ...)
        }
    }
}

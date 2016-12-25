######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################


###################### PsigniSetup ###################################
# Basic data information and fitting directions
#
PsigniSetup <- function (
    x, k, n,
    priors=list("","","Uniform(0,.1)"),
    sigmoid="logistic",
    core="mw0.1",
    number.of.alternatives=2,
    cuts=0.5) {
    # The basic data class that is used for all kinds of fitting

    data <- list (
        stimulus.intensities   = as.double(x),
        number.of.correct      = as.integer(k),
        number.of.trials       = as.integer(n),
        number.of.blocks       = as.integer(length(k)),
        sigmoid                = as.character(sigmoid),
        core                   = as.character(core),
        number.of.alternatives = as.integer(number.of.alternatives),
        priors                 = priors,
        cuts                   = cuts)

    # Store attributes
    attr(data,"class") <- "psignisetup"

    return (data)
}

print.psignisetup <- function ( data ) {

    nprm <- if (data$number.of.alternatives<2) 4 else 3

    # Get parameter names
    parnames <- c("alpha","beta","lambda","gamma")
    if ( substr(data$core,1,2)=="mw" ) {
        parnames[1] = "m    "
        parnames[2] = "w    "
    }

    # Basic
    cat ( paste ("Data from ", data$number.of.alternatives, "-AFC experiment with ", data$number.of.blocks," blocks\n", sep="") )
    cat ( paste ("Intended fit is with",data$sigmoid,"sigmoid and",data$core,"core\n") )
    cat ( "\nPriors:\n" )
    for ( i in 1:nprm ) {
        cat ( paste ( "  ", parnames[i], "\t", data$priors[i], "\n" ) )
    }
    cat ("\n")

    # Cuts
    cat ( paste ( "Requested cuts:", data$cuts[1] ) )
    for ( i in 2:length(data$cuts) ) {
        cat ( paste ( ",", data$cuts[i] ) )
    }
    cat ("\n\n")

    # Data
    print ( data.frame ( stimulus.intensities=data$stimulus.intensities, number.of.correct=data$number.of.correct, number.of.trials=data$number.of.trials ) )
}

############################################################
#            psiginference methods                         #
############################################################

print.psiginference <- function ( inference ) {
    ##############################################
    # Print most important inference information

    # Determine names of parameters
    parnames <- c("alpha","beta","lambda","gamma")
    if ( substr(inference$core,1,2)=="mw" ) {
        parnames[1] = "m    "
        parnames[2] = "w    "
    }

    # Write some header information
    if ( attr(inference,"inference")=="point" ) {
        cat ( "Point estimate of psychometric function\n" )
    } else if ( attr(inference,"inference")=="bootstrap" ) {
        cat ( "Bootstrap inference on psychometric function\n" )
    } else if ( attr(inference,"inference")=="mcmc" ) {
        cat ( "Bayesian inference (MCMC) on psychometric function\n" )
    }
    cat ( paste ( "Fit performed with", inference$sigmoid, "sigmoid and", inference$core, "core\n" ) )

    # Write goodness of fit measures
    # Deviance
    cat ( paste ( "\nDeviance:", inference$deviance ) )
    if ( attr(inference,"inference") == "point" ) {
        cat ( paste ( " (asymptotic upper 95% limit:", inference$number.of.parameters-1,")" ) )
    } else if ( attr(inference,"inference") == "bootstrap" ) {
        cat ( paste ( " (monte-carlo upper 95% limit:", quantile(inference$deviance.samples,.95, na.rm=TRUE), ")" ) )
    } else if ( attr(inference,"inference") == "mcmc" ) {
        cat ( paste ( " (Bayesian p-value for deviance:", mean (
            as.double(inference$deviance.samples<inference$deviance.predictions)[mcmc$burnin:mcmc$number.of.samples], na.rm=TRUE ),")" ) )
    }
    cat ( "\n" )

    # Correlations
    Rpdline <-  paste ( "Rpd =", inference$Rpd )
    Rkdline <-  paste ( "Rkd =", inference$Rkd )
    if ( attr(inference,"inference") == "point" ) {
        Rpdline <- paste ( Rpdline, "\n" )
        Rkdline <- paste ( Rkdline, "\n" )
    } else if ( attr(inference,"inference")=="bootstrap" ) {
        Rpdline <- paste ( Rpdline, "95%-CI: ", paste(quantile( inference$Rpd.samples, c(.025,.975), na.rm=TRUE ),collapse=","), "\n" )
        Rkdline <- paste ( Rkdline, "95%-CI: ", paste(quantile( inference$Rkd.samples, c(.025,.975), na.rm=TRUE ),collapse=","), "\n" )
    } else if ( attr(inference,"inference")=="mcmc") {
        st <- inference$burnin
        sp <- inference$number.of.samples
        Rpdline <- paste ( Rpdline,  "Bayesian p-value:", mean( as.double(inference$Rpd.samples[st:sp]<inference$Rpd.predictions[st:sp]), na.rm=TRUE ), "\n" )
        Rkdline <- paste ( Rkdline,  "Bayesian p-value:", mean( as.double(inference$Rkd.samples[st:sp]<inference$Rkd.predictions[st:sp]), na.rm=TRUE ), "\n" )
    }
    cat ( Rpdline )
    cat ( Rkdline )

    # Parameter estimates
    cat ( "\n\nParameter\tEstimate\t2.5%\t50%\t97.5%\n" )
    for ( i in 1:inference$number.of.parameters ) {
        cat ( paste ( "  ", parnames[i], "\t", round(inference$estimate[i],3), "\t" ) )
        if ( attr(inference,"inference") != "point" ) {
            for ( j in 1:3 ) {
                cat ( paste ( "\t", round(inference$parameter.ci[i,j], 3) ) )
            }
        }
        cat ( "\n" )
    }

    # Threshold estimates
    cat ( "\nThreshold\tEstimate\t2.5%\t50%\t97.5%\n" )
    for ( i in 1:inference$number.of.cuts ) {
        cat ( paste ( "  ", round(inference$cuts[i],2), " \t", round(inference$threshold[i],3), "\t" ) )
        if ( attr(inference,"inference") != "point" ) {
            for ( j in c(1,2,3) ) {
                cat ( paste ( "\t", round(inference$threshold.ci[i,j], 3) ) )
            }
        }
        cat ( "\n" )
    }

    # Influential observations
    if ( attr(inference,"inference") != "point" ) {
        cat ( "\nBlock data and influences\n" )
        d.and.infl <- data.frame (
            stimulus.intensities = inference$stimulus.intensities,
            number.of.correct    = inference$number.of.correct,
            number.of.trials     = inference$number.of.trials,
            influence            = inference$influential
            )
        print(d.and.infl)
    }

}

######################################################################

bootstrap.goodness.of.fit <- function ( inference ) {
    diag <- PsigDiagnostics ( inference$estimate, inference )
    m <- matrix(seq(1,6),2,3,TRUE)
    layout(m)

    ##############################
    # Psychometric function plot
    plot ( inference$stimulus.intensities, inference$number.of.correct/inference$number.of.trials,
        type="n",
        xlab="Stimulus intensity",
        ylab=if(inference$number.of.alternatives<2) "prob(YES)" else "prob(correct)",
        main="Psychometric function",
        ylim=if(inference$number.of.alternatives<2) c(0,1) else c(1./inference$number.of.alternatives,1)
        )
    # Data
    points ( inference$stimulus.intensities, inference$number.of.correct/inference$number.of.trials,
        cex=0.05*inference$number.of.trials )
    # Fitted curve
    psi <- PsigEvaluate ( inference$estimate, inference )
    lines ( psi$x, psi$Psi.x )
    # Confidence limits
    for ( i in 1:inference$number.of.cuts ) {
        lines ( inference$threshold.ci[i,],
            rep( PsigEvaluate(inference$estimate,inference,inference$threshold[i])$Psi.x, length(inference$threshold.ci[i,]) ))
    }

    ##############################
    # Rpd plot
    Psi.x <- PsigEvaluate(inference$estimate, inference, inference$stimulus.intensities)$Psi.x
    plot ( Psi.x, diag$deviance.residuals,
        type="n",
        xlab="Model prediction",
        ylab="Deviance residuals",
        main=paste ("Rpd =", signif(diag$Rpd)) )
    points ( Psi.x, diag$deviance.residuals )
    abline(reg=lm( diag$deviance.residuals~Psi.x ), lty=3 )

    ##############################
    # Rkd plot
    plot ( seq(1,inference$number.of.blocks), diag$deviance.residuals,
        type="n",
        xlab="Block index",
        ylab="Deviance residuals",
        main=paste ("Rkd =", signif(diag$Rkd)) )
    points ( seq(1,inference$number.of.blocks), diag$deviance.residuals )
    abline(reg=lm( diag$deviance.residuals~seq(1,inference$number.of.blocks) ), lty=3 )

    # Deviance histogram
    h <- hist ( inference$deviance.samples, xlab="deviance", main="Expected deviance" )
    abline ( v=inference$deviance, col="red" )
    abline ( v=quantile(inference$deviance.samples, .95 ), lty=3)

    # Rpd histogram
    h <- hist ( inference$Rpd.samples, xlab="Rpd", main="Expected Rpd", xlim=c(-1,1) )
    abline ( v=inference$Rpd, col="red" )
    abline ( v=quantile(inference$Rpd.samples,c(.025,.975),na.rm=TRUE), lty=3 )

    # Rpd histogram
    h <- hist ( inference$Rkd.samples, xlab="Rkd", main="Expected Rkd", xlim=c(-1,1) )
    abline ( v=inference$Rkd, col="red" )
    abline ( v=quantile(inference$Rkd.samples,c(.025,.975),na.rm=TRUE), lty=3 )

}

######################################################################

bayes.goodness.of.fit <- function ( inference ) {
    diag <- PsigDiagnostics ( inference$estimate, inference )
    m <- matrix(seq(1,6),2,3,TRUE)
    layout(m)

    ##############################
    # Psychometric function plot
    plot ( inference$stimulus.intensities, inference$number.of.correct/inference$number.of.trials,
        type="n",
        xlab="Stimulus intensity",
        ylab=if(inference$number.of.alternatives<2) "prob(YES)" else "prob(correct)",
        main="Psychometric function",
        ylim=if(inference$number.of.alternatives<2) c(0,1) else c(1./inference$number.of.alternatives,1)
        )

    # Some posterior samples
    for ( i in 1:20 ) {
        index.parameters <- round ( runif(1, min=mcmc$burnin,max=mcmc$number.of.samples) )
        psi <- PsigEvaluate ( inference$parameter.samples[index.parameters,], inference )
        lines ( psi$x, psi$Psi.x, col="gray70" )
    }

    # Data
    points ( inference$stimulus.intensities, inference$number.of.correct/inference$number.of.trials,
        cex=0.05*inference$number.of.trials )
    # Fitted curve
    psi <- PsigEvaluate ( inference$estimate, inference )
    lines ( psi$x, psi$Psi.x )
    for ( i in 1:inference$number.of.cuts ) {
        lines ( inference$threshold.ci[i,],
            rep( PsigEvaluate(inference$estimate,inference,inference$threshold[i])$Psi.x, length(inference$threshold.ci[i,]) ))
    }

    ################################
    # Rpd plot
    Psi.x <- PsigEvaluate(inference$estimate, inference, inference$stimulus.intensities)$Psi.x
    plot ( Psi.x, diag$deviance.residuals,
        type="n",
        xlab="Model prediction",
        ylab="Deviance residuals",
        main=paste ("Rpd =", signif(diag$Rpd)) )
    points ( Psi.x, diag$deviance.residuals )
    abline(reg=lm( diag$deviance.residuals~Psi.x ), lty=3 )

    ################################
    # Rkd plot
    plot ( seq(1,inference$number.of.blocks), diag$deviance.residuals,
        type="n",
        xlab="Block index",
        ylab="Deviance residuals",
        main=paste ("Rkd =", signif(diag$Rkd)) )
    points ( seq(1,inference$number.of.blocks), diag$deviance.residuals )
    abline(reg=lm( diag$deviance.residuals~seq(1,inference$number.of.blocks) ), lty=3 )


    # Deviance scatter
    st <- inference$burnin
    sp <- inference$number.of.samples
    plot ( inference$deviance.predictions[st:sp], inference$deviance.samples[st:sp],
        type="n",
        xlab="predicted deviance",
        ylab="observed deviance",
        main="Posterior predictions for deviance" )
    points ( inference$deviance.predictions[st:sp], inference$deviance.samples[st:sp] )
    abline(a=0,b=1,lty=3)

    # Rpd scatter
    plot ( inference$Rpd.predictions[st:sp], inference$Rpd.samples[st:sp],
        type="n",
        xlab="predicted Rpd",
        ylab="observed Rpd",
        main="Posterior predictions for Rpd" )
    points ( inference$Rpd.predictions[st:sp], inference$Rkd.samples[st:sp] )
    abline(a=0,b=1,lty=3)

    # Rkd scatter
    plot ( inference$Rkd.predictions[st:sp], inference$Rkd.samples[st:sp],
        type="n",
        xlab="predicted Rkd",
        ylab="observed Rkd",
        main="Posterior predictions for Rkd" )
    points ( inference$Rkd.predictions[st:sp], inference$Rkd.samples[st:sp] )
    abline(a=0,b=1,lty=3)
}

######################################################################

plot.psiginference <- function ( inference ) {
    if ( attr(inference,"inference") == "bootstrap" ) {
        bootstrap.goodness.of.fit ( inference )
    } else if ( attr(inference,"inference") == "bayes" ) {
        bayes.goodness.of.fit ( inference )
    }
}


############################################################
#             psiginference types                          #
############################################################

MAPestimation <- function ( psignidata ) {

    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    # Do the real shit
    map <- .C ( "mapestimate",
        stimulus.intensities=as.double(psignidata$stimulus.intensities),
        number.of.correct=as.integer(psignidata$number.of.correct),
        number.of.trials=as.integer(psignidata$number.of.trials),
        number.of.blocks=as.integer(psignidata$number.of.blocks),
        sigmoid=as.character(psignidata$sigmoid),
        core=as.character(psignidata$core),
        number.of.alternatives=as.integer(psignidata$number.of.alternatives),
        priors=as.character(psignidata$priors),
        number.of.parameters=as.integer(nprm),
        estimate=as.double( vector("numeric", nprm) ),
        deviance=as.double(0),
        Rpd=as.double(0),
        Rkd=as.double(0)
        )

    # Most of the items don't make sense for this
    map$data.samples <- NULL
    map$parameter.samples <- NULL
    map$deviance.samples <- NULL
    map$Rpd.samples <- NULL
    map$Rkd.samples <- NULL
    map$threshold.samples <- NULL
    map$influential <- NULL
    map$acceleration <- NULL
    map$bias <- NULL
    map$parameter.ci <- NULL

    # Determine thresholds
    map$cuts <- psignidata$cuts
    map$number.of.cuts <- length(psignidata$cuts)
    map$threshold <- PsigDiagnostics( map$estimate, psignidata )$threshold

    # Set attributes
    attr(map,"class") <- "psiginference"
    attr(map,"inference") <- "point"

    return (map)
}

######################################################################

PsigBootstrap <- function ( psignidata, number.of.samples=2000, generating=-999 ) {

    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    # Do the real shit
    boots <- .C ( "performbootstrap",
        stimulus.intensities   = as.double(psignidata$stimulus.intensities),
        number.of.correct      = as.integer(psignidata$number.of.correct),
        number.of.trials       = as.integer(psignidata$number.of.trials),
        number.of.blocks       = as.integer(psignidata$number.of.blocks),
        sigmoid                = as.character(psignidata$sigmoid),
        core                   = as.character(psignidata$core),
        number.of.alternatives = as.integer(psignidata$number.of.alternatives),
        priors                 = as.character(psignidata$priors),
        generating             = as.double(generating),
        number.of.parameters   = as.integer(nprm),
        number.of.samples      = as.integer(number.of.samples),
        cuts                   = as.double(psignidata$cuts),
        number.of.cuts         = as.integer(length(psignidata$cuts)),
        data.samples           = as.integer(vector("numeric",psignidata$number.of.blocks*number.of.samples)),
        parameter.samples      = as.double(vector("numeric",nprm*number.of.samples)),
        deviance.samples       = as.double(vector("numeric",number.of.samples)),
        Rpd.samples            = as.double(vector("numeric",number.of.samples)),
        Rkd.samples            = as.double(vector("numeric",number.of.samples)),
        threshold.samples      = as.double(vector("numeric",length(psignidata$cuts)*number.of.samples)),
        influential            = as.double(vector("numeric",psignidata$number.of.blocks)),
        acceleration           = as.double(vector("numeric",length(psignidata$cuts))),
        bias                   = as.double(vector("numeric",length(psignidata$cuts))),
        threshold.ci           = as.double(vector("numeric",length(psignidata$cuts)*3))
        )

    # Convert sample arrays to matrices that have meaningful dimensions
    boots$data.samples <- matrix(boots$data.samples, boots$number.of.samples,boots$number.of.blocks, TRUE)
    boots$parameter.samples <- matrix(boots$parameter.samples, boots$number.of.samples, boots$number.of.parameters, TRUE)
    boots$threshold.samples <- matrix(boots$threshold.samples, boots$number.of.samples, boots$number.of.cuts, TRUE)
    boots$threshold.ci <- matrix(boots$threshold.ci, boots$number.of.cuts, 3)

    # Some parameters are determined from constrained ML which is the same as MAP
    map <- MAPestimation ( psignidata )
    boots$estimate <- map$estimate
    boots$deviance <- map$deviance
    boots$Rpd <- map$Rpd
    boots$Rkd <- map$Rkd
    boots$threshold <- map$threshold

    if ( generating==-999 ) {
        # We did a nonparametric bootstrap but the goodness of fit data should be parametric in any case
        boots.parametric <- .C ( "performbootstrap",
            stimulus.intensities   = as.double(psignidata$stimulus.intensities),
            number.of.correct      = as.integer(psignidata$number.of.correct),
            number.of.trials       = as.integer(psignidata$number.of.trials),
            number.of.blocks       = as.integer(psignidata$number.of.blocks),
            sigmoid                = as.character(psignidata$sigmoid),
            core                   = as.character(psignidata$core),
            number.of.alternatives = as.integer(psignidata$number.of.alternatives),
            priors                 = as.character(psignidata$priors),
            generating             = as.double(map$estimate),
            number.of.parameters   = as.integer(nprm),
            number.of.samples      = as.integer(number.of.samples),
            cuts                   = as.double(psignidata$cuts),
            number.of.cuts         = as.integer(length(psignidata$cuts)),
            data.samples           = as.integer(vector("numeric",psignidata$number.of.blocks*number.of.samples)),
            parameter.samples      = as.double(vector("numeric",nprm*number.of.samples)),
            deviance.samples       = as.double(vector("numeric",number.of.samples)),
            Rpd.samples            = as.double(vector("numeric",number.of.samples)),
            Rkd.samples            = as.double(vector("numeric",number.of.samples)),
            threshold.samples      = as.double(vector("numeric",length(psignidata$cuts)*number.of.samples)),
            influential            = as.double(vector("numeric",psignidata$number.of.blocks)),
            acceleration           = as.double(vector("numeric",length(psignidata$cuts))),
            bias                   = as.double(vector("numeric",length(psignidata$cuts))),
            threshold.ci           = as.double(vector("numeric",length(psignidata$cuts)*3))
            )
        boots$deviance.samples <- boots.parametric$deviance.samples
        boots$Rpd.samples      <- boots.parametric$Rpd.samples
        boots$Rkd.samples      <- boots.parametric$Rkd.samples
    }

    # These are not used
    boots$logratios <- NULL
    boots$proposal <- NULL
    boots$deviance.predictions <- NULL

    # Determine parameter confidence intervals as quantiles
    boots$parameter.ci <- matrix ( nrow=boots$number.of.parameters, ncol=3 )
    for ( i in 1:nprm ) {
        boots$parameter.ci[i,] = quantile ( boots$parameter.samples[,i], c(.05,.5,.95) )
    }

    # Store attributes
    attr(boots,"class") <- "psiginference"
    attr(boots,"inference") <- "bootstrap"

    return (boots)
}

PsigBayes <- function ( psignidata, number.of.samples=2000, start=NULL, proposal=NULL ) {

    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    # Check proposal dist
    if (is.null(proposal)) {
        # Use Raftery/Lewis instead? Would have to be implemented...
        proposal <- if (nprm==4) c(.4,.4,.01,.01) else c(.4,.4,.01)
    } else {
        if (length(proposal)<nprm) {
            cat ("Error in PsigBayes: Wrong length of proposal argument\n")
            return (NULL)
        }
    }

    # Check starting values
    if (is.null(start)) {
        start <- MAPestimation(psignidata)$estimate
    }

    # Do the real work
    mcmc <- .C ( "performmcmc",
        stimulus.intensities   = as.double(psignidata$stimulus.intensities),
        number.of.correct      = as.integer(psignidata$number.of.correct),
        number.of.trials       = as.integer(psignidata$number.of.trials),
        number.of.blocks       = as.integer(psignidata$number.of.blocks),
        sigmoid                = as.character(psignidata$sigmoid),
        core                   = as.character(psignidata$core),
        number.of.alternatives = as.integer(psignidata$number.of.alternatives),
        priors                 = as.character(psignidata$priors),
        proposal               = as.double(proposal),
        generating             = as.double(start),
        number.of.parameters   = as.integer(nprm),
        number.of.samples      = as.integer(number.of.samples),
        cuts                   = as.double(psignidata$cuts),
        number.of.cuts         = as.integer(length(psignidata$cuts)),
        parameter.samples      = as.double(vector("numeric",nprm*number.of.samples)),
        deviance.samples       = as.double(vector("numeric",number.of.samples)),
        Rpd.samples            = as.double(vector("numeric",number.of.samples)),
        Rkd.samples            = as.double(vector("numeric",number.of.samples)),
        data.samples           = as.integer(vector("numeric",psignidata$number.of.blocks*number.of.samples)),
        Rpd.predictions        = as.double(vector("numeric",number.of.samples)),
        Rkd.predictions        = as.double(vector("numeric",number.of.samples)),
        deviance.predictions   = as.double(vector("numeric",number.of.samples)),
        logratios              = as.double(vector("numeric",psignidata$number.of.blocks*number.of.samples))
        )

    # Convert arrays to matrices with meaningful dimensions
    mcmc$data.samples <- matrix(mcmc$data.samples, mcmc$number.of.samples, mcmc$number.of.blocks, TRUE)
    mcmc$parameter.samples <- matrix(mcmc$parameter.samples, mcmc$number.of.samples, mcmc$number.of.parameters, TRUE)
    mcmc$logratios <- matrix(mcmc$logratios, mcmc$number.of.samples, mcmc$number.of.blocks, TRUE)

    # Set some default values for burnin and thinning
    mcmc$burnin <- as.integer(0.5*mcmc$number.of.samples)
    mcmc$thin   <- 1

    # Determine mean estimate
    mcmc$estimate <- apply ( mcmc$parameter.samples[mcmc$burnin:mcmc$number.of.samples,], 2, mean )

    # Get some diagnostics for the the mean estimate
    diag <- PsigDiagnostics( mcmc$estimate, psignidata )
    mcmc$deviance <- diag$deviance
    mcmc$threshold <- diag$threshold
    mcmc$Rpd <- diag$Rpd
    mcmc$Rkd <- diag$Rkd

    # Determine parameter confidence intervals as quantiles
    mcmc$parameter.ci <- matrix( nrow=mcmc$number.of.parameters, ncol=3 )
    for ( i in 1:nprm ) {
        mcmc$parameter.ci[i,] = quantile ( mcmc$parameter.samples[,i], c(.05,.5,.95) )
    }

    # Determine threshold confidence intervals
    prm2thres <- function ( parameters, psignidata ) {
        # Will perform the conversion to thresholds
        return ( PsigDiagnostics ( parameters, psignidata )$threshold )
    }
    # Get threshold samples
    mcmc$threshold.samples <- t(apply(mcmc$parameter.samples,1,prm2thres,D))
    # Now create the threshold matrix
    mcmc$threshold.ci <- matrix(nrow=mcmc$number.of.cuts,ncol=3)
    for ( i in 1:mcmc$number.of.cuts ) {
        mcmc$threshold.ci[i,] <- quantile ( mcmc$threshold.samples[mcmc$burnin:mcmc$number.of.samples,i], c(0.025, 0.5, 0.975) )
    }

    # Determine influential observations
    lpr <- mcmc$logratios[mcmc$burnin:mcmc$number.of.samples,]
    KL_determine <- function ( logratios ) {
        return ( -mean(logratios) + log(mean(exp(logratios))) )
    }
    mcmc$influential = apply ( lpr, 2, KL_determine )


    # Store attributes
    attr(mcmc,"class") <- "psiginference"
    attr(mcmc,"inference") <- "mcmc"

    return (mcmc)
}

###################################################
# Diagnostic tools                                #
###################################################

PsigDiagnostics <- function ( parameters, psignidata ) {

    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    # Work it
    diag <- .C ( "getdiagnostics",
        stimulus.intensities   = as.double(psignidata$stimulus.intensities),
        number.of.correct      = as.integer(psignidata$number.of.correct),
        number.of.trials       = as.integer(psignidata$number.of.trials),
        number.of.blocks       = as.integer(psignidata$number.of.blocks),
        sigmoid                = as.character(psignidata$sigmoid),
        core                   = as.character(psignidata$core),
        number.of.alternatives = as.integer(psignidata$number.of.alternatives),
        priors                 = as.character(psignidata$priors),
        number.of.parameters   = as.integer(nprm),
        cuts                   = as.double(psignidata$cuts),
        number.of.cuts         = as.integer(length(psignidata$cuts)),
        parameters             = as.double(parameters),
        deviance               = as.double(0),
        Rpd                    = as.double(0),
        Rkd                    = as.double(0),
        threshold              = as.double(vector("numeric",length(psignidata$cuts))),
        deviance.residuals     = as.double(vector("numeric",psignidata$number.of.blocks))
        )

    # Store the diagnostics
    attr(diag,"class") <- "psigdiagnostics"

    return (diag)
}

######################################################################

PsigEvaluate <- function ( parameters, psignidata, x=NULL ) {

    # Create x if non is supplied
    if (is.null(x)) x <- seq(min(psignidata$stimulus.intensities), max(psignidata$stimulus.intensities), length.out=100)

    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    # Work it
    Fx <- .C ( "pmfevaluate",
        stimulus.intensities   = as.double(x),
        number.of.intensities  = as.integer(length(x)),
        parameters             = as.double(parameters),
        number.of.parameters   = as.integer(nprm),
        sigmoid                = as.character(psignidata$sigmoid),
        core                   = as.character(psignidata$core),
        number.of.alternatives = as.integer(psignidata$number.of.alternatives),
        f.x                    = as.double(vector("numeric",length(x)))
        )

    return (list(x=x,Psi.x=Fx$f.x))
}


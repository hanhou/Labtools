#!/usr/bin/env python

from optparse import OptionParser
from pypsignifit import psigobservers
from pypsignifit import psigcorrect
import pypsignifit
from numpy import array,arange,mgrid,clip,zeros,sort, random,mean, asarray
import os,sys
import pylab
from pylab import prctile as ptile
import operator
import time
import subprocess

__helptext__ = """
Determine coverage of confidence intervals for a given combination of analysis/generating parameters.
"""

# Analyze command line options
parser = OptionParser ( usage=__helptext__ )

# Generating parameters
parser.add_option ( "--gen-observer", dest="observer", default="binomial",
        help="model for the variance of the observers responses. Valid choices are 'binomial', 'betabinomial', or 'learning'. Default: 'binomial'" )
parser.add_option ( "--gen-observer-params", dest="observer_params", default="",
        help="parameters of the variance model for the generating observer. If the observer is 'binomial', no parameter is required, the 'betabinomial' observer requires an additional parameter parameter m given the dispersion of the beta distribution used to model the observers response variance. If the 'learning' observer is selected, a string with the following parameters should be given: 'L,r', where L is the learning range both prm1 and prm2 vary from gen_prm+L to gen_prm-L. r is the learning rate. That is, both parameters vary from one trial to the next as prm(n+1)=prm(n)+(prm_gen-L-prm(n))/r. Thus, low values of r indicate faster learning. Default: '' (binomial observer), '10' (betabinomial observer), '0.7,40' (learning observer)" )
parser.add_option ( "--gen-sigmoid",  dest="gen_sigmoid", default="logistic",
        help="sigmoid used when generating the data. Default: 'logistic'" )
parser.add_option ( "--gen-core",     dest="gen_core",    default="mw0.1",
        help="core type used when generating the data. Default: 'mw0.1'" )
parser.add_option ( "--gen-prm1",     dest="gen_prm1",    default=4,    type="float",
        help="first parameter of the psychometric function used to generate the data (typically called alpha, or m). Default: 4" )
parser.add_option ( "--gen-prm2",     dest="gen_prm2",    default=2,    type="float",
        help="second parameter of the psychometric function used to generate the data (typically called beta, or w). Default: 2" )
parser.add_option ( "--gen-prm3",     dest="gen_prm3",    default=0.02, type="float",
        help="third parameter of the psychometric function used to generate the data (this is the lapse rate). Default: 0.02" )
parser.add_option ( "--gen-prm4",     dest="gen_prm4",    default=0.02, type="float",
        help="fourth parameter of the psychometric function used to generate the data (this is the guessing rate if a yes-no paradigm is employed). In nAFC, this parameter is not used. Default: 0.02" )
parser.add_option ( "--gen-nafc",     dest="gen_nafc",    default=2,    type="int",
        help="number of alternatives in the generating task (1 indicates a yes-no task). Default: 2" )

# Analyzing parameters
parser.add_option ( "--ana-sigmoid",  dest="ana_sigmoid", default="logistic",
        help="sigmoid used when analyzing the data. Default: logistic" )
parser.add_option ( "--ana-core",     dest="ana_core",    default="mw0.1",
        help="core type used when analyzing the data. Default: 'mw0.1'" )
parser.add_option ( "--ana-nafc",     dest="ana_nafc",    default=None,    type="int",
        help="number of alternatives assumed when analyzing the data. Default is to take the same as gen_nafc" )
parser.add_option ( "--nbootstrap",   dest="nbootstrap",  default=2000, type="int",
        help="number of bootstrap repetitions when determining confidence intervals or goodness of fit statistics. Default: 2000" )

# constraints for bootstrap
parser.add_option ( "--constraint-prm1", dest="constraint_prm1",
        default="unconstrained",
        help="Parameter constraint on first parameter: for bootstrap inference")

parser.add_option ( "--constraint-prm2", dest="constraint_prm2",
        default="unconstrained",
        help="Parameter constraint on second parameter: for bootstrap inference")

parser.add_option ( "--constraint-prm3", dest="constraint_prm3",
        default="Uniform(0,0.1)",
        help="Parameter constraint on third parameter: for bootstrap inference")

parser.add_option ( "--constraint-prm4", dest="constraint_prm4",
        default="Uniform(0,0.1)",
        help="Parameter constraint on fourth parameter: for bootstrap inference")

# priors for Bayes
parser.add_option ( "--prior-prm1", dest="prior_prm1",
        default="Gauss(0,100)",
        help="Parameter prior on first parameter: for Bayesian inference")

parser.add_option ( "--prior-prm2", dest="prior_prm2",
        default="Gamma(1.01, 2000)",
        help="Parameter prior on second parameter: for Bayesian inference")

parser.add_option ( "--prior-prm3", dest="prior_prm3",
        default="Beta(2,50)",
        help="Parameter prior on third parameter: for Bayesian inference")

parser.add_option ( "--prior-prm4", dest="prior_prm4",
        default="Beta(2,50)",
        help="Parameter prior on fourth parameter: for Bayesian inference")

# Simulation options
parser.add_option ( "--nsimulations", dest="nsimulations", default=1000, type="int",
        help="number of full simulation runs used to determine coverage information. Default: 1000" )
parser.add_option ( "--blocksize",    dest="blocksize",    default=10,   type="int",
        help="number of trials per block in the simulated experiment. Default: 10" )
parser.add_option ( "--nblocks",      dest="nblocks",      default=5,    type="int",
        help="number of blocks in the simulated experiment. Default: 5" )
parser.add_option ( "--fixed-levels", dest="fixed_levels", default=None,
        type="string", help="list of stimulus levels, if desired, if None,"+\
                "psychometric function will be sampled. Default: None")
parser.add_option ( "--fixed-sequence", dest="fixed_sequence", action="store_true",
        help="if this is set, the psychometric function will be sampled at fixed levels" )
parser.add_option ( "--fixed-pmf", dest="fixed_pmf", action="store_true",
        help="if this is set, samples will be generated from the true pmf" )
parser.add_option ( "--seed", dest="seed", default="fixed",
        type="string", help="seed for simulation, can be 'fixed, 'time', or an"+\
                "integer value.")

# output options
parser.add_option ( "-o", "--output", dest="outputfile", default="test",
        help="name of the output file in which the data should be stored,"+\
        "'.data' will be appended to the filename." )

parser.add_option ( "--metadata", action="store_true", default=False,
        help="generate metadata file, filename will the name of the "+\
        "outputfile with the suffix '.meta'.")

#parser.add_option ( "--datareduce", dest="datareduce", action="store_true",
#        help="reduce data based on the estimated nu parameter" )

parser.add_option ( "--disable-nonparametric", dest="nonparametric",
        action="store_false", default=True,
        help="do not run the non-paramteric bootstrap")

parser.add_option ( "--disable-parametric", dest="parametric",
        action="store_false", default=True,
        help="do not run the paramteric bootstrap")

parser.add_option ( "--disable-bayes", dest="bayes",
        action="store_false", default=True,
        help="do not run the bayesian analysis")




options,arguments = parser.parse_args()


############################################################
#                                                          #
# Create useful variables from command line arguments      #
#                                                          #
############################################################

# Observer fixed parameters
if options.gen_nafc == 1:
    gen_prm = [ options.gen_prm1, options.gen_prm2, options.gen_prm3, options.gen_prm4 ]
elif options.gen_nafc > 1:
    gen_prm = [ options.gen_prm1, options.gen_prm2, options.gen_prm3 ]
else:
    raise IOError, "gen_nafc should be > 0, but is %d" % ( options.gen_nafc, )

if options.ana_nafc is None:
    options.ana_nafc = options.gen_nafc

# Observer variable parameters
gen_kwargs = { "nafc": options.gen_nafc,
        "sigmoid": options.gen_sigmoid,
        "core":    options.gen_core }

# Create the desired Observer
if options.observer == "binomial":

    def create_new_observer ():
        return psigobservers.Observer ( *gen_prm, **gen_kwargs )

elif options.observer == "betabinomial":

    if options.observer_params=="":
        M = 10
    else:
        M = float(options.observer_params)
    def create_new_observer ():
        return psigobservers.BetaBinomialObserver ( *(gen_prm+[M]), **gen_kwargs )

elif options.observer == "learning":

    if options.observer_params == "":
        L,r = .7,40.
    else:
        L,r = [float(x) for x in options.observer_params.split(",")]
    end1 = options.gen_prm1 - L
    end2 = options.gen_prm2 - L
    start_prm = list(gen_prm)
    start_prm[0] += L
    start_prm[1] += L
    def create_new_observer ():
        return psigobservers.LinearSystemLearner ( *(start_prm+[{'a': (r,end1), 'b': (r,end2)}]), **gen_kwargs )

else:
    raise IOError, "Invalid observer model: %s" % ( options.observer, )

# Estimation related parameters
ana_kwargs = { "nafc": options.gen_nafc,
        "sigmoid":     options.ana_sigmoid,
        "core":        options.ana_core }

if not options.ana_nafc is None:
    ana_kwargs["nafc"] = options.ana_nafc

# Create observer
OO = psigobservers.Observer ( *gen_prm, **gen_kwargs )

# Create stimulus levels
if options.fixed_levels is not None:
    # Use the stimulus levels that were given on the command line
    message = "'fixed-levels' must be a sequence of numbers, of length nblocks."
    try:
        x = eval(options.fixed_levels)
    except Exception:
        raise ValueError(message + options.fixed_levels)
    if not operator.isSequenceType(x) or False in [operator.isNumberType(i) for i in x]:
        raise ValueError(message + options.fixed_levels)
    if not len(x) == options.nblocks:
        raise ValueError("Argument mismatch: You gave "+str(len(x))+" levels on the command line"+\
                " but specified "+str(options.nblocks)+" blocks.")
    print "levels:", x
else:
    # Create stimuli (roughly based on results from Hills PhD thesis in chapter 5)
    Fx = mgrid[.1:.99:1j*(options.nblocks)]
    # Fx -= Fx.mean()
    # Fx /= Fx.std()
    # Fx *= 0.22    # this corresponds to Hills PhD thesis (sigma_p ~ 0.11 in 2AFC)
    # if options.ana_nafc == 1:
    #     Fx += 0.5
    # else:
    #     Fx += 0.6   # This corresponds to Hills PhD thesis (p_bar ~ 0.8 in 2AFC)
    # Fx = array( Fx.tolist()+[.99] )
    # Fx = clip(Fx,.001,.999)
    x = OO.getlevels(Fx)
    y = mgrid[0.001:0.999:100j]
    print x,Fx

# Constraints and Priors
constraints = [options.constraint_prm1,
        options.constraint_prm2,
        options.constraint_prm3]
# priors      = ["Gauss(4,.1)", "Gamma(1,4)","Beta(2,50)"]   # Hilft auch nicht so viel
# priors      = ["Gauss(0,100)", "Gamma(1.01,2000)","Beta(2,50)"]
priors = [options.prior_prm1, options.prior_prm2, options.prior_prm3]
if options.ana_nafc < 2:
    constraints.append(options.constraint_prm4)
    priors.append(options.prior_prm4)
    #priors += ["Beta(1,10)"]
print priors


# Parse and set seed

def check_int(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

if options.seed not in ["fixed", "time"] and not check_int(options.seed):
    raise ValueError("'seed' must be either 'fixed', 'time' or an integer value.")
elif options.seed == 'fixed':
    print "Seed is default."
elif options.seed == 'time':
    seed = int(time.time())
    print "Seed is time since epoch in seconds: '%d'" % seed
    pypsignifit.set_seed(seed)
    options.seed = seed
else:
    seed = int(options.seed)
    print "Seeed is value given on command line: '%d'" % seed
    pypsignifit.set_seed(seed)

# check that there is something to do
bayes = options.bayes
nonparametric = options.nonparametric
parametric = options.parametric
if not nonparametric and not parametric and not bayes:
    raise ValueError("You must specify one of: 'nonparametric', 'parametric' "+\
            "'bayes' in order for this script to do anything!")

# set the version
options.version = pypsignifit.version

# Organize output
def get_yes_no():
    while True:
        decision = raw_input ( " Overwrite? [Y/N] " )
        if decision.upper() == "Y":
            return True
        elif decision.upper() == "N":
            return False
        else:
            print "'%s' is not a valid option, please answer [Y/N]." % decision
outputfile = options.outputfile + ".data"
metadatafile = options.outputfile + ".meta" if options.metadata else None
if os.path.exists( outputfile ):
    sys.stderr.write ( "Output file %s exists.\n" %(outputfile,) )
    if metadatafile and os.path.exists( metadatafile ):
        sys.stderr.write ( "Metadata file %s exists.\n" %(metadatafile,) )
    if get_yes_no():
       pass
    else:
        sys.stderr.write ("Terminating\n")
        sys.exit(0)
print "writing output to", outputfile
outfile = open(outputfile, "w")
if metadatafile:
    print "writing metadata to", metadatafile
    open(metadatafile, "w").write("%s\n" % str(options))

############################################################
#                                                          #
# Perform the simulation                                   #
#                                                          #
############################################################

def getcpe ( est, mc ):
    mc = mc.copy()
    mc.sort()
    return mean ( mc<est )

def check_ci ( observer, inference ):
    true_thres = observer.thres
    ci = inference.getCI ( 1, (.025, .975) )
    if ci[0]<true_thres and ci[1]>true_thres:
        return 1
    else:
        return 0

def write_header(f):
    outs = "run m.gen w.gen "
    if nonparametric:
        outs += ("m.npr.e m.npr.l m.npr.h "+
                 "w.npr.e w.npr.l w.npr.h "+
                 "d.npr d.npr.crit nu.npr "+
                 "rpd.npr rpd.npr.l rpd.npr.h "+
                 "rkd.npr rkd.npr.l rkd.npr.h "+
                 "infl.npr." \
                + " infl.npr.".join([str(x) for x in range(options.nblocks)]) + " " \
                + "npr.time ")
    if parametric:
        outs += ("m.par.e m.par.l m.par.h "+
                 "w.par.e w.par.l w.par.h "+
                 "d.par d.par.crit nu.par "+
                 "rpd.par rpd.par.l rpd.par.h "+
                 "rkd.par rkd.par.l rkd.par.h "+
                 "infl.par." \
                + " infl.par.".join([str(x) for x in range(options.nblocks)]) + " " \
                + "par.time ")
    if bayes:
        outs += ("m.bay.e "+
                "m.bay.map "+
                "m.bay.median "+
                "m.bay.l "+
                "m.bay.h "+
                "w.bay.e "+
                "w.bay.map "+
                "w.bay.median "+
                "w.bay.l "+
                "w.bay.h "+
                "d.bay "+
                "d.bay.p "+
                "nu.bay "+
                "rpd.bay "+
                "rpd.bay.p "+
                "rkd.bay "+
                "rkd.bay.p "+
                "conv.bay "+
                "Rhat.0 "+
                "Rhat.1 "+
                "Rhat.2 "+
                "infl.bay." \
                + " infl.bay.".join([str(x) for x in range(options.nblocks)])+" " +
                "bay.time ")
    outs += " stim." + " stim.".join([str(x) for x in range(options.nblocks)])
    outs += " resp." + " resp.".join([str(x) for x in range(options.nblocks)])
    outs += "\n"
    f.write(outs)
    f.flush()

def write_id_gen(f, simulation):
    outs = "%d %g %g " % (simulation, gen_prm[0], gen_prm[1] )
    f.write ( outs )

def write_nonparametric(f, Bnpr, runtime):
    outs  = "%g %g %g " % (Bnpr.estimate[0],Bnpr.getCI(1,(.025,)),Bnpr.getCI(1,(.975))) # m.npr.e m.npr.l m.npr.h
    outs += "%g %g %g " % (Bnpr.estimate[1],ptile(Bnpr.mcestimates[:,1],2.5),ptile(Bnpr.mcestimates[:,1],97.5)) # w.npr.e w.npr.l w.npr.h
    outs += "%g %g %g " % (Bnpr.deviance, ptile(Bnpr.mcdeviance,95), psigcorrect.estimate_nu (Bnpr)[0]) # d.npr d.npr.crit nu.npr
    outs += "%g %g %g " % (Bnpr.Rpd,ptile(Bnpr.mcRpd,2.5),ptile(Bnpr.mcRpd,97.5)) # rpd.npr rpd.npr.l rpd.npr.h
    outs += "%g %g %g " % (Bnpr.Rkd,ptile(Bnpr.mcRkd,2.5),ptile(Bnpr.mcRkd,97.5)) # rkd.npr rkd.npr.l rkd.npr.h
    outs += ("%g "*options.nblocks) % tuple(Bnpr.infl)
    outs += "%g " % runtime
    f.write ( outs )

def write_parametric(f, Bpar, runtime):
    outs  = "%g %g %g " % (Bpar.estimate[0],Bpar.getCI(1,(.025,)),Bpar.getCI(1,(.975))) # m.par.e m.par.l m.par.h
    outs += "%g %g %g " % (Bpar.estimate[1],ptile(Bpar.mcestimates[:,1],2.5),ptile(Bpar.mcestimates[:,1],97.5)) # w.par.e w.par.l w.par.h
    outs += "%g %g %g " % (Bpar.deviance, ptile(Bpar.mcdeviance,95), psigcorrect.estimate_nu (Bpar)[0]) # d.par d.par.crit nu.par
    outs += "%g %g %g " % (Bpar.Rpd,ptile(Bpar.mcRpd,2.5),ptile(Bpar.mcRpd,97.5)) # rpd.par rpd.par.l rpd.par.h
    outs += "%g %g %g " % (Bpar.Rkd,ptile(Bpar.mcRkd,2.5),ptile(Bpar.mcRkd,97.5)) # rkd.par rkd.par.l rkd.par.h
    outs += ("%g "*options.nblocks) % tuple(Bpar.infl)
    outs += "%g " % runtime
    f.write(outs)

def write_bayes(f, mcmc, mcmc_conv, runtime):
    outs  = "%g %g %g %g %g " % (
            mcmc.estimate[0],        # m.bay.m
            mcmc.mapestimate[0],     # m.bay.map
            mcmc.posterior_median[0],# m.bay.median
            mcmc.getCI(1,(.025,)),   # m.bay.l
            mcmc.getCI(1,(.975)))    # m.bay.h

    outs += "%g %g %g %g %g " % (
            mcmc.estimate[1],                  # w.bay.m
            mcmc.mapestimate[1],               # w.bay.map
            mcmc.posterior_median[1],          # w.bay.median
            ptile(mcmc.mcestimates[:,1],2.5),  # w.bay.l
            ptile(mcmc.mcestimates[:,1],97.5)) # w.bay.h
    outs += "%g %g %g " % (mcmc.deviance,mcmc.bayesian_p('deviance'), psigcorrect.estimate_nu (mcmc)[0]) # d.bay d.bay.p
    outs += "%g %g " % (mcmc.Rpd,mcmc.bayesian_p('Rpd')) # d.rpd d.rpd.p
    outs += "%g %g " % (mcmc.Rkd,mcmc.bayesian_p('Rkd')) # d.rkd d.rkd.p
    outs += "%d %g %g %g " % (mcmc_conv,mcmc.Rhat(0),mcmc.Rhat(1),mcmc.Rhat(2))
    outs += ("%g "*options.nblocks) % tuple(mcmc.infl)
    outs += "%g " % runtime
    f.write(outs)

def write_data(f, data):
    # this conversion is required since data[0,:] does not work on lists
    data = asarray(data)
    outs  = ("%g "*options.nblocks) % tuple(data[:,0])
    outs += ("%d "*options.nblocks) % tuple(data[:,1].astype("i"))
    f.write ( outs )

#if options.datareduce:
#    outfile_reduced = open ( options.outputfile + "reduce","w" )
#    write_header( outfile_reduced )

count_npr = 0.
count_par = 0.
count_bay = 0.
not_converged = 0

random.shuffle ( x )

# write header
write_header(outfile)

sys.stderr.write("\n")
for simulation in xrange ( options.nsimulations ):
    sys.stderr.write ( "\nSimulation %d is running" % ( simulation, ) )
    O = create_new_observer ()
    # print "\nO=",O
    if not options.fixed_sequence:
        random.shuffle ( x )
    data = O.DoAnExperiment ( x, ntrials=options.blocksize )
    print "\ndata =",data
    print constraints

    write_id_gen(outfile, simulation)

    if nonparametric:
        Bnpr = pypsignifit.BootstrapInference ( data, priors=constraints, parametric=False, **ana_kwargs )
        if options.fixed_pmf:
            Bnpr.estimate = OO.params
        tic = time.time()
        Bnpr.sample ( options.nbootstrap )
        toc = time.time()
        count_npr += check_ci ( O, Bnpr )
        write_nonparametric(outfile, Bnpr, toc-tic)
        print "Done npar"

    if parametric:
        Bpar = pypsignifit.BootstrapInference ( data, priors=constraints, parametric=True,  **ana_kwargs )
        if options.fixed_pmf:
            Bpar.estimate = OO.params
        tic = time.time()
        Bpar.sample ( options.nbootstrap )
        toc = time.time()
        count_par += check_ci ( O, Bpar )
        write_parametric(outfile, Bpar, toc-tic)
        print "Done par"

    ####################
    # How to make sure that in the end ALL chains have converged?
    # We can give upper and lower limits for m and w from our sampling positions.
    # m cannot be outside the sampled range and w should not be wider than the sampled range (or twice that)
    if bayes:
        tic = time.time()
        mcmc = pypsignifit.BayesInference ( data, sample=True, priors=priors, **ana_kwargs )
        for prm in [0,1,2]:
            if not mcmc.geweke(prm)[2] is None:
                for j in mcmc.geweke(prm)[2]:
                    mcmc.resample(j)
        N = mcmc.mcestimates.shape[0]
        mcmc.sample( start = mcmc.farstart )
        mcmc.sample( start = mcmc.farstart )
        for prm in [0,1,2]:
            if not mcmc.geweke(prm)[2] is None:
                for j in mcmc.geweke(prm)[2]:
                    mcmc.resample(j)
        print "Rhat:  ",mcmc.Rhat(0),mcmc.Rhat(1),mcmc.Rhat(2)
        print "Geweke:",mcmc.geweke(0)[2],mcmc.geweke(1)[2],mcmc.geweke(2)[2]
        mcmc_conv = 1
        if mcmc.Rhat (0)>1.1 or mcmc.Rhat (1)>1.1 or mcmc.Rhat (2)>1.1:
            not_converged += 1
            mcmc_conv = 0
            # pypsignifit.ConvergenceMCMC(mcmc,0)
            # pypsignifit.ConvergenceMCMC(mcmc,1)
            # pypsignifit.ConvergenceMCMC(mcmc,2)
            # pypsignifit.GoodnessOfFit(mcmc)
            # pypsignifit.show()
            # sys.exit()
        # pypsignifit.ConvergenceMCMC(mcmc,1)
        else:
            count_bay += check_ci ( O, mcmc )

        # print count_bay, mcmc.estimate, pylab.prctile(mcmc.mcestimates[:,0], (2.5,97.5)), pylab.prctile(mcmc.mcestimates[:,1], (2.5,97.5))
        print count_bay, mcmc.getCI(1,(.025,.975))
        toc = time.time()
        write_bayes(outfile, mcmc, mcmc_conv, toc-tic)
        print "Done Bayes"

    write_data(outfile, data)
    outfile.write ("\n")
    outfile.flush()

# datareduce section is not up-to-date
# current state is that it doesn't seem to help much for improving the results
# in case you want to use this you must first make it work again
    if False:
        data = array ( data )
        dataml = data.copy ()
        nu = psigcorrect.estimate_nu (Bpar)[0]
        print "==============", nu, "==============="
        dataml[:,1] = ( data[:,1] * nu ).astype("i")
        dataml[:,2] = ( data[:,2] * nu ).astype("i")
        print dataml
        try:
            Bnpr = pypsignifit.BootstrapInference ( dataml, sample=options.nbootstrap, priors=constraints, parametric=False, **ana_kwargs )
        except:
            sys.stderr ( "An error ocurred in simulation %d during nonparametric bootstrap\n" % ( simulation,) )
        print "Done npar"
        try:
            Bpar = pypsignifit.BootstrapInference ( dataml, sample=options.nbootstrap, priors=constraints, parametric=True,  **ana_kwargs )
        except:
            sys.stderr ( "An error ocurred in simulation %d during parametric bootstrap\n" % ( simulation,) )
        print "Done par"

        ####################
        # How to make sure that in the end ALL chains have converged?
        # We can give upper and lower limits for m and w from our sampling positions.
        # m cannot be outside the sampled range and w should not be wider than the sampled range (or twice that)
        datamcmc = data.copy ()
        # nu = psigcorrect.estimate_nu (mcmc)[0]
        print "==============", nu, "==============="
        datamcmc[:,1] = ( data[:,1] * nu ).astype("i")
        datamcmc[:,2] = ( data[:,2] * nu ).astype("i")
        print datamcmc
        try:
            mcmc = pypsignifit.BayesInference ( datamcmc, sample=True, priors=priors, **ana_kwargs )
            for prm in [0,1,2]:
                if not mcmc.geweke(prm)[2] is None:
                    for j in mcmc.geweke(prm)[2]:
                        mcmc.resample(j)
            N = mcmc.mcestimates.shape[0]
            mcmc.sample( start = mcmc.farstart )
            mcmc.sample( start = mcmc.farstart )
            for prm in [0,1,2]:
                if not mcmc.geweke(prm)[2] is None:
                    for j in mcmc.geweke(prm)[2]:
                        mcmc.resample(j)
        except:
            sys.stderr ( "An error ocurred in simulation %d during mcmc sampling\n" % ( simulation,) )
        print "Rhat:  ",mcmc.Rhat(0),mcmc.Rhat(1),mcmc.Rhat(2)
        print "Geweke:",mcmc.geweke(0)[2],mcmc.geweke(1)[2],mcmc.geweke(2)[2]
        mcmc_conv = 1
        if mcmc.Rhat (0)>1.1 or mcmc.Rhat (1)>1.1 or mcmc.Rhat (2)>1.1:
            not_converged += 1
            mcmc_conv = 0
        writelog ( outfile_reduced, Bnpr, Bpar, mcmc, mcmc_conv )

sys.stderr.write ( "\r"+50*" "+"\n" )

outfile.close()

print "Coverages:"
if nonparametric:
    print "  nonparametric bootstrap:", count_npr/options.nsimulations
if nonparametric:
    print "  parametric bootstrap:   ", count_par/options.nsimulations
if bayes:
    print "  MCMC (bayesian):        ", count_bay/(options.nsimulations-not_converged)
    print "  MCMC runs that did not converge:",not_converged

# pypsignifit.show()

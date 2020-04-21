# wuhan
The data and R codes of COVID-19.

DATA

"contact-tracing.csv": all contact-tracing data.

"forwardtime.csv": the forward time (mixture) in our cohort.

"serial.csv": the serial intervals.

CODE

"contacttracing.R": to obtain serial intervals using contact-tracing data.

"incubation_Sx.R": simulation study for gamma, weibull and lognormal models respectively.

"incubation_data.R": to estimate the incubation period distribution, using gamma, weibull and lognormal models.

"incubation_CI.R": to obtain the confidence intervals of incubation period by bootstrap.

"incubation_loglik.R": to draw log-likelihood of incubation period vs pi, using gamma, weibull and lognormal models.

"density.R": to estimate the generation time distribution in simulation studies. The bandwidth should be chosen ad hoc.

"deconvolution.R": to estimate the generation time distribution using real data by boundary kernel deconvolution. The bandwidth should be chosen ad hoc.


Modification and reproduction of these files requires permission from the author.

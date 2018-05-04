## Overview of PCM 

O'Meara, Brain C. Evolutionary Inferences from Phylogenies: A Review of Methods. Annu. Rev. Ecol. Evol. Syst. 2012. 43:267–85

Review of the paper:
This gives a summary of all PCM methods. 
Continous-time markov chains with finite state space (decrete traits)
Figure 1 is a rate matrix which is really similar to corralated trait matrix

Figure 2 shows all different types of matric 
GTR nucleotide 
Binary correlation 
covarion -- that nucleotides substitution is different -- non-independence of sites; ITS sites are under this since its job is to fold protiens 
ordered transition

The only way to make covarion a "small" problem is to make it mixture model


Figure 3 -- Dealing with heterogeneity on the tree
F -- glaciler event, humming bird with red flower, C4 with C02 in environment 
C -- discrete gamma -- select rates from one "bin" within the gamma distribution; this is an intermediate number of GTR; other extremes are each site gets one or each site gets a new rate
	-> 8 is a good number of bins 
** each site only has one model 

Marginal reconstruction is more common than joint (the pie are the nodes)


Lambda transformation 
	tree + data + phyl model = Results
		Makes more sense with regression 
	transformed tree + data + reg. model = Results 
		Can work with one trait, but what if the two traits like different kappa (param)
		
		
		
	

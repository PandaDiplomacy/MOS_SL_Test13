def retrieve_LF(z,line):
   if line == 'H_alpha':
      logLStar = 0.45*z + 41.87
      logPhiStar = -0.38*z + z - 3.18
      alpha = -1.60
   elif line == 'O_II': # fitted manually from Ly et al. 2007 http://arxiv.org/pdf/astro-ph/0610846v1.pdf
      logLStar = 0.54*z + 41.56
      logPhiStar = 0.37*z - 2.72
      alpha = 0.35*z -1.46 #alpha = 0.61*z - 1.81 (Fitted in excel)
   return alpha,logLStar,logPhiStar

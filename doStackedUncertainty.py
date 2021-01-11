import numpy as np
import gecorg as go
import glob


if __name__=='__main__':
    #will need to add the stuff for scale factor calc and so on, that will need to be propagated
    zptcut = 200.0
    hptcut = 250.0
    metcut = 250.0
    bkgerrfs = go.gatherBkg('analysis_output_ZpAnomalon/2021-01-11','selected_errors',zptcut,hptcut,metcut)
    print bkgerrfs

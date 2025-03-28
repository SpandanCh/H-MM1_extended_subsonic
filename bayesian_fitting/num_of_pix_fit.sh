#!/bin/bash

ek=0
for fldr in nested-sampling/chains/*

    do fle=$fldr"/3-summary.txt" 

        if [ -f "$fle" ]
        then ek=$(($ek+1))
        fi

    done
    
echo "number of pixels with three-component fit : "$ek

# total number of pixels with snr above 2 (cut given in bayesian fit) : from python
tot_pix_snr=13634 #11547 #13634 #10519 
perc1=$(($ek*100 / $tot_pix_snr))
echo "% complete : "$perc1

# difference
diff=$(($tot_pix_snr - $ek))
echo "number of pixels left :"$diff

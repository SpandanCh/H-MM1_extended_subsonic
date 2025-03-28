#!/bin/bash

for fldr in nested-sampling/chains/*
    do
        fle=$fldr"/3-.txt"
        fle_sze=$(stat -c%s "$fle")
 
        if  (( fle_sze<1000 ))
            then
            rm $fldr"/3-"*
            declare -a l   
            for ((i=0; i<${#fldr}; i++))
                do 
                    l[$i]="${fldr:$i:1}"
                done
            x=${l[-5]}${l[-4]}
            y=${l[-2]}${l[-1]}
            python innocent_script.py 1 $y $x
            echo $fldr
        fi
        
    done

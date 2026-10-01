#!/usr/bin/env bash
for run in 6077; do
  root -l -b -q "Plot_CDet_GoodPulseCandidates_AllTDC.C+(${run},\"${OUT_DIR}\",\"CDet_run${run}_good_pulse_tdc\",1.0,0.0,60.0,0.0,40.0,false,0.0,-26.0,0.020,5.0,2.0,0.0,-26.0,0.040,5.0,2.0,-10.0,10.0,3.0,4.5,true,0.51,0.08,0.0,0.17)"
done

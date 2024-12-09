 
copy initvalues_postjitter_%1.ctl initvalues.ctl
copy pop24_2_fixM.ctl pop24.ctl
pop24 -nox 

move pop24.rep .\results\pop24_%1.rep
move pop24.std .\results\pop24_%1.std
move pop24.par .\results\pop24_%1.par
move pop24.rdat .\results\pop24_%1.rdat
move pop24.cor  .\results\pop24_%1.cor

if NOT exist .\results mkdir .\results   

copy pop24_%1.dat pop24.dat
pop24 -nox 

move pop24.rep .\results\m_24_2_%1.rep
move pop24.std .\results\m_24_2_%1.std
move pop24.par .\results\m_24_2_%1.par
move pop24.rdat .\results\m_24_2_%1.rdat
move pop24.cor  .\results\m_24_2_%1.cor

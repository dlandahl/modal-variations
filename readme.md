
Work in progress implementation of a research paper for generating random procedural variations of impact sound effects by analysing resonant modes, randomising their amplitude, extracting the residue (signal without the most prominent modes), and then creating random variations of that using dip filters.
<br/>One key feature of this research is to extract an arbitrary time-amplitude envelope to allow resynthesising the modes more accurately.
<br/>You might expect sound energy in a system to diffuse exponentially (rate of energy leaving the system is proportional to the amount of energy in the system at that time), but this idealisation is not quite right due to natural phenomena such as energy transfer between modes and other non-linear effects.

<br/><br/> Right now we are doing modal analysis, but not extracting the residual or resynthesising the modes. I will also implement the random spectral dips using biquadratic filters, which will be useful for sounds which do not have strong resonances.
<br/><br/>[http://www.nikunjr.com/Projects/Crackdown/crackdown.pdf](Link) to the paper

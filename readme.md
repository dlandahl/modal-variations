
Work in progress implementation of a research paper for generating procedural variations of impact sound effects by analysing the resonant frequencies, extracting the residue (signal without the most prominent modes), and then creating random variations by resynthesising the modes with randomised amplitudes each time the sound is triggered. This can also save a lot of system memory, because the most prominent modes of resonance take up all of the tail of the sound in some classes of sounds, and only the residue has to be stored as audio data.

One key feature of this research is to extract an arbitrary time-amplitude envelope to allow resynthesising the modes more accurately. You might expect sound energy in a resonant mode to diffuse exponentially (rate of energy leaving the system is proportional to the amount of energy in the system at that time), but this idealisation is not quite right due to natural phenomena such as energy transfer between modes and other non-linear effects. Extracting, storing, and applying an arbitrary amplitude envelope is very cheap (a float per STFT frame, lerp between them), and makes the sound much more realistic.

Right now we are doing modal analysis, but not extracting the residual or resynthesising the modes. I have also implemented the random spectral dips using biquadratic filters, which will be useful for sounds which do not have strong resonances.
<br/><br/>[Link](http://www.nikunjr.com/Projects/Crackdown/crackdown.pdf) to the paper

# EccentricIMR

## About

A Mathematica package for generating gravitational waveforms for
nonspinning eccentric binary black hole mergers.

## Installation

### Using Git (recommended)

Change into your Mathematica applications directory

    cd ~/Library/Mathematica/Applications

Clone the repository

    git clone https://github.com/ianhinder/EccentricIMR.git

### Zip file download

Download
[master.zip](https://github.com/ianhinder/EccentricIMR/archive/master.zip),
extract it, rename the extracted directory EccentricIMR-master as
EccentricIMR, and move it into your Mathematica applications directory
(~/Library/Mathematica/Applications for Mac OS,
~/.Mathematica/Applications for Linux).

## Quick start

Open a new Mathematica notebook and enter the following:

    << EccentricIMR`;

    params = <|"q" -> 1, "x0" -> 0.07, "e0" -> 0.1,
               "l0" -> 0, "phi0" -> 0, "t0" -> 0|>;

    hEcc = EccentricIMRWaveform[params, {0, 10000}];

    ListLinePlot[Re[hEcc]]

![Eccentric waveform](ecc-waveform.png)

## Documentation

### EccentricIMRWaveform[_parameters_, {_t1_, _t2_}]

Generate an eccentric inspiral-merger-ringdown waveform with the parameters given in the time range {t1,t2}.  

_parameters_ is an Association with the following entries:

Parameter | Meaning
--------- | ---
q         | Mass ratio of the binary (q=m1/m2)
t0        | _Reference time_ at which the remaining parameters are quoted
x0        | Dimensionless frequency parameter (x = (M om_orb)^(2/3)) evaluated at the reference time
e0		   | Eccentricity at the reference time
l0		   | Mean anomaly at the reference time
phi0	   | Orbital phase at the reference time

See [arXiv:0806.1037](http://arxiv.org/abs/arXiv:0806.1037) for full details about the meaning of the parameters.

The returned waveform is expressed as a list of {t, h22} pairs, where t is the retarted time coordinate and h22 is the l=2, m=2 spin-weighted spherical harmonic coefficient of the waveform.

Example:

    In[1]:= hEcc = EccentricIMRWaveform[<|"q" -> 1, "x0" -> 0.07,
      "e0" -> 0.1, "l0" -> 0, "phi0" -> 0, "t0" -> 0|>,
      {0, 10000}];

    In[2]:= Take[hEcc, 10]

    Out[2]:= {{0., -0.127741 + 0.000381028 I}, {1., -0.127596 + 
    0.00619006 I}, {2., -0.127182 + 0.0119862 I}, {3., -0.126498 + 
    0.0177567 I}, {4., -0.125546 + 0.0234884 I}, {5., -0.12433 + 
    0.0291683 I}, {6., -0.122851 + 0.0347835 I}, {7., -0.121113 + 
    0.0403215 I}, {8., -0.119122 + 0.0457697 I}, {9., -0.11688 + 
    0.051116 I}}




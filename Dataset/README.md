# DATA SET FOR FLUIDIC PINBALL AT RE=80

Nan Deng

Last edited: Oct 22, 2021

Data link: [Dropbox link for dataset of Fluidic Pinball at Re80](https://www.dropbox.com/sh/eqsmi8o6zctsp5y/AAArortRZxmu-0gbXNtOdtqta?dl=0)

---
## Fluidic Pinball -- a sandbox to test different flow modelling and control strategies

In this folder is the data set of **transient and post-transient flow dynamics** associated with 
the **two successive supercritical bifurcations** of **Hopf** and **pitchfork** type for the fluidic pinball at Re=80.

### Data generator:

The instantaneous flow flild is calculated by **direct numerical simulation** (DNS) of the Navier-Stokes equations 
for the two-dimensional incompressible wake flow around the fluidic pinball.
The unsteady Navier-Stokes solver (UNS3 solver by Marek Morzyński, https://gitlab.com/PinBall) is 
based on a **second-order finite-element discretization method** of the Taylor-Hood type, 
on an **unstructured grid** of **4 225 triangular elements** and **8 633 vertices**, and an **implicit** integration of the **third-order in time**.

**Numerical details** to see:
    Fluidic Pinball - A Toolkit for Exploring Multiple-Input Multiple-Output Flow Control
    http://berndnoack.com/FlowControl.php

**Least-order mean-field model** to see:
    Low-order model for successive bifurcations of the fluidic pinball.
    Deng, N., Noack, B. R., Morzyński, M. and Pastur, L. R. 
    Journal of Fluid Mechanics, 884, A37.
    https://doi.org/10.1017/jfm.2019.959

**Transient and post-transient dynamics** and **Galerkin force model** to see:
    Galerkin force model for transient and post-transient dynamics of the fluidic pinball.
    Deng, N., Noack, B. R., Morzyński, M. and Pastur, L. R. 
    Journal of Fluid Mechanics, 918, A4.
    https://ddoi:10.1017/jfm.2021.299

---
## Flow fields (u,v,p) at Re=80
    
### Mesh Info:

1. Visulalization.m -- Plot the snapshots

2. Mesh/Grid2.dat -- x and y coordinates of the 8633 nodes.
    
3. elem.dat -- Connectivity matrix of the 4225 TIRA6 elements

4. Comp_Vorticity.m -- Function for calculating the vorticity field from the velocity field and plotting the vortex shedding.

5. redblueTecplot.m -- Colorbar from blue to red

6. find.dat -- A matrix for mirror-refelecting data

### Three steady solutions:

SS_SYM, SS_ASYM_UP, SS_ASYM_DOWN.

### Three transitions from different steady solutions are provided:
1. Dataset starting with symmetric** steady solution, T = 0.0 : 0.1 : 1500.

    UALL_tran_Re80_sym.mat
    
    VALL_tran_Re80_sym.mat
    
    PALL_tran_Re80_sym.mat
    
    15001 snapshots of u and v velocity as well as pressure (UALL, VALL, PALL).

2. Dataset starting with the **upwards** deflected asymmetric steady solution, T = 0.0 : 0.1 : 1000.

    UALL_tran_ASYM80_UP.mat
    
    VALL_tran_ASYM80_UP.mat
    
    PALL_tran_ASYM80_UP.mat
    
    10001 snapshots of u and v velocity as well as pressure (UALL, VALL, PALL).

3. Dataset starting with the **downwards** deflected asymmetric steady solution, T = 0.0 : 0.1 : 1000.

    UALL_tran_ASYM80_DOWN.mat
    
    VALL_tran_ASYM80_DOWN.mat
    
    PALL_tran_ASYM80_DOWN.mat
    
    10001 snapshots of u and v velocity as well as pressure (UALL, VALL, PALL). 
    
---
## References
[1] Noack, B. R. & Morzynski, M. 2017 The fluidic pinball — a toolkit for multiple-input multiple-output flow control (version 1.0). Tech. Rep. 02/2017. Chair of Virtual Engineering, Poznan University of Technology, Poland.

[2] Deng, N., Noack, B. R., Morzyski, M. & Pastur, L. R. 2020 Low-order model for successive bifurcations of the fluidic pinball. Journal of Fluid Mechanics 884, A37.

---
## Contacts:
For further information please contact:

Email adress: 

[nan.deng@ensta-paris.fr](nan.deng@ensta-paris.fr)

[j.nan.deng@gmail.com](j.nan.deng@gmail.com)

Homepage: [https://sites.google.com/view/nandeng](https://sites.google.com/view/nandeng)



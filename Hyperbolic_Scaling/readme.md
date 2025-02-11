# **Hyperbolic_Scaling Directory**  

This directory contains the MATLAB scripts and functions used for simulating the **2D glioma model** on a square grid. The results correspond to **subsection 2.5.5** of the thesis *Diss_Kumar_Pawan.pdf* (available in the parent directory).  

## **Main Files/Scripts**  

- **`main_hyp_all.m`** – The main MATLAB script that calls all necessary functions and simulates the **2D tumor model** on a square grid. It generates the majority of the results presented in **subsection 2.5.5** of the thesis.  

## **Functions**  

- **`tissue_Q_macro.m`** – Computes the **macroscopic tissue Q** described in **subsection 2.5.2** of the thesis.  
- **`function_g.m`** – Computes **g(S(x))**, which appears in **Equation 2.4.46**.  
- **`tissue_unimodal.m`** – Computes the **directed tissue q**, given by **Equation 2.5.5**.  
- **`set_tumor_diff.m`** – Assembles the **tumor diffusion matrix** and stores related variables.  
- **`set_acidity_diff.m`** – Assembles the **acid diffusion matrix**.  

## **Scripts**  

Each script simulates the model under different conditions and produces results corresponding to specific figures in the thesis.  

- **`main_SU_delta1.m`** – Simulates the model (**Equations 2.4.46 & 2.3.7**) with **anisotropic tissue**, **two occlusion sites**, and **delta = 1** in **Equation 2.5.5**. (*Figures 2.12 & 2.13*)  
- **`main_SU_exp2_delta0p2.m`** – Simulates the model with **anisotropic tissue**, **two occlusion sites**, and **delta = 0.2** in **Equation 2.5.5**. (*Figures 2.14 & 2.15*)  

## **Numerical Method**  

### **Time Discretization**  
- **IMEX Method** (Implicit-Explicit scheme):  
  - **Diffusion terms** are solved **implicitly** (Implicit Euler).  
  - **Source and taxis/advection terms** are solved **explicitly** (Euler).  

### **Space Discretization**  
- **Weickert’s 3×3 stencil discretization** *(9-point central difference, as described in Table 2.3 of the thesis)* is used for **tumor diffusion**, based on Weickert, J. *Anisotropic Diffusion in Image Processing*. Teubner Stuttgart, 1998. [Available here](http://www.lpi.tel.uva.es/muitic/pim/docus/anisotropic_diffusion.pdf).  
- **Standard 5-point stencil** *(central difference, as shown in Table 2.2 of the thesis)* is used for **acid diffusion**.  
- **Second-order upwind scheme with Van Leer flux limiter** is used for the **taxis term** in both x and y directions (*subsection 2.5.3 of the thesis*).  

## **Performance**  

- **Simulation Time:**  
  - **30-35 minutes** on an **8-core, 16GB RAM, i7-processor @ 3.60GHz (Linux)**.  
  - **35-40 minutes** on a **dual-core, 8GB RAM, i5-processor @ 2.3GHz (Mac)**.  

## **Requirements**  

- **MATLAB Version:** 2018a or later (older versions may also work).  


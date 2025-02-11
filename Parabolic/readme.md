
# **Parabolic Directory**  

This directory contains the MATLAB scripts and functions used for simulating the **2D glioma model** on a square grid. The results correspond to **subsection 2.5.4** of the thesis *Diss_Kumar_Pawan.pdf* (available in the parent directory).  

Additionally, the anisotropic diffusion method used in this work is based on the approach presented in:  

**Weickert, J.** *Anisotropic Diffusion in Image Processing.* Teubner Stuttgart, 1998. [Available here](http://www.lpi.tel.uva.es/muitic/pim/docus/anisotropic_diffusion.pdf).  

## **Main Files/Scripts**  

- **`main_all.m`** – The main MATLAB script that calls all necessary functions and simulates the **2D tumor model** on a square grid. It generates the majority of the results presented in **subsection 2.5.4** of the thesis.  

## **Functions**  

- **`tissue_Q_macro.m`** – Computes the **macroscopic tissue Q** described in **subsection 2.5.2** of the thesis.  
- **`function_g.m`** – Computes **g(S(x))**, which appears in **Equation 2.4.16**.  
- **`tissue_q_un.m`** – Computes the **undirected tissue small q**, given by **Equation 2.5.4**.  
- **`set_tumor_diff.m`** – Assembles the **tumor diffusion matrix** and stores related variables.  
- **`set_acidity_diff.m`** – Assembles the **acid diffusion matrix**.  

## **Scripts**  

Each script simulates the model under different conditions and produces results corresponding to specific figures in the thesis.  

- **`main_2D_aniso_I.m`** – Simulates the model (**Equations 2.4.16 & 2.3.7**) with **anisotropic tissue** and **three occlusion sites** as initial conditions. (*Figure 2.7*)  
- **`main_2D_aniso_II.m`** – Simulates the model with **anisotropic tissue** and **two occlusion sites** as initial conditions. (*Figure 2.8*)  
- **`main_2D_iso_I.m`** – Simulates the model with **isotropic tissue** and **three occlusion sites** as initial conditions. (*Figure 2.5*)  
- **`main_2D_iso_II.m`** – Simulates the model with **isotropic tissue** and **two occlusion sites** as initial conditions. (*Figure 2.6*)  
- **`main_2D_lower_grades.m`** – Simulates the model with **anisotropic tissue**, **two occlusion sites**, and **stronger proton buffering**. (*Figure 2.9*)  
- **`main_2D_modified_source.m`** – Simulates the model with **anisotropic tissue**, **two occlusion sites**, and a **modified source term** for the glioma cell equation. (*Figure 2.10*)  

## **Numerical Method**  

### **Time Discretization**  
- **IMEX Method** (Implicit-Explicit scheme):  
  - **Diffusion terms** are solved **implicitly** (Implicit Euler).  
  - **Source and taxis/advection terms** are solved **explicitly** (Euler).  

### **Space Discretization**  
- **A 3×3 stencil discretization**, based on **Weickert, J.** *Anisotropic Diffusion in Image Processing*, Teubner Stuttgart, 1998. [Available here](http://www.lpi.tel.uva.es/muitic/pim/docus/anisotropic_diffusion.pdf), is used for **tumor diffusion** (*9-point central difference, as described in Table 2.3 of the thesis*).  
- **A standard 5-point stencil** (central difference) is used for **acid diffusion** (*Table 2.2 in the thesis*).  
- **An upwind scheme** (first-order) is applied for the **taxis term** in both the x and y directions.  

---

This way, the reference blends naturally into the explanation without breaking the flow. Let me know if you'd like any more refinements! 🚀

## **Performance**  

- **Simulation Time:**  
  - **4-5 minutes** on an **8-core, 16GB RAM, i7-processor @ 3.60GHz (Linux)**.  
  - **5-6 minutes** on a **dual-core, 8GB RAM, i5-processor @ 2.3GHz (Mac)**.  

## **Requirements**  

- **MATLAB Version:** 2018a or later (older versions may also work).  

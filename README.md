# **Repository Overview - PseudoPalisadesMultiScale**

Welcome to this GitHub repository, which houses the code associated with the article **"Multiscale modeling of glioma pseudopalisades: contributions from the tumor microenvironment"**, available at [Springer](https://link.springer.com/article/10.1007/s00285-021-01599-x) and in **Chapter 2** of the thesis at [KLUEDO](https://kluedo.ub.rptu.de/frontdoor/index/index/docId/6573).  

This repository primarily focuses on the **thesis results**, which provide a more detailed and comprehensive presentation of the work compared to the article. The thesis includes everything covered in the paper, along with additional insights and extended analyses.  

## Key Features  

This repository provides implementations for:  

- **Numerical simulations** of tumor dynamics  
- **Data analysis** related to glioma modeling  
- **Visualization of results** from the simulations  

## Project Structure  

### **Parabolic**  
Contains codes and results corresponding to **Figures 2.3 – 2.10** in the thesis.  

### **Difference**  
Contains codes and results for **Figure 2.11** in the thesis.  

### **Hyperbolic**  
Contains codes and results for **Figures 2.12 – 2.15** in the thesis.  

### **Pattern**  
Contains codes and results for **Figure 2.16** in the thesis.  

### **TWA**  
Includes codes used for **Appendix 2.B** of the thesis.  

### **MS_Phase_plane.nb**  
A **Mathematica notebook** for plotting the **MS phase portrait** included in **Figure 2.17** of the thesis. The corresponding plot is available as **MS_phase_plot.eps**.  

### **Others**  
This subdirectory contains implementations of the **Discontinuous Galerkin Method (DGM)** for solving the model on **real brain geometry**. The DGM implementation is written in **C++** and utilizes the [deal.II library](https://www.dealii.org).  

DGM is applied to real brain geometry by integrating imaging data—particularly **Diffusion Tensor Imaging (DTI)**—to facilitate analysis and incorporate it into the model.  

To run this part, **data.ncdf** is required in the `data` folder. Alternatively, synthetic DTI data can be used.  

---
Each subdirectory contains a separate README file explaining its contents in detail.

Feel free to explore the directories to dive deeper into specific aspects of the project. If you have any questions or need further clarification, don't hesitate to reach out.  

## **Author**  
**Pawan Kumar**: [@its-Pa1](https://github.com/its-Pa1)  

© 2025, Pawan Kumar. All Rights Reserved.  


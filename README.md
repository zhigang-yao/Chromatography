## MGDG and MALG sampler
This is a MATLAB implementation of our sampling method described in our paper *A Statistical Approach to Estimating Adsorption-Isotherm Parameters in Gradient-Elution Preparative Liquid Chromatography*, [arXiv:2201.00958](	arXiv:2201.00958).

#### Description
There are several functions
- **solver.m**: the numerical solver one get $\mathbf{Y}$ with $\mathbf{\xi}$ and time $\mathcal{T_n}$;
- **compress.m**: function getting $\eta$ from $\xi$;
- **recover.m**: function getting $\xi$ from $(\nu, \eta)$;
- **Grad.m**: function calculating the numerical gradient;
- **GD2D.m**: function running the gradient descent algorithm in 2D case;
- **MALG.m**: main function running our MALG algorithm;
- **MGDG.m**: main function running our MGDG algorithm.

To get most of the result in our paper, you can run the scripts provided by running the following code in matlab:

    Gen_data_get_init;
    Comp_MH_MGDG_MALG;
    Run_MALG;
    Run_MGDG;
    Draw_Summary_plots;

Note that the above code calls the *parfor* function and can take up a lot of computation time and resources.

#### Customization

If you want to apply our samplers to other problems, we recommend that you modify the simple functions in *solver.m*, *compress.m*, *recover.m*, *Grad.m*, *GD2D.m* and the sampling dimension in *MALG.m* and *MGDG.m* accordingly, and then tunning the proposal distributions as we discussed in our paper.

#### Citation
Please cite our paper with the following bib information if you use this code in your own work:

    @article{su2023statistical,
    title={A statistical approach to estimating adsorption-isotherm parameters in gradient-elution preparative liquid chromatography},
    author={Su, Jiaji and Yao, Zhigang and Li, Cheng and Zhang, Ye},
    journal={The Annals of Applied Statistics},
    volume={17},
    number={4},
    pages={3476--3499},
    year={2023},
    publisher={Institute of Mathematical Statistics}
    }
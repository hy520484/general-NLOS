# Demo code

This repository contains the code for the paper *Non-line-of-sight imaging with arbitrary illumination and detection pattern*.

# System Requirements

## Hardware Requirements

+ We recommend a parallel computing platform with enough RAM to support the in-memory operations. The code has been tested on an AMD EPYC 7452 server with the following settings.

  > RAM: 256 GB
  >
  > CPU: 64 cores, 2.35 GHz/core

## Software Requirements

### OS Requirements

+ The code has been tested on a server with the following setting.

  > CentOS 7
  >
  
+ The code runs with Matlab R2020a or higher versions.

### Installation

+  The code only depends on the Matlab software and does not need additional installation.

# Instructions

### CC-SOCR 

+ **How to run the code** To reproduce the reconstruction results in the article, choose `id` from 1 to 25, and run `main.m`. 

  > Each `id` is uniquely related to the `Instance_name` and the `Relay_setting`.

+ **Expected output** A directory named *results*\_<font color = blue>*Instance_name*</font>\_<font color = blue>*Relay_setting*</font> will be created. The following files will be saved in the folder.

  > `d1.mat`, `d2.mat`,... : The updated virtual confocal signal in each iteration. The number of  virtual confocal signals equals `num_main_loop-1`.
  >
  > `u1.mat`, `u2.mat`,... : The updated directional albedo in each iteration. The number of the reconstructions equals `num_main_loop`.
  >
  > `albedo.png` : The reconstructed albedo value viewed from the relay surface.
  >
  > `three view.png` : The three views of the reconstructed albedo.
  >
  > `x-component`, `y-component`, `z-component`: The three components of the directional albedo in the depth, horizontal, and vertical directions. These results indicate the surface normal of the reconstructed targets.
  >
  > `illumination_pattern.png` : Coordinates of the illumination points.

  > **An example** By default, `id = 1`. This is the case with `Instance_name = 'bunny'` and `Relay_setting = 'full'`. After running `main.m`, the directory *results_bunny_full* is created. For this instance, we set `num_main_loop = 2`. The files in the directory are
  >
  > `d1.mat`, `u1.mat`, `u2.mat`, `albedo.png`, `three_view.png`, `x-component.png`, `y-component.png`, `z-component.png`, and `illumination_pattern.png`. Some of these results are shown in Supplementary Figure 1.

+ **Run time** We recommend to run the code on parallel computing platforms. The execution time for the instance of the statue with signals measured at 200 randomly distributed focal points （`id=22`） is 281 s for 2 full iterations with an AMD EPYC 7452 server with 64 CPU cores. See *Supplementary Table 3* for more details.

+ **Run the code on new datasets**  To run the code on new datasets, the following data should be provided

  > `C_i.mat` The coordinates of the illumination points
  >
  > `C_d.mat` The coordinates of the detection points
  >
  > `sub_Sig.mat` The input signal
  >
  > `experimental_setup` The experimental setup including *temporal resolution*, *reconstruction domain*, *voxel size*, and *virtual focal points*.
  >
  > `parameters` The parameters used for the reconstruction. 
  >
  > `save_results_dir` The name of the directory created for saving results.

### PF-BP

For the phasor field method with back-projection implementation [1], we provide the implemented code.

Run  `main_phasor_BP.m` in the folder `phasor_BP_example`. The output is shown in the third row and second column of Figure 7. 

# License

This project is covered under the **Apache 2.0 License**.



## Reference

1. Liu, X. *et al.* Non-line-of-sight imaging using phasor-field virtual wave optics. *Nature* **572**, 620–623 (2019).
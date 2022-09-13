# Test computation using GUI

Provided are brief manuals to compute a functional of the geopotential (here,
the disturbing potential) and its commission error.  All paths in this section
are relative to the top directory of the repository.

## Functional of the geopotential

* Panel `Geopotential model and reference system selection`

  * Import the gravity field model file `./data/input/EGM96.mat` using the 
    `Browse...` button.


* Panel `Point type selection`

  * Select the checkbox `Load data` and click the `Browse...` button (next to
    the checkbox) to import the computation points from the
    `./data/input/sctr-points.txt` file.


* Panel `Calculated parameters and output selection`

  * From the first pop-up menu, choose `Disturbing_potential`.

  * Enter name of your output file using the `Output folder and file` button
    (enter only the name of the file, i.e. without any suffix).

  * Click the `OK` button to start the computation.


* Compare your results with the ones from the enclosed sample output file
  `./data/output/output_EGM96.txt`.


## Commission error

All paths in this section are relative to the top directory of the repository.

* Panel `Geopotential model and reference system selection`

  * Import the error variance--covariance matrix file 
    `./data/input/GRIM5C1_covmat.mat` using the `Browse...` button.

  * Set `nmin` to `2`.


* Panel `Point type selection`

  * Select the checkbox `Load data` and click the `Browse...` button (next to
    the checkbox) to import the computation points from the
    `./data/input/sctr-points.txt` file.

* Panel `Calculated parameters and output selection`

  * Select the checkbox `Commission error`.

  * From the first pop-up menu, choose `Disturbing_potential`.

  * Enter name of your output file using the `Output folder and file` button
    (enter only the name of the file, i.e. without any suffix).

  * Click the `OK` button to start the computation.

*  Compare your results with the ones from the enclosed sample output file
   `./data/output/output_GRIM5C1.txt`.






# Test computation using the command line

All paths in this section are relative to the `src` directory of the
repository.

* In MATLAB, change your current working directory to `src`.

* Execute the following command:

  ```
  GrafLab('OK', ...
          3986004.415E+8, ...
          6378136.3, ...
          0, ...
          360, ...
          1, ...
          '../data/input/EGM96.mat', ...
          0, ...
          1, ...
          [], ...
          [], ...
          [], ...
          [], ...
          [], ...
          [], ...
          [], ...
          '../data/input/sctr-points.txt', ...
          [], ...
          [], ...
          [], ...
          '../data/output/my_output_EGM96', ...
          0, ...
          5, ...
          1, ...
          [], ...
          1, ...
          1, ...
          0, ...
          [], ...
          [], ...
          [], ...
          [], ...
          [], ...
          1)
  ```

  The description of the input variables can be found in the source code of
  GrafLab (look for the text: `function GrafLab` in `./GrafLab.m`) or in the
  GrafLab [cookbook](https://github.com/blazej-bucha/graflab-cookbook).

* Compare your output file, `../data/output/my_output_EGM96.txt`, with the
  reference file `../data/output/output_EGM96.txt`.


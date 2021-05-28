[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SebastianHaben/icc-pyir/master)
# icc-pyir
This project's aim is to analyse Pyridine IR data by caluculating the BAS and LAS density of the material. icc-pyr takes a pre- and post-adsorption spectrum as input and can create a PDF report.   

## Installation
Prerequisites:

* Python >=3.7

Dependencies:
The code was tested under Win10 using the following dependencies with their respective versions.

* jinja2 = 2.10
* lmfit = 0.9.15
* matplotlib = 3.1.1
* numpy = 1.17.2
* scipy = 1.2.0
* xhtml2pdf = 0.2.5

To install using git:

`git clone https://github.com/SebastianHaben/icc-pyir.git`

## Usage
icc-pyir is centered around the `Analyser` class. Firstly, load the import it from the `analyser` module and create an `pyr_analyser` object:
```
from analyser import Analyser
pyr_analyser = Analyser()
```

Afterwards, we add a pre-adsorption and postadrosption spectrum.
```
    pyr_analyser.add_pre_adsorption_spectrum("..\\test\\data\\SH-00.002.00\\SH-00.002.00_pre_adsorption.SPA")
    pyr_analyser.add_post_adsorption_spectrum("..\\test\\data\\SH-00.002.00\\SH-00.002.00_post_adsorption.SPA")
```

If you want you can add a normailzation step by using the `normalize_spectra` function and providing a band guess. This step is optional.
Next we need we create the difference spectrum and find the pyridine bands and a search window where they should roughly be. Since the pyridine bands are usually located at the same position the respective band values were choosen as the default values. 

    pyr_analyser.calculate_difference_spectra()
    pyr_analyser.find_py_adsorption_peaks()

After finding the peaks we will fit them using Gaussian functions and also integrating them.

    pyr_analyser.fit_py_adsorption_peaks()
    pyr_analyser.integrate_fitted_peaks()

Afterwards we can calculate the acid site density by providing the sample mass (mg) and the sample disc radius (cm).

    pyr_analyser.calculate_acid_sites(17.38)

Finally we can create our PDF-report by creating an PDFReport object and parsing the analyser object as well as the your Material ID:

    from reports import PDFReport
    report = PDFReport(pyr_analyser, 'SH-00.002.00')
    report.build_report()
    report.save_report()

## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)

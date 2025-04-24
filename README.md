# GCProcessing
## Scope
This program is designed to process gas chromatography data. It is capable of reading the data from a CSV file, 
performing calibration, and calculating the amount of compounds in the sample. Although the program is designed with 
two-dimensional chromatography in mind, it can also be used for one-dimensional chromatography despite the different terminology.
For the optimal performance of the code, the data should be in a specific format. Failure to follow the format may result in 
inconsistencies, mistakes, malperformance, or even crashes. The program is designed to be as user-friendly as possible without 
hiding behind a sophisticated interface. It is therefore recommended to read the documentation and the code carefully before 
using the program. The code generally contains the main functions necessary to process GC data, but the main objective of this
program is to provide a framework for the user to build upon. The code is not intended to be a complete solution, but rather a starting 
point for user-specific developments. 

Please contact me if you have any questions or suggestions. I am open to any feedback and will be happy to help you with the code.

## Technical Requirements
The code is written in Python 3.12 and requires the following libraries:

* pandas
* numpy
* matplotlib
* xlsxwriter
* rdkit
* openpyxl
* scikit-learn
* tqdm
* pywin32
* chemicals

Make sure to install the latest version of all the libraries before running the code. No specific functionality is used from the 
latest versions of python or the libraries to creat incompatibilities with older versions, but it is recommended to use the latest.

## Input data
The code requires the following input data all in comma-delimited csv format:
### Blob table
This is where the data of the sample (after peak detection and naming) is stored. The format is the default format provided by
GC Image software. The table should contain the following columns:
* _Compound Name_: The name of the compound (blob). Multiple blobs are allowed to have the same name. Their volumes will be summed up and
their retention times is averaged.
* _Retention Time I (min)_: The retention time of the compound in the first dimension. This is used to sort the blobs.
* _Retention Time II (sec)_: The retention time of the compound in the second dimension. This is used as the secondary dimension for sorting.
* _Volume_: The volume of the blob (or area of the peak) in the chromatogram.
* _Inclusion_: A boolean value (i.e., TRUE/FALSE) indicating whether the blob/peak is included in the quantification or not.
* _Amount_: The amount of the compound in the sample. This is **ONLY** used for internal standard blobs as this is the program's
method for determining whether a blob is an internal standard or not. For other compounds, this column should be empty or 0.
* _Internal Standard_: The internal standard used for the compound. This is an integer assigned by GC Image software to indicate which 
blob is the internal standard used for the quantification of the specific compound. The actual values can be arbitrary, as long as all the
compounds with the same internal standard (including the internal standard peak itself) have the same value.

The correct spelling of each column heading is important, but the order is not. Ensure to provide all the columns for optimal performance.

### Database
This is a database file containing the information of the compounds such as their formula, grouping, molecular weight, etc.
An example (_db.csv_) is provided with the code, but the user can create their own database file. The database should contain all
the columns included in the example file. 

### NIST library
This is a library file containing the information of the compounds such as their formula, molecular weight, InChI, etc. Use the 
file _nist_compounds.csv_ provided with the code.

### Calibration data
This is only for external calibration data. The csv file should contain only two columns. The first column must contain the amount
(wt.% for liquid and bar for gases used in the µ-Pyrolyzer unit) and the second column the peak volume. The headings are not important,
but the order is.
## Formulas

To calculate the compound amount, use the following formula:

self.mol = (((self.volume - calibrant.curve[1]) / calibrant.curve[0]) * calibrant.compound.mrf) / (
                calibrant.compound.mol_wt * 1e6 * self.compound.mrf)

$n_{compound} = \frac{(V_{compound} - b_{cal})}{a_{cal}} \times \frac{MRF_{cal}}{MRF_{compound} * MW_{compound} * 10^6}$

Where,
$n_{compound}$ is the compound amount in $mol$,
$V_{compound}$ is the peak volume of the compound,
$b_{cal}$ is the y-intercept of the calibration curve ($b=0$ for internal calibration),
$a_{cal}$ is the slope of the internal calibration curve,
$MRF_{cal}$ is the molecular response factor of the calibrant,
$MRF_{compound}$ is the molecular response factor of the compound,
and $MW_{compound}$ is the molecular weight of the compound ($g.mol^{-1}$).

The molecular response factor used is described in details by de Saint Laumer et al., 2010 and is calculated using the following equations:

$\Delta H_{comb} = 11.06 + 103.57n_C + 21.85 n_H -48.18n_O +7.46n_N + 74.67n_S - 23.57n_F - 27.43n_{Cl} - 11.90n_{Br} - 2.04n_I +46.5n_{Si}$

$MRF = -0.0708 + 8.57 \times 10^{-4} \Delta H_{comb} + 1.27 \times 10^{-1} n_{Benz} + 6.18 \times 10^{-2}n_{Br}$

Where,
$\Delta H_{comb}$ is the combustion enthalpy of the compound,
$n_X$ is the number of element X in the chemical formula of the compound,
and $n_{Benz}$ is the number of benzene rings in the structure.

## Code structure

### Classes
In total, there are three classes (all of type _dataclass_) defined in the code for the better organization of the objects. 
The use of class-objects (instead of lists) negatively impacts the efficiency of the code but results in more readability. 
The classes used in the code are the following:

* Compound
* Calibrant
* Blob



```python
@dataclass(slots=True)
class Compound:
    name: str
    formula: str | None = None
    mol_wt: float | None = None
    cas: str | None = None
    combustion: float | None = None
    mrf: float | None = None
    grouping: float | None = None
    cnumber: int | None = None
    found_in_db: bool | None = None
    n_benzene: int | None = None
    elements: pd.DataFrame | None = None
    Tb: float | None = None
    dipole: float | None = None
```
The _Compound_ class is used for all the other objects since every calibrant and every blob also represents a compound. 
The only mandatory input for a Compound object is the name. If the name is not provided, a ValueError will be raised warning the user.
The compound object can fill up the necessary information using the _.search()_ method. To get the information, _compound_search()_ 
function (defined outside the scope of the class) is used. The function searches first through a database file (i.e., db) and 
if failed, in the NIST library file (nist).

#### _NEW_
The _Compound_ now also stores the boiling point (Tb) and the dipole moment of the compound. The boiling point is used to calculate the retention time of the compound.
This information is used to verify the peak assignments using the retention time plots. For this feature to work, please use 
the new database file (i.e., _db_Tb_dipole.csv_) containing the relevant information and install the chemicals package.
```python
@dataclass(slots=True)
class Calibrant:
    compound: Compound
    cal_type: str | None = None
    cal_quantity: np.ndarray | float | None = None
    cal_volume: np.ndarray | float | None = None
    path: str = 'calibration.csv'
    curve: list[float] = []
```
The _Calibrant_ class is used to store the information of the calibrants. Calibrants inherit their properties from the _Compound_ class.
By default, the internal calibrants (all of the same kind defined in the _extract_calibrants()_ method below) are limitless
and extracted from the blob table, but the external calibrants need to be added manually to the calibrants dictionary with 
keys that are not already used. In the main code, key 0 is already assigned to the external calibrant (uncomment the lines to use).
```python
blobs, calibrants, mass_closure = extract_calibrants(blobs, 'Internal Liquid', db, nist)
```

The class contains the following methods:
* _calibration_method()_: Can be used to ask for the calibration method in case it is not provided already.
* _calibration_curve()_: Calculates the slope and the intercept of the calibration curve. The curve is stored at self.curve.
By default, the calibration curve is assumed to be a line, but other functions can be set in the code, if necessary.
* _plot_calibration()_: Plots the calibration curve using the given amounts and responses and calibration types. Requires the curve to be already calculated.
```python
@dataclass(slots=True)
class Blob:
    compound: Compound
    volume: float
    retI: float = 0.0
    retII: float = 0.0
    inclusion: bool = field(default=True)
    internal_standard: int = field(default=0)
    intensity: float | None = None
    mol: float | None = None
    mass: float | None = None
    wt_yield: float | None = None
    elements: pd.DataFrame | None = None
```
The _Blob_ class is used to store the information of the blobs. The class also contains the _.process()_ method to calculate
the amount of the compound in the sample as well as static methods to normalize a blob list and to create a dataframe from the
blobs for the grouping of the results and presentation.

## References
de Saint Laumer, Jean-Yves, et al. "Quantification in gas chromatography: prediction of flame ionization detector response factors from combustion enthalpies and molecular structures." Analytical chemistry 82.15 (2010): 6457-6462.

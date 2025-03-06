# GCProcessing
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
In total, there are three classes (all of type *dataclass*) defined in the code for the better organization of the objects. The use of class-objects (instead of lists) negatively impacts the efficiency of the code but results in more readability. The classes used in the code are the following:

* Compound
* Calibrant
* Blob



```python
@dataclass(slots=True)
class Compound:
    name: str | None = None
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
```
The *Compound* class is used for all of the other objects since every calibrant and every blob also representes a compound. The only mandatory input for a Compound object is the name. If the name is not provided, a ValueError will be raised warning the user.

The compound object can fill up the necessary information using the *.search* method. To get the information, *compound_search* function (defined outside the namespace of the class) is used. The function searches first through a database file (i.e., db) and if failed, in the NIST library file (nist)

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

```python
@dataclass(slots=True)
class Blob:
    compound: Compound
    retI: float
    retII: float
    volume: float
    inclusion: bool = field(default=True)
    intensity: float | None = None
    mol: float | None = None
    mass: float | None = None
    wt_yield: float | None = None
    elements: pd.DataFrame | None = None
```

## References
de Saint Laumer, Jean-Yves, et al. "Quantification in gas chromatography: prediction of flame ionization detector response factors from combustion enthalpies and molecular structures." Analytical chemistry 82.15 (2010): 6457-6462.

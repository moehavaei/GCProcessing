from collections import defaultdict
from typing import Optional, Self
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from tkinter import filedialog as fd
from dataclasses import dataclass, field

from pandas import DataFrame
from xlsxwriter import Workbook
from rdkit import Chem
from rdkit.Chem import Descriptors, rdmolops, GetPeriodicTable
import re
from sklearn.linear_model import LinearRegression

# Configuring the logging settings
logging.basicConfig(filename='log.txt', level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Initializing the periodic table for MW calculation
periodic_table = GetPeriodicTable()


def plot_calibration(amounts: list | pd.Series, responses: list | pd.Series, cal_curve: list[float], cal_type: str) -> None:
    """
    Plots the calibration curve using the given amounts and responses and calibration types. Requires the curve to be already calculated.

    Example
    ========
    plot_calibration(amounts=[1, 2, 3], responses=[1, 2, 3], cal_type='External Gas', cal_curve=[1, 0])

    :param cal_type: Type of the calibration, i.e., "External Gas" or "External Liquid"
    :param amounts: The array containing the amounts for the calibration points (for gas, in barg).
    :param responses: The array containing the peak volumes for the calibration points.
    :param cal_curve: The slope and the intercept of the linear regression, i.e., [a, b] in y = ax + b.
    :return: None
    :raises: ValueError if the calibration type is not recognized.
    """
    plt.scatter(amounts, responses)
    match cal_type:
        case "External Liquid":
            x = np.linspace(0, max(amounts) * 1.05, 100)
            y = cal_curve[0] * 10 * x + cal_curve[1]
            plt.xlabel("Concentration\n[wt.%]")
        case "External Gas":
            x = np.linspace(0, max(amounts) * 1.05, 100)
            y = cal_curve[0] * x + cal_curve[1]
            plt.xlabel("Amount\n[µg]")
        case _:
            logging.error(f'Unknown calibration type: "{cal_type}" entered in the plot_calibration function.')
            raise ValueError(f'Unknown calibration type: "{cal_type}" entered in the plot_calibration function.')
    plt.plot(x, y, label=f'y = {cal_curve[0]:.2f} × x + {cal_curve[1]:.1f}', linestyle='--',
             linewidth=2,
             color='#EE964B')
    # plt.xlim(0, max(cal['amount']*1.05))
    plt.title("Calibration curve")

    plt.ylabel("Detector response")
    plt.legend()
    plt.savefig('Calibration.png', format='png', bbox_inches='tight')
    plt.show()


def load_data(path: str, file_name: str) -> tuple[DataFrame, str] | None:
    """
    Tries to open a file and if the file was not found, prompts the user (with a browse window) to choose the correct file.

    Example
    ========
    data, path = load_data('data.csv', 'Data')

    :param path: The default path where the file is located.
    :param file_name:The name of the file. This does not have to be limited to the actual name of the file!
    :return: A dataframe containing the content of the csv file, the path of the file, in case the file was not found.
    :raises: FileNotFoundError if the file was not found and the user did not choose a file.
    :raises: pd.errors.EmptyDataError if the file was empty.
    """
    try:
        return pd.read_csv(path), path
    except FileNotFoundError:
        logging.info(f'{file_name} was not found.')
        print(f'{file_name} was not found. Choose the correct path.')
        new_path: str = fd.askopenfilename(title=f'Select a file for the {file_name}',
                                           filetypes=(("CSV Files", "*.csv"), ("All", "*.*")))
        return pd.read_csv(new_path), new_path
    except pd.errors.EmptyDataError:
        logging.critical(f'{file_name} was empty.')
        pd.errors.EmptyDataError(f'{file_name} was empty.')


def blob_cleanup(raw_blobs: pd.DataFrame) -> pd.DataFrame | None:
    """
    Cleans up a blob table. First, the function removes excluded blobs, then it combines the same compounds by
    adding their volume and averaging their retention times. Finally, the function sorts the table by retention times
    and returns the cleaned-up blob table. When two internal standards are chosen for a single compound in two different peaks
    (by mistake), the larger value will be chosen.

    If the dataframe lacks the retention time columns, a warning will be given and the run proceeds.
    If the dataframe lacks the inclusion column, a warning will be given and the run proceeds assuming all compounds to be included.
    If the dataframe lacks the remaining columns, a KeyError will be raised.

    Example
    ========
    cleaned_blobs = blob_cleanup(raw_blobs)

    :param raw_blobs:
    :return: A cleaned-up blob table.
    :raises: KeyError if the dataframe lacks the Volume, or Internal Standard columns.
    """
    try:
        blobs = raw_blobs[raw_blobs['Inclusion'] == True]
        blobs = blobs.groupby('Compound Name', as_index=False).agg({'Retention I (min)': 'mean',
                                                                    'Retention II (sec)': 'mean',
                                                                    'Volume': 'sum',
                                                                    'Internal Standard': 'max'})
        blobs = blobs.sort_values(by=['Retention I (min)', 'Retention II (sec)'], ignore_index=True)
        return blobs
    except KeyError as e:
        if 'Retention I (min)' not in raw_blobs.columns or 'Retention II (sec)' not in raw_blobs.columns:
            logging.warning('Retention time(s) not provided in the blob table. The run proceeded assuming all retention times are 0.0.')
            print('Retention time(s) not provided in the blob table. The run will proceeded, but this is not recommended.')
            blobs = raw_blobs.groupby('Compound Name', as_index=False).agg({'Volume': 'sum'})
            return blobs
        if 'Inclusion' not in raw_blobs.columns:
            logging.warning('Inclusion column not provided in the blob table. The run proceeded assuming all blobs are included.')
            print('Inclusion column not provided in the blob table. The run will proceeded, but this is not recommended.')
            blobs = raw_blobs.groupby('Compound Name', as_index=False).agg({'Retention I (min)': 'mean',
                                                                    'Retention II (sec)': 'mean',
                                                                    'Volume': 'sum',
                                                                    'Internal Standard': 'max'})
            return blobs
        logging.error(e)
        raise KeyError(f'KeyError: {e} in blob_cleanup function. The run cannot proceed without the necessary columns.')


def formula_to_dataframe(formula: str) -> pd.DataFrame:
    """
    This function opens a string formula into a dictionary for easier access to the elements and the quantities

    Example
    ========
    formula_to_dataframe("C6H12")    # "C6H12" -> {'C': 6, 'H': 12}

    :param formula: A string containing the formula of the compound
    :return: A dictionary with the elements as the keys and the quantities as the values
    """
    element_quantity_pairs = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    data = {element: int(quantity) if quantity else 1 for element, quantity in element_quantity_pairs}
    return pd.DataFrame(data, index=[0])


def grouping(elements: pd.DataFrame, n_benzene: int) -> float:
    """
    This function determines the grouping of a compound based on its elemental composition and number of benzene rings.
    The function works according to a priority system. For instance, a chlorinated aromatic compound will be labeled
    as a halogenated compound (i.e., grouping 7). Any saturated compound will be labeled as a i-paraffin (i.e., grouping 0.9).
    For one degree of unsaturation, the compound is assumed to be an i-olefin (i.e., grouping 1.55). For two degrees of
    unsaturation, the compound is assumed to be a diene (i.e., grouping 1.7). These assumptions are made assuming that the typical
    compounds already exist in the database. If this does not apply to specific compounds, the grouping can be manually changed in
    the database after updating it and the code should be re-run. Additionally, bi+-phenyls and similar compounds cannot be
    distinguished from polyaromatic compounds.

    :param elements: A dictionary containing the elemental composition of the compound.
    :param n_benzene: The number of benzene rings in the compound.
    :return: Grouping of the compound.
    """

    conditions = {
        'Cl': 7,
        'F': 7,
        'Br': 7,
        'I': 7,
        'O': 3,
        'N': 4,
        'S': 5
    }

    for element, value in conditions.items():
        if elements[element].iloc[0] > 0:
            return float(value)

    if n_benzene == 1:
        return 2
    elif n_benzene == 2:
        return 2.2
    elif n_benzene > 2:
        return 2.3

    unsaturation = (elements['C'].iloc[0] * 2 + 2 - elements['H'].iloc[0]) / 2
    unsaturation_map = {0: 0.9, 1: 1.55, 2: 1.7}

    return unsaturation_map.get(unsaturation, None)


def count_benzene_rings(inchi: str) -> int:
    """
    This function takes an InChI string as input and returns the number of benzene rings in the molecule.

    :param inchi: A string containing the InChI string
    :return: The number of benzene rings in the molecule
    """
    # Convert InChI to molecule object
    mol = Chem.MolFromInchi(inchi)

    if mol is None:
        return 0  # Return 0 if the InChI string cannot be parsed

    # Get all aromatic rings in the molecule
    aromatic_rings = rdmolops.GetSymmSSSR(mol)

    benzene_rings_count = 0

    # Loop through the rings to find benzene rings (6-membered aromatic rings)
    for ring in aromatic_rings:
        # Check if the ring has 6 atoms and all atoms are aromatic carbons
        if len(ring) == 6 and all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            # Check if all atoms in the ring are carbon
            if all(mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C' for atom_idx in ring):
                benzene_rings_count += 1

    return benzene_rings_count


def calculate_mrf(element_composition: pd.DataFrame, n_benzene: int) -> tuple[float, float]:
    """
    Calculates the enthalpy of combustion and molecular response factor for a given compound using the elemental
    composition of the compound and the number of benzene rings in the molecule.
    ref: https://doi.org/10.1002/jssc.201500106

    :param element_composition: DataFrame containing the elemental composition of the compound, i.e., elements as
    columns and number of the elements as the value.
    :param n_benzene: Number of benzene in the compound as an integer.
    :return: enthalpy of combustion (kJ/mol), molecular response factor (1/mol)
    """
    combustion: float = (11.06 + 103.57 * element_composition['C'].iloc[0] + 21.85 * element_composition['H'].iloc[0]
                         - 48.18 * element_composition['O'].iloc[0] + 7.46 * element_composition['N'].iloc[0]
                         + 74.67 * element_composition['S'].iloc[0] - 23.57 * element_composition['F'].iloc[0]
                         - 27.43 * element_composition['Cl'].iloc[0] - 11.90 * element_composition['Br'].iloc[0]
                         - 2.04 * element_composition['I'].iloc[0] + 46.5 * element_composition['Si'].iloc[0])
    mrf: float = (-0.0708 + 8.57e-4 * combustion + 1.27e-1 * n_benzene
                  + 6.18e-2 * element_composition['Br'].iloc[0])
    return combustion, mrf


def translate_grouping(grouping_number: float) -> str:
    """
    Translates the number grouping into text (group names). For instance, 1.0 -> Paraffin
    The values and the numbers can be changed based on the specific characteristics of the samples.
    This function can be discarded if the database contains group names from scratch. However, modifications
    need to be made wherever this function is called, and where the grouping is read from the database as a
    float number.

    Example
    =======
    translate_grouping(1.0)    # 1.0 -> 'Paraffin'

    :param grouping_number: Grouping number
    :return: Grouping text
    """
    group_names: dict = {
        1.0: 'Paraffin',
        1.1: 'Naphthene',
        1.11: 'Dinaphthene',
        1.12: 'Naphthenoaromatic',
        1.5: 'Olefin',
        2: 'MAH',
        2.2: 'DAH',
        2.3: 'PAH',
        2.02: 'Di-phenyl',
        2.03: 'Tris+-phenyl',
        3: 'Oxygenated',
        4: 'Nitrogenated',
        5: 'Sulfurinated',
        7: 'Halogenated',
        0.9: 'i-Paraffin',
        1.55: 'i-Olefin',
        1.7: 'Diene',
        0.0: 'Others'
    }
    try:
        return group_names[grouping_number]
    except KeyError:
        return 'Others'


def calculate_mol_wt(formula_broken: pd.DataFrame) -> float:
    """
    Calculates the molecular weight of the compound using a DataFrame of the formula containing the elements and their number.

    :param formula_broken: Pandas DataFrame containing elements as the column and their number as the value in index 0.
    :return: Molecular weight of the compound.
    """
    mol_wt: float = 0.0
    for element in formula_broken.columns:
        mol_wt += periodic_table.GetAtomicWeight(element) * formula_broken.loc[0, element]
    return mol_wt


def prompt_compound(name: str) -> pd.DataFrame:
    """
    Asks the user for compounds not found in NIST library or the database.
    The user is prompted for the chemical formula, number of benzene rings.

    :param name: Name of the compound.
    :return: Relevant information of the compound.
    """
    print(f"Compound {name} was not found in the database or the NIST library.")
    logging.info(f"Compound {name} was not found in the database or the NIST library.", )
    print("What is the chemical formula of the compound? Enter the formula in the standard format, i.e., C5H9Cl\n")
    formula: str = str(input())
    formula_broken: pd.DataFrame = formula_to_dataframe(formula)
    mol_wt: float = calculate_mol_wt(formula_broken)
    print(f"How many benzene rings does {name} have?\n")
    try:
        n_benzene: int = int(input())
    except ValueError:
        logging.info(f"ValueError: entered value for {name} benzene rings was not a number.")
        print("Please enter a valid number.")
        n_benzene: int = int(input())
    elements: list[str] = ['C', 'H', 'O', 'N', 'Cl', 'S', 'F', 'Si', 'Br', 'I']
    element_composition: pd.DataFrame = pd.DataFrame(columns=elements)
    element_composition.loc[0] = 0
    for element in (set(elements) & set(formula_broken.columns)):
        element_composition.loc[0, element] = formula_broken.loc[0, element]
    combustion, mrf = calculate_mrf(element_composition, n_benzene)
    group: float = grouping(element_composition, n_benzene)
    compound_info: pd.DataFrame = pd.DataFrame({'formula': [formula],
                                                'grouping': [group],
                                                'mol_wt': [mol_wt],
                                                'n_benzene': [n_benzene],
                                                'combustion': [combustion],
                                                'mrf': [mrf]})
    compound_info = pd.concat([compound_info, element_composition], axis=1)
    return compound_info


def compound_search(*, name, db: pd.DataFrame, nist: pd.DataFrame) -> tuple[
    pd.DataFrame, bool] | None:
    """
    Searches a compound first in a database, and then in the NIST library (both locally available).
    If the compound was not found in either, the user will be prompted for the compound info.
    To use this function, all the keywords must be passed to avoid confusing the nist file and the db file.

    Example
    =========
    compound_search(name='Butane', db=db, nist=nist)

    :param name: Name of the compound. i.e., 'Butane'. The search is case-insensitive.
    :param db: DataFrame containing the database.
    :param nist: DataFrame containing the NIST database.
    :return: If found, the DataFrame containing the relevant information about the compound and
    a boolean indicating if the compound was found.
    """
    found_in_db: bool = False
    db_search: pd.DataFrame = db[db['compound'].str.lower() == name.lower()]
    elements: list[str] = ['C', 'H', 'O', 'N', 'Cl', 'S', 'F', 'Si', 'Br', 'I']
    element_quantities: pd.DataFrame = pd.DataFrame([[0] * len(elements)], columns=elements)
    formula: str = ''
    mol_wt: float = 0.0
    n_benzene: int = 0
    combustion: float = 0.0
    mrf: float = 0.0
    group: float = 0.0
    if not db_search.empty:
        try:
            found_in_db = True
            formula = db_search['formula'].iloc[0]
            formula_broken: pd.DataFrame = formula_to_dataframe(formula)
            mol_wt = float(db_search['MW'].iloc[0])
            n_benzene = int(db_search['n_Benz'].iloc[0])
            combustion = float(db_search['Combust'].iloc[0])
            mrf = float(db_search['MRF'].iloc[0])
            group = float(db_search['grouping'].iloc[0])
            for element in (set(elements) & set(formula_broken.columns)):
                element_quantities[element] = db_search[element].iloc[0]
        except ValueError as e:
            logging.warning(f"Compound {name} had unacceptable values in the database")
            raise ValueError(f"Compound {name} had unacceptable values in the database: {e}")
    else:
        nist_search = nist[nist['name'].str.lower() == name.lower()]
        if nist_search.empty:
            nist_search = nist[nist['synonyms'].str.lower() == name.lower()]
        if nist_search.empty:
            nist_search = nist[nist['synonyms'].str.contains(fr"\b{name}\b", case=False, na=False)]
        if nist_search.empty:
            nist_search = nist[nist['synonyms'].str.contains(fr"{name}\b", case=False, na=False)]
        if nist_search.empty:
            nist_search = nist[nist['synonyms'].str.contains(fr"\b{name}", case=False, na=False)]
        if nist_search.empty:
            compound_info = prompt_compound(name)
            return compound_info, found_in_db
        else:
            formula = nist_search['formula'].iloc[0]
            mol_wt = float(nist_search['mol_weight'].iloc[0])
            n_benzene = count_benzene_rings(nist_search['inchi'].iloc[0])
            formula_broken = formula_to_dataframe(formula)
            for element in (set(elements) & set(formula_broken.columns)):
                element_quantities[element] = formula_broken[element].iloc[0]
            group = grouping(element_quantities, n_benzene)
            combustion, mrf = calculate_mrf(element_quantities, n_benzene)
    compound_info: pd.DataFrame = pd.DataFrame({'formula': [formula],
                                                'grouping': [group],
                                                'mol_wt': [mol_wt],
                                                'n_benzene': [n_benzene],
                                                'combustion': [combustion],
                                                'mrf': [mrf]})
    compound_info = pd.concat([compound_info, element_quantities], axis=1)
    return compound_info, found_in_db


def append_to_csv(file_path: str, df_to_append: pd.DataFrame) -> None:
    """
    Appends the pandas DataFrame with the new compounds to the database csv file.

    :param file_path: CSV file path
    :param df_to_append: The pandas DataFrame to append
    :return:
    """

    # Reads the existing data
    existing_data = pd.read_csv(file_path)

    # Ensure columns in df_to_append match the existing Excel columns
    matching_df = df_to_append.reindex(columns=existing_data.columns)
    matching_df = matching_df.dropna(axis=1, how='all')

    # Combines the new and existing data
    combined_data = pd.concat([existing_data, matching_df], ignore_index=True)

    # Write the updated data back to the CSV file
    combined_data.to_csv(file_path, index=False)


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

    def __post_init__(self):
        if not isinstance(self.name, str):
            logging.critical('Name of the compound needs to be a string. Please provide a valid name.')
            raise ValueError('Name of the compound needs to be a string. Please provide a valid name.')

    def search(self, db: pd.DataFrame, nist: pd.DataFrame) -> None:
        """
        Searches for the compound using the 'compound_search()' function in order to find the compound info.

        :param db: DataFrame containing the database.
        :param nist: DataFrame containing the NIST database.
        :return:
        """
        compound_info: pd.DataFrame | None
        found_in_db: bool
        compound_info, found_in_db = compound_search(name=self.name, db=db, nist=nist)
        if compound_info is not None:
            self.formula = compound_info['formula'].iloc[0]
            self.mol_wt = compound_info['mol_wt'].iloc[0]
            # self.cas = compound_info['cas'].iloc[0]
            self.combustion = compound_info['combustion'].iloc[0]
            self.mrf = compound_info['mrf'].iloc[0]
            self.grouping = compound_info['grouping'].iloc[0]
            self.elements = pd.DataFrame(columns=['C', 'H', 'O', 'N', 'Cl', 'S', 'F', 'Si', 'Br', 'I'])
            for element in self.elements.columns:
                self.elements.loc[0, element] = compound_info[element].iloc[0]
            self.n_benzene = compound_info['n_benzene']
            self.cnumber = compound_info['C'].iloc[0]
            self.found_in_db = found_in_db


@dataclass(slots=True)
class Calibrant:
    compound: Compound
    cal_type: str | None = None
    cal_quantity: pd.Series | float | None = None
    cal_volume: pd.Series | float | None = None
    path: str = 'calibration.csv'
    curve: list = field(default_factory=list)

    def calibration_method(self) -> None:
        """
        Can be used to ask for the calibration method in case it is not provided already.
        :return:
        """
        if self.cal_type is None:
            print(
                'You have not chosen a calibration method. By default, internal calibration using a gas is assumed.\n')
            internal_verification: str = str(input('Would you like to choose another method? (y, n) \n')).lower()
            while internal_verification not in ['y', 'n', '']:
                internal_verification = str(input('Invalid value. Please, enter "y" or "n". \n')).lower()
            if internal_verification == 'n':
                self.cal_type = 'Internal Gas'
            else:
                print(
                    'Which one of the following methods would you like to use? (a, b, c)\n ',
                    'a) External calibration with a gas (points provided in a csv file titled "calibration.csv")\n',
                    'b) External calibration with a liquid (points provided in a csv file titled "calibration.csv")\n',
                    'c) Internal calibration in a liquid (values will be prompted)')
                calibration_prompt: str = str(input()).lower()
                while calibration_prompt not in ['a', 'b', 'c']:
                    print('Invalid value. Please, enter "a", "b", or "c". \n')
                    calibration_prompt: str = str(input()).lower()

                match calibration_prompt:
                    case 'a':
                        self.cal_type = 'External Gas'
                    case 'b':
                        self.cal_type = 'External Liquid'
                    case 'c':
                        self.cal_type = 'Internal Liquid'
                    case _:
                        logging.critical('Calibration type was not recognized in the "calibration_method()" function. Please, choose a valid calibration type.')
                        raise ValueError('Calibration type was not recognized. Please, choose a valid calibration type.')

    def calibration_curve(self, intercept: bool = False) -> None:
        """
        Calculates the slope and the intercept of the calibration curve. The curve is stored at self.curve.

        :param intercept: Whether an intercept should be calculated or not. By default, the line
        passes through the origin.
        :return:
        """
        match self.cal_type:
            case 'Internal Gas':
                if self.cal_quantity is None: self.cal_quantity = float(
                    input(f'What is the quantity of the calibrant {self.compound.name}?\n'))
                if self.cal_volume is None: self.cal_volume = float(
                    input(f'What is the volume of the calibrant peak {self.compound.name}?\n'))
                try:
                    quantity = 0.05 * (self.cal_quantity + 1.01325) * 1e5 * 0.000000118278 / (
                        8.314 * (273.15 + 50)) * self.compound.mol_wt * 1e6
                    self.curve = [self.cal_volume / quantity, 0.0]
                except ZeroDivisionError:
                    print(f'Quantity of the calibrant {self.compound.name} cannot be zero!')
                    logging.error(f'Quantity of the calibrant {self.compound.name} cannot be zero!')
            case 'Internal Liquid':
                if self.cal_quantity is None: self.cal_quantity = float(
                    input(f'What is the concentration (wt.%) of the calibrant {self.compound.name}?\n'))
                if self.cal_volume is None: self.cal_volume = float(
                    input(f'What is the volume of the calibrant peak {self.compound.name}?\n'))
                try:
                    self.curve = [self.cal_volume / self.cal_quantity, 0.0]
                except ZeroDivisionError:
                    print('Concentration of the calibrant cannot be zero!')
                    logging.error(f'Quantity of the calibrant {self.compound.name} cannot be zero!')
            case 'External Gas':
                calibration_table, _ = load_data(self.path, 'Calibration')
                cal_curve = LinearRegression(fit_intercept=intercept)
                self.cal_quantity = 0.05 * (calibration_table.iloc[:, 0] + 1.01325) * 1e5 * 0.000000118278 / (
                        8.314 * (273.15 + 50)) * self.compound.mol_wt * 1e6
                self.cal_volume = calibration_table.iloc[:, 1]
                cal_curve.fit(self.cal_quantity.values.reshape(-1, 1), self.cal_volume)
                self.curve = [float(cal_curve.coef_[0]), cal_curve.intercept_]
            case 'External Liquid':
                calibration_table = pd.read_csv(self.path)
                cal_curve = LinearRegression(fit_intercept=intercept)
                self.cal_quantity = calibration_table.iloc[:, 0]
                self.cal_volume = calibration_table.iloc[:, 1]
                cal_curve.fit(self.cal_quantity.values.reshape(-1, 1) * 10.0,
                              self.cal_volume.values)  # Multiplied by 10 to convert wt.% to µg.
                self.curve = [float(cal_curve.coef_[0]), cal_curve.intercept_]


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

    def process(self, calibrant: Calibrant, sample_amount: float = 0.1) -> None:
        """
        This function uses the calibration curve of a specific calibrant and calculates the quantity of each blob
        as well as the elemental composition. In the end, the database csv file is appended and
        the -possibly- updated database is returned.

        :param calibrant: Calibrant as a Calibrant object. The calibration curve needs to be calculated using
        the 'calibration_curve' method before running this function.
        :param sample_amount: Amount of samples in mg. For liquid injections, must be set on 0.1, which is the default value.
        :return:
        """

        if calibrant.cal_type == ('Internal Liquid' or 'External Liquid'):
            sample_amount = 0.1
        self.mol = (((self.volume - calibrant.curve[1]) / calibrant.curve[0]) * calibrant.compound.mrf) / (
                calibrant.compound.mol_wt * 1e6 * self.compound.mrf)
        self.mass = self.mol * self.compound.mol_wt * 1000
        self.wt_yield = self.mass / sample_amount * 100
        self.elements = self.compound.elements.copy()
        for element in self.elements.columns:
            self.elements.loc[0, element] = self.compound.elements.loc[0, element] * periodic_table.GetAtomicWeight(
                element) / self.compound.mol_wt * self.wt_yield

    @staticmethod
    def normalize_blob_list(blob_list: list['Blob'], mass_closure: float, calibrants: dict[int, Calibrant]) -> tuple[list['Blob'], bool]:
        """
        Normalizes the yields of the blobs in a blob list.
        :param calibrants: A dictionary containing the calibrants.
        :param blob_list: List of Blob objects.
        :param mass_closure: The mass closure of the blobs before normalization.
        :return: The normalized blob list and a True boolean declaring that the yields have been normalized.
        """

        for calibrant in calibrants.values():
            if calibrant.cal_type == 'Internal Liquid':
                mass_closure = mass_closure - calibrant.cal_quantity
        for blob in blob_list:
            blob.wt_yield = blob.wt_yield / mass_closure * 100
            for element in blob.elements.columns:
                blob.elements.loc[0, element] = blob.elements.loc[0, element] / mass_closure * 100
        return blob_list, True

    @staticmethod
    def blob_list_to_dataframe(blob_list: list['Blob']) -> pd.DataFrame:
        """
        Populates a pandas dataframe from a list of Blob objects. It extracts the name, C#, group name, yield
        and the elemental composition of the blobs.
        :param blob_list: List of Blob objects.
        :return: Pandas DataFrame containing the important information of the blobs.
        """
        elements: list = ['C', 'H', 'O', 'N', 'Cl', 'S', 'F', 'Si', 'Br', 'I']
        blob_df = pd.DataFrame(columns=['Compound Name', 'C#', 'Group Name', 'Yield [wt.%]'] + elements)
        for i, blob in enumerate(blob_list):
            blob_df.loc[i, 'Compound Name'] = blob.compound.name
            blob_df.loc[i, 'C#'] = blob.compound.cnumber
            blob_df.loc[i, 'Group Name'] = translate_grouping(blob.compound.grouping)
            blob_df.loc[i, 'Yield [wt.%]'] = blob.wt_yield
            blob_df.loc[i, 'Volume'] = blob.volume
            for element in elements:
                blob_df.loc[i, element] = blob.elements.loc[0, element]

        return blob_df


def update_db(db: pd.DataFrame, path: str, blob_list: list[Blob]) -> None:
    """
    Updates the database from the compounds not found in the database (using 'found_in_db' flag).
    :param db: DataFrame containing the database to avoid re-reading.
    :param path: Path of the database.
    :param blob_list: List of Blob objects in order to find 'not-found' compounds.
    :return:
    """
    new_compounds: pd.DataFrame = pd.DataFrame(columns=db.columns)
    for i, blob in enumerate(blob_list):
        if not blob.compound.found_in_db:
            new_compounds.loc[i, 'compound'] = blob.compound.name
            new_compounds.loc[i, 'grouping'] = blob.compound.grouping
            new_compounds.loc[i, 'MW'] = blob.compound.mol_wt
            new_compounds.loc[i, 'formula'] = blob.compound.formula
            new_compounds.loc[i, 'n_Benz'] = blob.compound.n_benzene[0]
            new_compounds.loc[i, 'Combust'] = blob.compound.combustion
            new_compounds.loc[i, 'MRF'] = blob.compound.mrf
            elements: list[str] = ['C', 'H', 'O', 'N', 'Cl', 'S', 'F', 'Si', 'Br', 'I']
            for element in elements:
                new_compounds.loc[i, element] = blob.compound.elements.loc[0, element]

    append_to_csv(path, new_compounds)

def piona_table(blob_df: pd.DataFrame) -> pd.DataFrame:
    """
    This function creates a table with the yields of the different groups of compounds for a PIONA table (paraffins, naphthenes, etc.)

    :param blob_df: DataFrame containing the blob information.
    :return: PIONA table as a DataFrame.
    """
    piona_cols = ['Paraffins', 'i-Paraffins', 'Olefins', 'Naphthenes', 'Aromatics']
    piona_labels = pd.DataFrame(columns=piona_cols)
    piona_labels[piona_cols[0]] = (blob_df['Group Name'] == "Paraffin")
    piona_labels[piona_cols[1]] = (blob_df['Group Name'] == "i-Paraffin")
    piona_labels[piona_cols[2]] = (blob_df['Group Name'] == "Olefin")
    piona_labels[piona_cols[3]] = ((blob_df['Group Name'] == "Naphthene") | (blob_df['Group Name'] == "Dinaphthene"))
    piona_labels[piona_cols[4]] = ((blob_df['Group Name'] == "MAH") | (blob_df['Group Name'] == "DAH")
                                   | (blob_df['Group Name'] == "PAH") | (blob_df['Group Name'] == "di-phenyl")
                                   | (blob_df['Group Name'] == "tris+-phenyl") | (
                                           blob_df['Group Name'] == "Naphthenoaromatic"))

    z = np.zeros((6, 5), float)
    piona = pd.DataFrame(z, columns=piona_cols)

    conditions = [(blob_df['C#'] > 1) & (blob_df['C#'] < 5),
                  (blob_df['C#'] > 4) & (blob_df['C#'] < 12),
                  (blob_df['C#'] > 11) & (blob_df['C#'] < 21),
                  (blob_df['C#'] > 20) & (blob_df['C#'] < 36),
                  (blob_df['C#'] > 35) & (blob_df['C#'] < 66)]
    for label in piona_cols:
        for i in range(5):
            piona.loc[i, label] = blob_df[(conditions[i]) & (piona_labels[label])].agg(
                {'Yield [wt.%]': 'sum'}).iloc[0]
        piona.loc[5, label] = piona[label].sum()

    piona_labels = ['C2-C4', 'C5-C11', 'C12-C20', 'C21-C35', 'C36-C65', 'Total']
    piona.index = piona_labels

    return piona

def extract_calibrants(blobs: pd.DataFrame, db: pd.DataFrame, nist: pd.DataFrame) -> tuple[pd.DataFrame, dict[int, Calibrant], float]:
    try:
        ISTDs_df: pd.DataFrame = blobs[(blobs['Amount'] > 0) & (blobs['Internal Standard'] > 0)]
        indexes = ISTDs_df.index
        blobs.drop(indexes, inplace=True)
        total_calibrant: float = ISTDs_df['Amount'].sum()
    except KeyError as e:
        print(
            'Amount column or Internal Standard column not included in the blob table or entered with a different name '
            '(e.g., Amount (wt.%)). For using this feature, please add the columns to the blob table.')
        logging.error(
            'Amount column or Internal Standard column not included in the blob table or entered with a different name '
            '(e.g., Amount (wt.%)). For using this feature, please add the columns to the blob table.')
        raise e
    calibrants: dict[int, Calibrant] = defaultdict(None)
    for i, ISTD in ISTDs_df.iterrows():
        cal_compound: Compound = Compound(ISTD['Compound Name'])
        cal_compound.search(db, nist)
        calibrants[ISTD['Internal Standard']] = (Calibrant(cal_compound, cal_type='Internal Liquid',
                                                           cal_volume=ISTD['Volume'], cal_quantity=ISTD['Amount']))
        calibrants[ISTD['Internal Standard']].calibration_method()
        calibrants[ISTD['Internal Standard']].calibration_curve()

    return blobs, calibrants, total_calibrant

if __name__ == "__main__":
    ...
    # comp: Compound = Compound()
    # db, path = load_data('db.csv', 'database')
    # nist, path = load_data('nist_compounds.csv', 'NIST library')
    # comp.search(db, nist)
    # print(comp)

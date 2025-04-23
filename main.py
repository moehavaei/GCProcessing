import os
import logging
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from blob_processing import *
import win32com.client as win32
from itertools import permutations
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

pd.set_option('future.no_silent_downcasting', True)
# Configuring the logging settings
logging.basicConfig(filename='log.txt', level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')



def main() -> None:
    # Loading the files.
    current_dir: str = os.getcwd()
    path_db, path_blobs, path_nist = current_dir + r'\db_Tb_dipole.csv', current_dir + r'\blob_table5.csv', current_dir + r'\nist_compounds.csv'
    db, path_db = load_data(path_db, 'Database')
    blobs, path_blobs = load_data(path_blobs, 'Blobs')
    blob_columns = ['Compound Name', 'Retention I (min)', 'Retention II (sec)', 'Volume', 'Inclusion']
    for col in blob_columns:
        if col not in blobs.columns:
            if col in ['Retention I (min)', 'Retention II (sec)']:
                blobs.drop(columns=['Retention I (min)', 'Retention I'])
            else:
                logging.error(f'Column {col} not found in the blob table.')
            raise KeyError(f'Column {col} not found in the blob table.')
    nist, path_nist = load_data(path_nist, 'NIST')

    """
    Provide the calibration information below. You can use one of the two options:

    1. First blob in the table is your internal calibrant with its "Inclusion" set to False:
        Use the first three lines and comment the next three lines.

    2. Manual input of the calibrant:
        Use the second line and comment the first line.
    """

    # Finding internal standards
    calibrants: dict[int, Calibrant] = defaultdict(Calibrant)
    blobs, calibrants, mass_closure = extract_calibrants(blobs, 'Internal Liquid', db, nist)


    """
    For external calibration (as a validation for the internal calibration), use the lines below. You may also use the external
    calibration curve as the primary calibration curve by adding '0' in the Internal Standard column of the blob table.
    """

    # external_calibrant_comp: Compound = Compound('n-Hexane')
    # external_calibrant_comp.search(db, nist)
    # external_calibrant: Calibrant = Calibrant(external_calibrant_comp, cal_type='External Liquid', path='calibration.csv')
    # external_calibrant.calibration_method()
    # external_calibrant.calibration_curve()
    # calibrants[0] = external_calibrant
    #
    # # Plotting the calibration curve (for external calibrations)
    # calibrants[0].plot_calibration()


    # Agglomerating redundant blobs and removing blobs not intended for inclusion
    blobs = blob_cleanup(blobs)

    # Creating a list of processed blobs:
    blob_list: list[Blob] = []

    # Processing the blobs:
    """
    Enter the sample amount in the line below:
    """
    sample_amount: float = 1.0

    for i in tqdm.tqdm(blobs.index):
        blob = blobs.loc[i]
        compound = Compound(blob['Compound Name'])
        try:
            compound.search(db, nist)
        except Exception as e:
            print(f"Error searching for compound {compound}: {e}")
            logging.error(f"Error searching for compound {compound}: {e}")
        processed_blob = Blob(compound, retI=blob['Retention I (min)'], retII=blob['Retention II (sec)'],
                              volume=blob['Volume'], inclusion=True, internal_standard=blob['Internal Standard'],)
        processed_blob.process(calibrants[processed_blob.internal_standard], sample_amount=sample_amount)
        print(processed_blob.compound.Tb)
        mass_closure += processed_blob.wt_yield
        blob_list.append(processed_blob)

    # Updating the database with the blobs that were not found in the database:
    update_db(db, path_db, blob_list)


    """
    Normalizing the blobs in the list. Calibrant needs to be entered so that if Internal Liquid has been used, 
    the quantity of the calibrant is taken into account:
    """

    normalized = True
    blob_list, normalized = Blob.normalize_blob_list(blob_list=blob_list, mass_closure=mass_closure, calibrants=calibrants)
    # Populating a DataFrame with the useful information from the blobs:
    blob_df = Blob.blob_list_to_dataframe(blob_list)
    blob_df.plot(kind='scatter', x='Retention I (min)', y='Tb', color='red', legend=False)
    blob_df.plot(kind='scatter', x='Retention II (sec)', y='Dipole', color='blue', legend=False)

    # Colors used in the graphs:
    colors: list[str] = [
        '#008080',  # Teal
        '#FF6F61',  # Coral
        '#FFD700',  # Gold
        '#6A5ACD',  # Slate Blue
        '#DC143C',  # Crimson
        '#3CB371',  # Medium Sea Green
        '#4169E1',  # Royal Blue
        '#FF1493',  # Deep Pink
        '#B8860B',  # Dark Goldenrod
        '#4B0082',  # Indigo
        '#999999',  # Gray
        '#70163C',  # Tyrian purple
        '#95B2B8',  # Cadet gray
        '#F4D1AE'  # Light orange
    ]

    """
    The section below checks the mass closure of the sample. If the mass closure is within 5% of 100%, a ✅ is printed,
    otherwise a ❌ is printed.
    """
    if 110 > mass_closure > 90:
        print(f'Mass closure before normalization: {mass_closure:.2f} wt.%  \u2705')
    else:
        print(f'Mass closure before normalization: {mass_closure:.2f} wt.%  \u274C')

    """
    Verification of the ISTDs.
    """
    if len(calibrants) > 1:
        plt.rcParams["font.family"] = "Times New Roman"
        ISTD_validation: dict[tuple[int, int], float] = {}
        for i, j in list(permutations(calibrants.keys(), 2)):
            if calibrants[i].cal_type not in ['External Liquid', 'External Gas']:
                blob_ISTD: Blob = Blob(calibrants[i].compound, volume=calibrants[i].cal_volume)
                blob_ISTD.process(calibrant=calibrants[j], sample_amount=sample_amount)
                error: float = (blob_ISTD.wt_yield - calibrants[i].cal_quantity) / calibrants[i].cal_quantity * 100
                ISTD_validation[(i, j)] = error
        labels: list = []
        for i, j in ISTD_validation.keys():
            labels.append(f'{calibrants[i].compound.name} using {calibrants[j].compound.name}')
        error_max = np.max(np.abs(list(ISTD_validation.values())))

        # Custom diverging colormap: red-green-red
        color_gradient: list[tuple[int, int, int]] = [(1, 0, 0), (0, 1, 0), (1, 0, 0)]  # red → green → red
        cmap: LinearSegmentedColormap = LinearSegmentedColormap.from_list("red-green-red", color_gradient, N=256)
        norm: TwoSlopeNorm = TwoSlopeNorm(vmin=-50, vcenter=0, vmax=50)

        fig, ax = plt.subplots(figsize=(9, len(ISTD_validation) * 1.0))
        bars = ax.barh(list(map(str, labels)), list(ISTD_validation.values()), color=cmap(norm(list(ISTD_validation.values()))))
        ax.set_xlabel('Error [%]', fontsize='x-large')
        ax.axvline(0, color='gray', linestyle='--', linewidth=1)
        plt.tight_layout()
        plt.savefig(f'{path_blobs.removesuffix(".csv")}_ISTD_validation.svg', bbox_inches='tight', format='svg')
        plt.show()

    """
    The section below generates the output of the code including the graphs and the Excel file.
    """

    # Plotting the combined data using subplots (2×2):
    fig, axs = plt.subplots(2, 2, figsize=(15, 12))

    # Pie chart for elemental composition of the sample:
    elemental_composition = blob_df[['C', 'H', 'O', 'N', 'Cl', 'S', 'F', 'Si', 'Br', 'I']].sum()
    non_zero_elements = elemental_composition[elemental_composition > 0]
    explode = [0.1 * i if s < 2 else 0 for i, s in enumerate(non_zero_elements.values)]
    axs[0, 0].pie(non_zero_elements, autopct='%1.1f%%', labels=non_zero_elements.keys(), colors=colors,
                  radius=0.8, labeldistance=1.5, pctdistance=1.3, explode=explode)
    axs[0, 0].set_title('Elemental Composition', fontsize='xx-large', fontweight='bold')

    # Stacked bar plot of the group-type C# distribution:
    group_name_yields = blob_df.pivot_table(index='C#', columns='Group Name', values='Yield [wt.%]', aggfunc='sum')
    group_name_yields.plot(kind='bar', stacked=True, ax=axs[0, 1], rot=0, color=colors)
    axs[0, 1].legend(title='Group', title_fontsize='x-large')
    axs[0, 1].set_title('Grouped carbon number distribution', fontsize='xx-large', fontweight='bold')
    axs[0, 1].set_xlabel('C#', fontsize='x-large')
    axs[0, 1].set_ylabel('Yield [wt. %]', fontsize='x-large')
    axs[0, 1].set_xticks(np.arange(1, max(group_name_yields.index) + 1, 2))

    # Carbon number distribution
    grouped = blob_df.groupby('C#').sum()
    axs[1, 0].bar(grouped.index, grouped['Yield [wt.%]'], label='Yield [wt. %]', color=colors[6])
    axs[1, 0].legend(title_fontsize='x-large')
    axs[1, 0].set_title('Carbon number distribution', fontsize='xx-large', fontweight='bold')
    axs[1, 0].set_xlabel('C#', fontsize='x-large')
    axs[1, 0].set_ylabel('Yield [wt. %]', fontsize='x-large')
    axs[1, 0].set_xticks(np.arange(1, max(grouped.index) + 1, 2))

    # PIONA bubble chart
    piona: pd.DataFrame = piona_table(blob_df)
    y_labels: list[str] = piona.index[0:5]
    x_labels: list[str] = piona.columns
    x, y = np.meshgrid(range(len(x_labels)), range(len(y_labels)))
    bubble_color = piona.iloc[0: 5].values.flatten()
    bubble_size = (bubble_color) / max(bubble_color) * 10
    bubble_chart = axs[1, 1].scatter(x.flatten(), y.flatten(), s=bubble_size * 50, c=bubble_color, alpha=0.8,
                                     edgecolors='w', cmap='viridis')
    axs[1, 1].grid(True, which='both', linestyle='--', linewidth=0.5)
    axs[1, 1].set_xticks(range(len(x_labels)))
    axs[1, 1].set_yticks(range(len(y_labels)))
    axs[1, 1].set_xticklabels(x_labels, rotation=45, ha="right")
    axs[1, 1].set_yticklabels(y_labels)
    axs[1, 1].set_xlabel('Group Name')
    cbar = plt.colorbar(bubble_chart, ax=axs[1, 1])
    cbar.set_label('Yield [wt.%]')
    axs[1, 1].grid(True, which='both', linestyle='--', linewidth=0.5)
    axs[1, 1].set_title('PIONA', fontsize='xx-large', fontweight='bold')

    # Saving the combined charts
    plt.savefig(f'{path_blobs.removesuffix(".csv")}_graphs.svg', bbox_inches='tight', format='svg')
    plt.show()

    # Saving individual plots

    plt.pie(non_zero_elements, autopct='%1.1f%%', labels=non_zero_elements.keys(), colors=colors,
            radius=0.8, labeldistance=1.5, pctdistance=1.3, explode=explode)
    plt.title('Elemental Composition', fontsize='xx-large', fontweight='bold')
    plt.savefig(f'{path_blobs.removesuffix(".csv")}_elemental_composition.svg', bbox_inches='tight', format='svg')

    plt.clf()

    group_name_yields.plot(kind='bar', stacked=True, rot=0, color=colors)
    plt.legend(title='Group', title_fontsize='x-large')
    plt.title('Grouped carbon number distribution', fontsize='xx-large', fontweight='bold')
    plt.xlabel('C#', fontsize='x-large')
    plt.ylabel('Yield [wt. %]', fontsize='x-large')
    plt.xticks(np.arange(1, max(group_name_yields.index) + 1, 2))
    plt.savefig(f'{path_blobs.removesuffix(".csv")}_grouped_cnumber.svg', bbox_inches='tight', format='svg')

    plt.clf()

    plt.bar(grouped.index, grouped['Yield [wt.%]'], label='Yield [wt. %]', color=colors[6])
    plt.legend(title_fontsize='x-large')
    plt.title('Carbon number distribution', fontsize='xx-large', fontweight='bold')
    plt.xlabel('C#', fontsize='x-large')
    plt.ylabel('Yield [wt. %]', fontsize='x-large')
    plt.xticks(np.arange(1, max(grouped.index) + 1, 2))
    plt.savefig(f'{path_blobs.removesuffix(".csv")}_cnumber.svg', bbox_inches='tight', format='svg')

    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(7.5, 6))
    bubble_chart = ax.scatter(x.flatten(), y.flatten(), s=bubble_size * 50, c=bubble_color, alpha=0.8,
                              edgecolors='w', cmap='viridis')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_xticks(range(len(x_labels)))
    ax.set_yticks(range(len(y_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_yticklabels(y_labels)
    ax.set_xlabel('Group Name')
    cbar = plt.colorbar(bubble_chart, ax=ax)
    cbar.set_label('Yield [wt.%]')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_title('PIONA', fontsize='xx-large', fontweight='bold')
    plt.savefig(f'{path_blobs.removesuffix(".csv")}_piona.svg', bbox_inches='tight', format='svg')
    plt.clf()

    now = datetime.now()
    date_time_str = now.strftime("%Y-%m-%d_%H-%M-%S")
    blob_data_to_export = blob_df.loc[:, ['Compound Name', 'C#', 'Group Name', 'Volume', 'Yield [wt.%]']]
    group_name_yields.fillna(0, inplace=True)
    grouped.fillna(0, inplace=True)
    grouped = grouped.loc[:, 'Yield [wt.%]']
    non_zero_elements.name = 'Share [wt.%]'
    overview = pd.DataFrame({'Date': [now.date().__str__()] + (len(calibrants) - 1) * [''],
                             'Calibrant(s)': [calibrant.compound.name for calibrant in calibrants.values()],
                             'Calibration type': [calibrant.cal_type for calibrant in calibrants.values()],
                             'Calibration curve': [calibrant.curve for calibrant in calibrants.values()],
                             'Normalized': [normalized] + (len(calibrants) - 1) * [''],
                             'Mass closure': [mass_closure] + (len(calibrants) - 1) * [''], })
    save_path: str = f"{path_blobs.removesuffix(".csv")}_{date_time_str}_output.xlsx"
    with pd.ExcelWriter(save_path, engine='xlsxwriter') as writer:
        overview.to_excel(writer, sheet_name='Overview', index=False, startrow=0, startcol=0)
        blob_data_to_export.to_excel(writer, sheet_name='Compound yields', index=False, startrow=0, startcol=0)
        group_name_yields.to_excel(writer, sheet_name='Grouped Yields', index=True, startrow=0, startcol=0)
        grouped.to_excel(writer, sheet_name='C# Distribution', index=True, startrow=0, startcol=0)
        non_zero_elements.to_excel(writer, sheet_name='Elemental Composition', index=True, startrow=0, startcol=0)
        piona.to_excel(writer, sheet_name='PIONA', index=True, startrow=0, startcol=0)
    # Open Excel and the file
    excel = win32.Dispatch("Excel.Application")
    excel.Visible = True  # Ensure it opens in a new window
    workbook = excel.Workbooks.Open(save_path)


if __name__ == '__main__':
    logging.info('Main program started')
    main()
    logging.info('Main program finished')

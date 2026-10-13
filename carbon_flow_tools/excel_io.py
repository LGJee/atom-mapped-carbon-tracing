import shutil
from copy import copy
from pathlib import Path

import pandas as pd


def format_equation_column(worksheet):
    """Limit equation width and wrap long equations."""
    for header in worksheet[1]:
        if header.value != "Equation":
            continue
        dimension = worksheet.column_dimensions[header.column_letter]
        dimension.width = 50
        dimension.bestFit = False
        for row in worksheet.iter_rows(min_col=header.column, max_col=header.column):
            cell = row[0]
            alignment = copy(cell.alignment)
            alignment.wrap_text = True
            cell.alignment = alignment
        if worksheet.protection.sheet:
            worksheet.protection.formatColumns = False
        break


def clear_formatting(workbook):
    """Keep plain output tables with wrapped equations."""
    for worksheet in workbook:
        for row in worksheet:
            for cell in row:
                cell._style = None
        for dimension in worksheet.row_dimensions.values():
            dimension.height = None
            dimension._style = None
        for dimension in worksheet.column_dimensions.values():
            dimension._style = None
        worksheet.sheet_format.defaultRowHeight = 15
        worksheet.freeze_panes = None
        worksheet.auto_filter.ref = None
        worksheet.sheet_view.showGridLines = True
        worksheet.sheet_properties.tabColor = None
        for rule in list(worksheet.conditional_formatting):
            del worksheet.conditional_formatting[rule]
        format_columns(worksheet)


def copy_input(source: Path, destination: Path) -> Path:
    """Copy input before analysis; never write results into the input."""
    if source.resolve() == destination.resolve():
        raise ValueError("Output cannot overwrite the MATLAB input.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(destination))
    print(f"Input copied; analysis output: {destination}", flush=True)
    return destination


def write_sheet(workbook, name: str, table: pd.DataFrame) -> None:
    if name in workbook.sheetnames:
        del workbook[name]
    worksheet = workbook.create_sheet(name)
    worksheet.sheet_state = "hidden" if name == "Mapping_Info" else "visible"
    worksheet.append(list(table.columns))
    for row in table.itertuples(index=False, name=None):
        worksheet.append([None if pd.isna(value) else value for value in row])
    format_columns(worksheet)


def format_columns(worksheet):
    """Reuse readable plain column widths for copied inputs and generated tables."""
    widths = {
        "ReactionID": 18,
        "MetaboliteID": 18,
        "MetaboliteName": 38,
        "SMILES_Yeast9": 65,
        "SMILES_Manual": 65,
        "ReactionName": 42,
        "Flux": 14,
        "Equation": 50,
        "ReactionSMILES": 65,
        "MappedReactionSMILES": 65,
        "Step": 10,
        "Source": 38,
        "SourceIDs": 20,
        "Production": 38,
        "ProductionID": 15,
        "CarbonLabels": 30,
        "LabeledSMILES": 70,
        "SelectedMapper": 18,
        "ReviewRequired": 18,
        "SelectionReason": 60,
        "SupportingMappers": 35,
        "ValidMappers": 35,
        "MappingTemplate": 55,
        "SelectedMappedReactionSMILES": 65,
        "MappingStatus": 40,
        "MappingDirection": 18,
        "HeavyAtomAgreement": 22,
    }
    for header in worksheet[1]:
        column, letter = str(header.value), header.column_letter
        width = (
            44
            if column.endswith(("_ScoreType", "_Validation"))
            else widths.get(column, 18)
        )
        if "_LabeledCarbon" in column or column == "TotalCarbonMetabolites":
            width = min(80, max(width, len(column) + 2))
        worksheet.column_dimensions[letter].width = width
    format_equation_column(worksheet)

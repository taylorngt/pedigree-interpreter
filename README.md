# Pedigree Interpreter Package

PedInterp is a python package used to interpret images of clinical pedigrees
into PED file formatting.


## Prerequisites

This package depends on the [Tesseract OCR Engine](https://github.com/tesseract-ocr/tesseract), which must be installed separately from the Python dependencies listed in 'pyproject.toml'.

Installation instruction for Tesseract OCR engin are available in the official Teserract documentation:

[Tesseract Installation Guide](https://tesseract-ocr.github.io/tessdoc/Installation.html)


## Installation

All required Python dependencies will be installed automatically with this package using pip.

```bash
pip install PedInterp
```

## Usage
This package has two main functions:
1. **`pedigree_processing(FamID, raw_image_dir, export, output_dir)`** - translates the clinical pedigree image (.png) with the name given as FamID as the file name into a Pandas dataframe with field corresponding to PED file formatting. If desired, the interpretted pedigree can automatically exported as a PED file to the desired output directory (output_dir)
2. **`PED_export(PED_df, output_dir)`** - exports a given PED file-formatted Pandas dataframe (PED_df, constructed via pedigree_processing) to a desired output directory. File will be named according to the Family ID listed in first individual entry in PED_df (i.e. output_dir/<FamilyID>.ped)

## License

[MIT]
(https://choosealicense.com/licenses/mit/)
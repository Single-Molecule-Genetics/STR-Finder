# STR Finder for Genome Sequences

## Overview
This Python script is designed to **detect Short Tandem Repeats (STRs)** in a given genome FASTA file. STRs are repetitive DNA sequences that appear in various regions of the genome and can be crucial for genetic studies, forensic analysis, and evolutionary biology.

## Features
- Identifies **mononucleotide, dinucleotide, and trinucleotide repeats**.
- Uses **regular expressions (regex)** to efficiently search for STR patterns.
- Processes genome sequences from **FASTA files**.
- Saves the results as a **pickle file**, which can be loaded for further analysis.
- Uses **command-line arguments** to allow flexibility in defining repeat detection thresholds.

## How It Works
### 1. Input:
- The script requires a **FASTA** file containing the genome sequence.
- The user specifies **minimum repeat lengths** for mono-, di-, and trinucleotide STRs.

### 2. Processing:
- The genome is read and parsed using **pyfasta**.
- **Regex patterns** are used to identify STRs based on defined thresholds.
- STRs are stored in dictionaries for easy access and further analysis.

### 3. Output:
- The script saves a `.pickle` file containing:
  - `ref_str`: STR coordinates for each chromosome.
  - `ref_info`: The genome sequence.
  - `chromosome_STRdict`: Annotated STRs mapped to their chromosome positions.
  - `str_region_annotation`: STR metadata for further interpretation.

## Installation
Before running the script, install dependencies:
```bash
pip install pyfasta
```

## Usage
Run the script using the command line:
```bash
python script.py -g /path/to/genome.fa -o /path/to/output/ --min_mono 5 --min_di 6 --min_tri 4
```

### Command-Line Arguments:
- `-g, --genome` → Path to the genome FASTA file.
- `-o, --output` → Directory where results will be saved.
- `--min_mono` → Minimum mononucleotide STR length (default: 5).
- `--min_di` → Minimum dinucleotide STR length (default: 6).
- `--min_tri` → Minimum trinucleotide STR length (default: 4).

## Example Run
```bash
python script.py -g data/genome.fa -o results/ --min_mono 6 --min_di 8 --min_tri 5
```

## Output File Format
The script outputs a **pickle file** named:
```
CuteCV_STRs_{min_mono}m_{min_di}d_{min_tri}t.pickle
```
which contains structured data for STR analysis.

## Why Use This Script?
- **Efficient:** Uses optimized regex patterns to detect STRs quickly.
- **Flexible:** Users can set custom STR length thresholds.
- **Reusable:** Outputs are structured for further processing in bioinformatics pipelines.

## Future Improvements
- Support for additional STR lengths (tetra- and pentanucleotide repeats).
- Integration with **bed files** for genome annotation.
- Output in **CSV or JSON** for broader usability.

## License
This project is licensed under the **MIT License**. Feel free to contribute, modify, or use it for your research.

## Contributions
If you'd like to contribute:
1. Fork the repository.
2. Create a new branch for your feature.
3. Submit a pull request with a description of your changes.

---
**Developed for Bioinformatics Research and STR Analysis!**

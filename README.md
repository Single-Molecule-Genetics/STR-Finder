# STR Finder for Genome Sequences

## Overview
This Python script is specifically designed to identify Short Tandem Repeats (STRs) within a given genome FASTA file. 
STRs are repetitive DNA sequences that occur throughout the genome and play a crucial role in genetic studies, 
forensic analysis, and evolutionary biology. At present, the script supports only perfect repeats.

The detection process employs regular expressions (regex) to identify continuous repeat stretches composed of 
mononucleotide (1-bp), dinucleotide (2-bp), and trinucleotide (3-bp) motifs. By default, STRs are defined based on
the following minimum repeat thresholds, which can be adjusted as needed:
* 	Five consecutive repeats for mononucleotides 
*   Three consecutive repeats for dinucleotides and trinucleotides

Once detected, STR-Finder leverages regex-based pattern matching to efficiently extract key STR features, including:
* Repeat motif (sequence of the repeating unit)
* Start and end positions in the genome 
* Repeat length, calculated based on the motif size and genomic coordinates

Finally, the identified STRs are annotated with chromosome number, genomic coordinates, and repeat length. 
All results are stored in structured Python dictionaries, ensuring efficient data handling for downstream analysis.

This tool provides user-configurable threshold settings, allowing flexibility to tailor the detection criteria based 
on specific research requirements.
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
Genome_{min_mono}m_{min_di}d_{min_tri}t.pickle
```
The output pickle file contains a Python dictionary with structured STR data. Each key represents a different method 
of organizing the STR coordinates:
	•	ref_str: A dictionary where chromosomes serve as keys, and the associated values are lists of STR coordinate ranges present on each chromosome.
	•	ref_info: A dictionary with chromosomes as keys and their corresponding FASTA-formatted sequences as values.
	•	chromosome_STRdict: A dictionary with chromosomes as keys, each containing a nested dictionary where the start position of an STR is the key, and the corresponding annotation—including the end position—is the value. For example: I-430-435-6A.
	•	str_region_annotation: Similar to chromosome_STRdict, but instead of using only the start position as a key, it utilizes ranges (start-end) in the nested dictionary.

This structured format allows for easy retrieval, processing, and further analysis of STR data.


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
## About Us  
This tool was developed by **Shehab Moukbel ALi Aldawla** as part of the **Single Molecule Genetics (SMG) team** at JKU.  
For more information, visit our [team page](https://www.jku.at/institut-fuer-biophysik/ueber-uns/team/single-molecule-genetics-irene-tiemann-boege/) or join our community discussions.
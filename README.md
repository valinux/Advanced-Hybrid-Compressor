# Advanced Hybrid Compressor

This repository contains an experimental Python-based compressor that combines LZ77-style dictionary tokenization with high-precision arithmetic coding. The design is inspired by seminal academic work in data compression, including:

- **Witten, I. H., Neal, R. M., & Cleary, J. G. (1987).** [Arithmetic Coding for Data Compression](https://dl.acm.org/doi/10.1145/318271.318279)
- **Duda, J. (2009).** [Asymmetric Numeral Systems: Entropy Coding Combining Speed of Huffman Coding with Compression Rate of Arithmetic Coding](https://arxiv.org/abs/1311.2540)

## Overview

This compressor uses a two-step process:

1. **Tokenization (LZ77-style):**  
   The input file is scanned using a sliding window. Data is divided into tokens that represent:
   - **Literal Tokens:** Single bytes from the input.
   - **Match Tokens:** Repeated sequences represented by an offset (how far back the match starts) and a match length.

2. **Arithmetic Coding:**  
   The token stream is then encoded using a high-precision arithmetic coder. Frequency tables for the token types (literal vs. match), literal byte values, match offsets, and match lengths are built adaptively and stored in a separate JSON mapping file. This mapping file is later used during decompression to rebuild the probability models accurately.

## Files

- **advanced_compressor.py**: The main Python script containing the compressor and decompressor.
- **mapping.json**: (Generated during compression) Contains frequency tables and compression parameters required for decompression.
- **README.md**: This file.

## Requirements

- Python 3.x (No third-party libraries are required.)

## Usage

### Compressing a File

To compress a file (e.g., `1.txt`), run the following command:

```bash
python advanced_compressor.py compress 1.txt mapping.json compressed.bin
```

- **1.txt**: Input file to compress.
- **mapping.json**: Output mapping file that stores the frequency models and parameters.
- **compressed.bin**: Output file containing the compressed data.

### Decompressing a File

To decompress the file, run:

```bash
python advanced_compressor.py decompress compressed.bin mapping.json decompressed.txt
```

- **compressed.bin**: The compressed binary file.
- **mapping.json**: The mapping file generated during compression.
- **decompressed.txt**: The output file, which should match the original `1.txt`.

## Code Overview

- **Tokenization**:  
  The `tokenize_data(data: bytes)` function implements LZ77-style tokenization using a sliding window, producing tokens for literals and matches.

- **Frequency Table Construction**:  
  The `build_frequency_tables(tokens)` function builds frequency models for token types, literal values, match offsets, and match lengths. These tables are later converted to cumulative frequency tables for arithmetic coding.

- **Arithmetic Coding**:  
  The `ArithmeticEncoder` and `ArithmeticDecoder` classes implement high-precision arithmetic coding (inspired by Witten et al. and Duda) to compress the token stream.

- **Mapping File Handling**:  
  Frequency tables and other compression parameters are stored in a JSON file. During decompression, a helper function converts the JSON keys (stored as strings) back to integers.

- **Command-Line Interface**:  
  The script supports two commands: `compress` and `decompress`. The corresponding functions process the file based on the provided arguments.

## Limitations and Future Improvements

This project is a research-inspired prototype. Potential improvements include:

- **Optimization:**  
  Enhancing the token search algorithm to improve performance.

- **Adaptive Models:**  
  Incorporating context mixing or on-the-fly adaptation of the probability models.

- **Alternative Coders:**  
  Experimenting with Asymmetric Numeral Systems (ANS) as an alternative to arithmetic coding.

- **Robustness:**  
  Improving error handling, especially for large files, and exploring multi-threading for performance gains.

## References

- Witten, I. H., Neal, R. M., & Cleary, J. G. (1987). [Arithmetic Coding for Data Compression](https://dl.acm.org/doi/10.1145/318271.318279)
- Duda, J. (2009). [Asymmetric Numeral Systems: Entropy Coding Combining Speed of Huffman Coding with Compression Rate of Arithmetic Coding](https://arxiv.org/abs/1311.2540)

## License

This project is licensed under the MIT License.

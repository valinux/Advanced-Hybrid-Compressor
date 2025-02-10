#!/usr/bin/env python3
"""
Advanced Hybrid Compressor
References:
  - Witten, I. H., Neal, R. M., & Cleary, J. G. (1987). Arithmetic coding for data compression.
  - Duda, J. (2009). Asymmetric numeral systems: entropy coding combining speed of Huffman coding with compression rate of arithmetic coding.

This compressor combines a dictionary‐based (LZ77-style) parsing with adaptive arithmetic coding.
It tokenizes the input into literal and match tokens, builds frequency models from a first pass,
and then encodes the token stream using high-precision arithmetic coding.
Mapping (frequency tables and parameters) is stored in a separate JSON file.
This design is meant to be a robust experimental framework—there is plenty of room to further
incorporate ideas like context mixing or ANS.
"""

import argparse
import json
import math
import os
import struct
import sys

# ---------------------------
# Global Parameters
# ---------------------------
WINDOW_SIZE = 1024         # Sliding window size for LZ77
MIN_MATCH_LENGTH = 3       # Minimum match length to consider a match token
MAX_MATCH_LENGTH = 255     # Maximum match length (to allow one-byte representation of length offset)

# ---------------------------
# Arithmetic Coding Classes
# ---------------------------
PRECISION = 32
MAX_RANGE = (1 << PRECISION) - 1
HALF = 1 << (PRECISION - 1)
QUARTER = 1 << (PRECISION - 2)
THREE_QUARTER = 3 * QUARTER

class ArithmeticEncoder:
    def __init__(self, precision=PRECISION):
        self.precision = precision
        self.low = 0
        self.high = MAX_RANGE
        self.underflow = 0
        self.output_bits = []

    def encode_symbol(self, symbol, cum_freq, total):
        """
        Encode a symbol given its cumulative frequency table (cum_freq) and total count.
        'cum_freq' should be a list where cum_freq[i] is the cumulative count up to symbol i.
        """
        range_width = self.high - self.low + 1
        self.high = self.low + (range_width * cum_freq[symbol + 1] // total) - 1
        self.low = self.low + (range_width * cum_freq[symbol] // total)

        while True:
            if self.high < HALF:
                self._output_bit(0)
            elif self.low >= HALF:
                self._output_bit(1)
                self.low -= HALF
                self.high -= HALF
            elif self.low >= QUARTER and self.high < THREE_QUARTER:
                self.underflow += 1
                self.low -= QUARTER
                self.high -= QUARTER
            else:
                break
            self.low *= 2
            self.high = self.high * 2 + 1

    def _output_bit(self, bit):
        self.output_bits.append(bit)
        for _ in range(self.underflow):
            self.output_bits.append(1 - bit)
        self.underflow = 0

    def finish(self):
        self.underflow += 1
        if self.low < QUARTER:
            self._output_bit(0)
        else:
            self._output_bit(1)
        return self.output_bits

class ArithmeticDecoder:
    def __init__(self, bit_stream, precision=PRECISION):
        self.precision = precision
        self.low = 0
        self.high = MAX_RANGE
        self.code = 0
        self.bit_stream = bit_stream
        self.bit_index = 0
        for _ in range(precision):
            self.code = (self.code << 1) | self._read_bit()

    def _read_bit(self):
        if self.bit_index < len(self.bit_stream):
            bit = self.bit_stream[self.bit_index]
            self.bit_index += 1
            return bit
        return 0

    def decode_symbol(self, cum_freq, total):
        range_width = self.high - self.low + 1
        value = ((self.code - self.low + 1) * total - 1) // range_width

        # Find the symbol s such that cum_freq[s] <= value < cum_freq[s+1]
        symbol = 0
        while cum_freq[symbol+1] <= value:
            symbol += 1

        self.high = self.low + (range_width * cum_freq[symbol+1] // total) - 1
        self.low = self.low + (range_width * cum_freq[symbol] // total)

        while True:
            if self.high < HALF:
                pass
            elif self.low >= HALF:
                self.code -= HALF
                self.low -= HALF
                self.high -= HALF
            elif self.low >= QUARTER and self.high < THREE_QUARTER:
                self.code -= QUARTER
                self.low -= QUARTER
                self.high -= QUARTER
            else:
                break
            self.low *= 2
            self.high = self.high * 2 + 1
            self.code = (self.code << 1) | self._read_bit()
        return symbol

# ---------------------------
# Bit Utilities
# ---------------------------
def bits_to_bytes(bit_list):
    total_bits = len(bit_list)
    num_bytes = math.ceil(total_bits / 8)
    padded_length = num_bytes * 8
    padded = bit_list + [0] * (padded_length - total_bits)
    b = bytearray()
    for i in range(0, padded_length, 8):
        byte = 0
        for bit in padded[i:i+8]:
            byte = (byte << 1) | bit
        b.append(byte)
    valid_bits_last = total_bits % 8 if total_bits % 8 != 0 else 8
    return bytes(b), valid_bits_last

def bytes_to_bit_list(b, valid_bits_last):
    bit_list = []
    for byte in b[:-1]:
        for i in range(7, -1, -1):
            bit_list.append((byte >> i) & 1)
    for i in range(7, 7 - valid_bits_last, -1):
        bit_list.append((b[-1] >> i) & 1)
    return bit_list

# ---------------------------
# Helper Function to Convert Frequency Dictionary Keys
# ---------------------------
def convert_freq_dict(d):
    """
    Convert a dictionary with string keys (from JSON) to integer keys.
    """
    return {int(k): v for k, v in d.items()}

# ---------------------------
# LZ77 Tokenization
# ---------------------------
# Tokens are tuples:
#   ('L', literal_value) for literal tokens
#   ('M', offset, length) for match tokens
def tokenize_data(data: bytes):
    tokens = []
    i = 0
    n = len(data)
    while i < n:
        best_length = 0
        best_offset = 0
        window_start = max(0, i - WINDOW_SIZE)
        # Search for the longest match in the sliding window.
        for j in range(window_start, i):
            length = 0
            while (i + length < n and j + length < n and
                   data[j + length] == data[i + length] and
                   length < MAX_MATCH_LENGTH):
                length += 1
            if length >= MIN_MATCH_LENGTH and length > best_length:
                best_length = length
                best_offset = i - j  # offset: how far back the match starts
        if best_length >= MIN_MATCH_LENGTH:
            tokens.append(('M', best_offset, best_length))
            i += best_length
        else:
            tokens.append(('L', data[i]))
            i += 1
    return tokens

# ---------------------------
# Frequency Table Construction
# ---------------------------
def build_frequency_tables(tokens):
    # For token type: use 0 for literal and 1 for match.
    token_type_freq = {0: 1, 1: 1}  # Laplace smoothing
    literal_freq = {i: 1 for i in range(256)}
    offset_freq = {i: 1 for i in range(WINDOW_SIZE)}  # We encode offset as (offset - 1), so valid symbols: 0..WINDOW_SIZE-1
    length_range = MAX_MATCH_LENGTH - MIN_MATCH_LENGTH + 1
    length_freq = {i: 1 for i in range(length_range)}  # encode length as (length - MIN_MATCH_LENGTH)
    
    for token in tokens:
        if token[0] == 'L':
            token_type_freq[0] += 1
            literal_freq[token[1]] += 1
        elif token[0] == 'M':
            token_type_freq[1] += 1
            offset_symbol = token[1] - 1
            length_symbol = token[2] - MIN_MATCH_LENGTH
            offset_freq[offset_symbol] += 1
            length_freq[length_symbol] += 1
    return token_type_freq, literal_freq, offset_freq, length_freq

def build_cumulative(freq_dict, max_symbol):
    cum = [0]
    total = 0
    for symbol in range(max_symbol + 1):
        total += freq_dict[symbol]
        cum.append(total)
    return cum, total

# ---------------------------
# Compression: Encoding the Token Stream
# ---------------------------
def compress_tokens(tokens, freq_tables):
    # Unpack frequency tables
    token_type_freq, literal_freq, offset_freq, length_freq = freq_tables
    # Build cumulative frequency tables for each alphabet.
    cum_tt, total_tt = build_cumulative(token_type_freq, 1)  # token type: 0 or 1
    cum_lit, total_lit = build_cumulative(literal_freq, 255)
    cum_off, total_off = build_cumulative(offset_freq, WINDOW_SIZE - 1)
    length_range = MAX_MATCH_LENGTH - MIN_MATCH_LENGTH
    cum_len, total_len = build_cumulative(length_freq, length_range)
    
    encoder = ArithmeticEncoder(precision=PRECISION)
    # For each token, encode its fields.
    for token in tokens:
        if token[0] == 'L':
            # Encode token type = 0
            encoder.encode_symbol(0, cum_tt, total_tt)
            # Encode literal value (0..255)
            encoder.encode_symbol(token[1], cum_lit, total_lit)
        elif token[0] == 'M':
            # Encode token type = 1
            encoder.encode_symbol(1, cum_tt, total_tt)
            # Encode offset (store as offset - 1)
            encoder.encode_symbol(token[1] - 1, cum_off, total_off)
            # Encode length (store as length - MIN_MATCH_LENGTH)
            encoder.encode_symbol(token[2] - MIN_MATCH_LENGTH, cum_len, total_len)
    bits = encoder.finish()
    return bits

# ---------------------------
# Compressor: Putting It All Together
# ---------------------------
def compress_file(input_filename: str, mapping_filename: str, output_filename: str):
    try:
        with open(input_filename, "rb") as infile:
            data = infile.read()
    except Exception as e:
        print(f"Error reading input file '{input_filename}': {e}", file=sys.stderr)
        return

    if not data:
        print("Input file is empty. Nothing to compress.", file=sys.stderr)
        return

    # Tokenize the data using our LZ77-style parser.
    tokens = tokenize_data(data)
    # Build frequency tables (with Laplace smoothing).
    freq_tables = build_frequency_tables(tokens)
    
    # Encode the token stream using arithmetic coding.
    bit_list = compress_tokens(tokens, freq_tables)
    encoded_bytes, valid_bits_last = bits_to_bytes(bit_list)

    # Write the compressed file. Prepend a one-byte header indicating how many bits in the last byte are valid.
    try:
        with open(output_filename, "wb") as outfile:
            outfile.write(struct.pack("B", valid_bits_last))
            outfile.write(encoded_bytes)
    except Exception as e:
        print(f"Error writing compressed file '{output_filename}': {e}", file=sys.stderr)
        return

    # Save mapping/model information to a JSON mapping file.
    mapping_data = {
        "window_size": WINDOW_SIZE,
        "min_match_length": MIN_MATCH_LENGTH,
        "max_match_length": MAX_MATCH_LENGTH,
        "original_length": len(data),
        "token_type_freq": freq_tables[0],
        "literal_freq": freq_tables[1],
        "offset_freq": freq_tables[2],
        "length_freq": freq_tables[3]
    }
    try:
        with open(mapping_filename, "w") as map_file:
            json.dump(mapping_data, map_file, indent=2)
    except Exception as e:
        print(f"Error writing mapping file '{mapping_filename}': {e}", file=sys.stderr)
        return

    print(f"Compression successful:\n  Input size: {len(data)} bytes\n  Compressed size: {os.path.getsize(output_filename)} bytes\n  Mapping file: {mapping_filename}")

# ---------------------------
# Decompression: Decoding the Token Stream and Reconstructing Data
# ---------------------------
def decompress_file(compressed_filename: str, mapping_filename: str, output_filename: str):
    # Load mapping/model file.
    try:
        with open(mapping_filename, "r") as map_file:
            mapping_data = json.load(map_file)
    except Exception as e:
        print(f"Error reading mapping file '{mapping_filename}': {e}", file=sys.stderr)
        return

    # Extract parameters.
    window_size = mapping_data["window_size"]
    min_match_length = mapping_data["min_match_length"]
    max_match_length = mapping_data["max_match_length"]
    original_length = mapping_data["original_length"]
    
    # Convert frequency dictionaries keys from strings back to integers.
    token_type_freq = convert_freq_dict(mapping_data["token_type_freq"])
    literal_freq = convert_freq_dict(mapping_data["literal_freq"])
    offset_freq = convert_freq_dict(mapping_data["offset_freq"])
    length_freq = convert_freq_dict(mapping_data["length_freq"])

    # Rebuild cumulative frequency tables.
    cum_tt, total_tt = build_cumulative(token_type_freq, 1)
    cum_lit, total_lit = build_cumulative(literal_freq, 255)
    cum_off, total_off = build_cumulative(offset_freq, window_size - 1)
    length_range = max_match_length - min_match_length
    cum_len, total_len = build_cumulative(length_freq, length_range)

    # Read compressed file.
    try:
        with open(compressed_filename, "rb") as comp_file:
            header = comp_file.read(1)
            if len(header) != 1:
                print("Compressed file is missing header.", file=sys.stderr)
                return
            (valid_bits_last,) = struct.unpack("B", header)
            compressed_bytes = comp_file.read()
    except Exception as e:
        print(f"Error reading compressed file '{compressed_filename}': {e}", file=sys.stderr)
        return

    if not compressed_bytes:
        print("Compressed file contains no data after header.", file=sys.stderr)
        return

    bit_list = bytes_to_bit_list(compressed_bytes, valid_bits_last)
    decoder = ArithmeticDecoder(bit_list, precision=PRECISION)

    # Decode tokens until we have reconstructed original_length bytes.
    output_data = bytearray()
    while len(output_data) < original_length:
        # Decode token type.
        token_type = decoder.decode_symbol(cum_tt, total_tt)
        if token_type == 0:
            # Literal token: decode literal value.
            literal_val = decoder.decode_symbol(cum_lit, total_lit)
            output_data.append(literal_val)
        elif token_type == 1:
            # Match token: decode offset and length.
            offset_symbol = decoder.decode_symbol(cum_off, total_off)
            length_symbol = decoder.decode_symbol(cum_len, total_len)
            offset = offset_symbol + 1
            length = length_symbol + min_match_length
            start = len(output_data) - offset
            if start < 0:
                print("Error during decompression: invalid match offset.", file=sys.stderr)
                return
            for i in range(length):
                output_data.append(output_data[start + i])
        else:
            print("Error: unknown token type encountered.", file=sys.stderr)
            return

    try:
        with open(output_filename, "wb") as outfile:
            outfile.write(output_data)
    except Exception as e:
        print(f"Error writing decompressed file '{output_filename}': {e}", file=sys.stderr)
        return

    print(f"Decompression successful. Output file size: {len(output_data)} bytes")

# ---------------------------
# Command-Line Interface
# ---------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Advanced Hybrid Compressor\n"
                    "This tool compresses/decompresses files using a novel combination of LZ77-style tokenization\n"
                    "and arithmetic coding. Mapping (frequency models and parameters) are stored separately.\n"
                    "References:\n"
                    "  - Witten, I. H., Neal, R. M., & Cleary, J. G. (1987)\n"
                    "  - Duda, J. (2009)"
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Commands: compress, decompress")

    parser_compress = subparsers.add_parser("compress", help="Compress a file")
    parser_compress.add_argument("input", help="Input file to compress")
    parser_compress.add_argument("mapping", help="Output mapping file (JSON)")
    parser_compress.add_argument("output", help="Output compressed file")

    parser_decompress = subparsers.add_parser("decompress", help="Decompress a file")
    parser_decompress.add_argument("input", help="Input compressed file")
    parser_decompress.add_argument("mapping", help="Mapping file (JSON) used during compression")
    parser_decompress.add_argument("output", help="Output decompressed file")

    args = parser.parse_args()

    if args.command == "compress":
        compress_file(args.input, args.mapping, args.output)
    elif args.command == "decompress":
        decompress_file(args.input, args.mapping, args.output)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

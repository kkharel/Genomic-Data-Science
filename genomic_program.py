#Write a Python program that takes as input a file containing DNA sequences in multi-FASTA format,
#and computes the answers to the following questions. 

#While developing program, please use the following example file to test  work: dna.example.fasta

#(1) How many records are in the file? 
#A record in a FASTA file is defined as a single-line header, followed by lines of sequence data. 
#The header line is distinguished from the sequence data by a greater-than (">") symbol in the 
#first column. The word following the ">" symbol is the identifier of the sequence, 
#and the rest of the line is an optional description of the entry. 
#There should be no space between the ">" and the first letter of the identifier. 

#(2) What are the lengths of the sequences in the file? 
#What is the longest sequence and what is the shortest sequence? 
#Is there more than one longest or shortest sequence? What are their identifiers? 

#(3) In molecular biology, a reading frame is a way of dividing the DNA sequence of nucleotides 
#into a set of consecutive, non-overlapping triplets (or codons). Depending on where we start, 
#there are six possible reading frames: three in the forward (5' to 3') direction and 
#three in the reverse (3' to 5'). 
#For instance, the three possible forward reading frames for the sequence 
#AGGTGACACCGCAAGCCTTATATTAGC are: 

#AGG TGA CAC CGC AAG CCT TAT ATT AGC

#A GGT GAC ACC GCA AGC CTT ATA TTA GC

#AG GTG ACA CCG CAA GCC TTA TAT TAG C 

#These are called reading frames 1, 2, and 3 respectively. An open reading frame (ORF) is 
#the part of a reading frame that has the potential to encode a protein. 
#It starts with a start codon (ATG), and ends with a stop codon (TAA, TAG or TGA). 
#For instance, ATGAAATAG is an ORF of length 9.

#Given an input reading frame on the forward strand (1, 2, or 3), program should be able to identify 
#all ORFs present in each sequence of the FASTA file, and answer the following questions: 
#  what is the length of the longest ORF in the file? 
#  What is the identifier of the sequence containing the longest ORF? 
#  For a given sequence identifier, what is the longest ORF contained in the sequence represented by 
#  that identifier? 
#  What is the starting position of the longest ORF in the sequence that contains it? 
#  The position should indicate the character number in the sequence. 
#  For instance, the following ORF in reading frame 1:

#>sequence1

#ATGCCCTAG

#starts at position 1.

#Note that because the following sequence:

#>sequence2

#ATGAAAAAA

#does not have any stop codon in reading frame 1, we do not consider it to be an ORF in reading frame 1. 

#(4) A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) 
#somewhere in the sequence. Although repeats can occur on both the forward and reverse strands 
#of the DNA sequence, we will only consider repeats on the forward strand here. 
#Also we will allow repeats to overlap themselves. For example, the sequence ACACA contains 
#two copies of the sequence ACA - once at position 1 (index 0 in Python), and once at position 3. 
#Given a length n, your program should be able to identify all repeats of length n in all sequences 
#in the FASTA file. The program should also determine how many times each repeat occurs in the file, 
#and which is the most frequent repeat of a given length.


from Bio import SeqIO
from collections import defaultdict
import os

class SequenceAnalyzer:
  def __init__(self, filepath, filename):
    self.filepath = filepath
    self.filename = filename
    self.sequences = self.read_fasta()

  def read_fasta(self):
    full_path = f"{self.filepath}/{self.filename}"
    sequences = {}
    for seq_record in SeqIO.parse(full_path, "fasta"):
      sequences[seq_record.id] = str(seq_record.seq)
    return sequences

  def count_records(self):
    return len(self.sequences)

  def sequence_lengths(self):
    return {seq_id: len(seq) for seq_id, seq in self.sequences.items()}

  def longest_and_shortest_sequences(self):
    lengths = self.sequence_lengths()
    max_length = max(lengths.values())
    min_length = min(lengths.values())
    longest_ids = [seq_id for seq_id, length in lengths.items() if length == max_length]
    shortest_ids = [seq_id for seq_id, length in lengths.items() if length == min_length]
    return max_length, min_length, longest_ids, shortest_ids

  def find_orfs(self, frame):
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    results = {}
    for seq_id, sequence in self.sequences.items():
      orfs = []
      for i in range(frame - 1, len(sequence) - 2, 3):
        if sequence[i:i+3] == start_codon:
          for j in range(i+3, len(sequence) - 2, 3):
            if sequence[j:j+3] in stop_codons:
              orfs.append((i + 1, j + 3, j + 3 - i))
              break
      results[seq_id] = orfs
    return results

  def longest_orf(self, frame):
    max_length = 0
    max_orf = None
    max_seq_id = None

    for seq_id, orfs in self.find_orfs(frame).items():
      for start, end, length in orfs:
        if length > max_length:
          max_length = length
          max_orf = (start, end)
          max_seq_id = seq_id
    return max_length, max_seq_id, max_orf

  def longest_orf_in_sequence(self, seq_id, frame):
    orfs = self.find_orfs(frame).get(seq_id, [])
    if not orfs:
      return None, 0

    max_length = 0
    longest_start = 0

    for start, end, length in orfs:
      if length > max_length:
        max_length = length
        longest_start = start

    return longest_start, max_length

  def find_repeats(self, sequence, n):
    repeat_counts = defaultdict(int)
    for i in range(len(sequence)-n+1):
      repeat = sequence[i:i+n]
      repeat_counts[repeat] += 1
    return repeat_counts

  def get_count(self, item):
    return item[1]
      
  def most_frequent_repeat(self, n):
    total_repeats = defaultdict(int)

    for sequence in self.sequences.values():
      repeat_counts = self.find_repeats(sequence, n)
      for repeat, count in repeat_counts.items():
        total_repeats[repeat] += count
    most_frequent = max(total_repeats.items(), key=self.get_count)
    return most_frequent

if __name__ == "__main__":
  filepath = input("Enter the path to the directory containing the FASTA file: ").strip()
  filename = input("Enter the filename of the FASTA file: ").strip()
  
  while not (os.path.isdir(filepath) and os.path.isfile(os.path.join(filepath, filename))):
    print("Invalid path or filename. Please check and enter again.")
    filepath = input("Enter the path to the directory containing the FASTA file: ").strip()
    filename = input("Enter the filename of the FASTA file: ").strip()

  analyzer = SequenceAnalyzer(filepath, filename)
  print("Reading Sequences from file...")
  #for seq_id, sequence in analyzer.sequences.items():
  #  print(f"ID: {seq_id}, Sequence: {sequence}")
  
  num_records = analyzer.count_records()
  print(f"Number of records: {num_records}")
  
  lengths = analyzer.sequence_lengths()
  print("Lengths of sequences:", lengths)
  
  max_length, min_length, longest_ids, shortest_ids = analyzer.longest_and_shortest_sequences()
  print(f"Longest sequence length: {max_length}, IDs: {longest_ids}")
  print(f"Shortest sequence length: {min_length}, IDs: {shortest_ids}")
  
  frame = None
  while frame not in [1, 2, 3]:
    try:
      frame = int(input("Enter the reading frame (1, 2, or 3): "))
    except ValueError:
      print("Invalid input. Please enter a number.")
  
  longest_orf_length, longest_orf_seq_id, longest_orf_coords = analyzer.longest_orf(frame)
  print(f"Longest ORF length: {longest_orf_length}")
  print(f"Sequence ID with longest ORF: {longest_orf_seq_id}")
  if longest_orf_coords:
    print(f"Starting position of longest ORF: {longest_orf_coords[0]}")
  
  seq_id = input("Enter the sequence identifier to find the longest ORF: ").strip()
  while seq_id not in analyzer.sequences:
    print("Sequence ID not found. Please enter a valid sequence identifier.")
    seq_id = input("Enter the sequence identifier to find the longest ORF: ").strip()
  
  longest_start, longest_length = analyzer.longest_orf_in_sequence(seq_id, frame)
  if longest_start is not None:
    print(f"Longest ORF in sequence {seq_id} starts at position {longest_start} with length {longest_length}")
  
  repeat_length = None
  while repeat_length is None:
    try:
      repeat_length = int(input("Enter the repeat length to find the most frequent repeat: "))
      if repeat_length <= 0:
        raise ValueError("Repeat length must be a positive integer.")
    except ValueError as e:
      print(f"Invalid input: {e}")
      repeat_length = None
  
  most_frequent = analyzer.most_frequent_repeat(repeat_length)
  print(f"Most frequent repeat of length {repeat_length}: {most_frequent[0]} occurs {most_frequent[1]} times")
